import multiprocessing
import sys
from typing import Optional

import libcasm.configuration as casmconfig
import libcasm.enumerate as casmenum
import libcasm.xtal.prims as xtal_prims


def _worker(
    queue: multiprocessing.Queue,
    prim_data: dict,
    supercells_data: dict,
):
    """Enumerates all occupations and adds results to queue, for each supercell.

    To pass data between processes, parameters must be pickleable. Here we
    convert CASM objects to dicts which only contain plain Python data structures
    and are therefore automatically pickleable.

    Results are ConfigurationSet, which are converted to dict and put in `queue`. When
    the worker is done, it puts a sentinel value (None) in the queue to signal the end
    of data. If an exception occurs, the worker puts the exception message in the queue
    as a string.

    Parameters
    ----------
    queue : multiprocessing.Queue
        Queue to send results to. Queues are thread and process safe. Any object put
        into a multiprocessing queue will be serialized. In this example, each object
        is a dict representation of a ConfigurationSet.
    prim_data : dict
        Dict representation of a Prim.
    supercells_data : dict
        Dict representation of a SupercellSet.

    """

    try:
        # Read Prim and SupercellSet from the dict paramters
        prim = casmconfig.Prim.from_dict(data=prim_data)
        supercells = casmconfig.SupercellSet.from_dict(data=supercells_data, prim=prim)

        # Enumerate configurations for each supercell,
        # sending results to queue one supercell at a time
        for record in supercells:
            # Make enumerator
            supercell_set = casmconfig.SupercellSet(prim=prim)
            configuration_set = casmconfig.ConfigurationSet()
            config_enum = casmenum.ConfigEnumAllOccupations(
                prim=prim,
                supercell_set=supercell_set,
            )

            # Enumerate configs for this supercell
            for configuration in config_enum.by_supercell_list(
                supercells=[record.supercell]
            ):
                configuration_set.add(configuration)

            # Sends ConfigurationSet as a dict to queue
            print(
                f"Sending {len(configuration_set)} configurations "
                f"from worker {multiprocessing.current_process().name} "
                f"for supercell {record.supercell_name}..."
            )
            data = configuration_set.to_dict()
            queue.put(data)

        # Sentinel value to signal end of data from this worker
        queue.put(None)
    except Exception as e:
        print(f"Error in worker {multiprocessing.current_process().name}: {e}")
        queue.put(None)


def parallel_enum(
    prim: casmconfig.Prim,
    min_vol: int,
    max_vol: int,
    n_workers: int,
    supercells: Optional[casmconfig.SupercellSet] = None,
    configurations: Optional[casmconfig.ConfigurationSet] = None,
):
    """Enumerate all configurations for a given Prim and supercell volume range,
    using multiple processes to parallelize the work by creating one task to enumerate
    configurations for each supercell

    Parameters
    ----------
    prim : casmconfig.Prim
        The prim.
    min_vol : int
        Minimum supercell volume.
    max_vol : int
        Maximum supercell volume.
    n_workers : int
        Number of worker processes to use.
    supercells : Optional[casmconfig.SupercellSet] = None
        A SupercellSet to collect enumerated supercells in. If None, a new SupercellSet
        will be created.
    configurations : Optional[casmconfig.ConfigurationSet] = None
        A ConfigurationSet to collect enumerated configurations in. If None, a new
        ConfigurationSet will be created.

    Returns
    -------
    supercells: casmconfig.SupercellSet
        The SupercellSet containing all enumerated supercells.
    configurations: casmconfig.ConfigurationSet
        The ConfigurationSet containing all enumerated configurations.
    """

    # Allocate supercells to workers
    allocated_supercells = [
        casmconfig.SupercellSet(prim=prim) for _ in range(n_workers)
    ]
    scel_enum = casmenum.ScelEnum(
        prim=prim,
    )
    for i, supercell in enumerate(scel_enum.by_volume(min=min_vol, max=max_vol)):
        allocated_supercells[i % n_workers].add_supercell(supercell)

    # Convert to plain-old-data for multiprocessing
    prim_data = prim.to_dict()
    allocated_supercells = [s.to_dict() for s in allocated_supercells]

    # Queue for collecting generated ConfigurationSets as dict
    queue = multiprocessing.Queue()

    # Generate configurations - in child processes
    processes = []
    for i in range(n_workers):
        process = multiprocessing.Process(
            target=_worker,
            args=(queue, prim_data, allocated_supercells[i]),
        )
        print("Starting worker", i)
        process.start()
        processes.append(process)

    # Collect generated configurations - in parent process
    if supercells is None:
        supercells = casmconfig.SupercellSet(prim=prim)
    if configurations is None:
        configurations = casmconfig.ConfigurationSet()
    print("Waiting for configurations...")
    sys.stdout.flush()
    n_workers_completed = 0
    while True:
        # Get data from the queue
        data = queue.get()

        # Sentinel value to signal end of data
        if data is None:
            n_workers_completed += 1
            if n_workers_completed == n_workers:
                break
            else:
                continue

        # Check for errors
        if not isinstance(data, dict):
            print(data)
            raise ValueError(f"Expected dict, got {type(data)}")

        # Convert data to ConfigurationSet - saving supercells in `supercells`
        dataset = casmconfig.ConfigurationSet.from_dict(
            data=data,
            supercells=supercells,
        )

        # Add configurations to the ConfigurationSet
        scel = None
        for record in dataset:
            if scel != record.configuration.supercell:
                scel = record.configuration.supercell
                scel_record = casmconfig.SupercellRecord(scel)
                scel_name = scel_record.supercell_name
                print(f"Adding {len(dataset)} configurations from {scel_name}...")
                sys.stdout.flush()
            configurations.add_record(record)

    # Wait for workers to finish
    for i_process, process in enumerate(processes):
        process.join()

    # Return the generated supercells and configurations
    return (supercells, configurations)


def test_parallel_enum():
    """Example parallelizing enumeration by supercell."""
    xtal_prim = xtal_prims.FCC(
        r=0.5,
        occ_dof=["A", "B", "C"],
    )
    prim = casmconfig.Prim(xtal_prim)
    supercells, configurations = parallel_enum(
        prim=prim,
        min_vol=1,
        max_vol=10,
        n_workers=4,
    )
    assert len(supercells) == 87
    assert len(configurations) == 79260

import numpy as np

import libcasm.configuration as casmconfig
import libcasm.mapping.info as mapinfo
import libcasm.mapping.methods as mapmethods
import libcasm.xtal as xtal
import libcasm.xtal.prims as xtal_prims
import libcasm.xtal.structures as xtal_structures


def make_configuration(
    supercell: casmconfig.Supercell,
    occupation: list[int],
):
    config = casmconfig.Configuration(supercell)
    config.set_occupation(occupation)
    return config


def make_unique_mapped_structures(
    unmapped_structure: xtal.Structure,
    structure_mappings: mapinfo.StructureMappingResults,
    prim_factor_group: list[xtal.SymOp],
) -> list[xtal.Structure]:
    ### Find unique mapped structures
    unique_mapped_structures = []
    for i, smap in enumerate(structure_mappings):
        # print(f"~~~ {i} ~~~")
        # print(xtal.pretty_json(smap.to_dict()))
        # print()

        # Construct mapped_strucutre from structure mapping
        mapped_structure = mapmethods.make_mapped_structure(
            structure_mapping=smap,
            unmapped_structure=unmapped_structure,
        )

        # Append symmetrically unique mapped structures to list
        if len(unique_mapped_structures) == 0:
            unique_mapped_structures.append(mapped_structure)
            continue
        found = False
        for S in prim_factor_group:
            transformed_mapped_structure = S * mapped_structure
            if mapped_structure.is_equivalent_to(transformed_mapped_structure):
                found = True
                break
        if not found:
            unique_mapped_structures.append(mapped_structure)
    return unique_mapped_structures


def test_FCC_binary_prim_1(FCC_binary_prim):
    prim = casmconfig.Prim(FCC_binary_prim)
    T_prim = np.eye(3, dtype=int)
    T_conventional = np.array(
        [  # conventional FCC cubic cell
            [-1, 1, 1],
            [1, -1, 1],
            [1, 1, -1],
        ],
        dtype=int,
    )

    ### FCC primitive cell ###

    supercell = casmconfig.make_canonical_supercell(casmconfig.Supercell(prim, T_prim))
    default_configuration = casmconfig.Configuration(supercell)

    # convert Configuration to Structure
    structure = default_configuration.to_structure()
    assert structure.lattice() == supercell.superlattice
    assert np.allclose(
        structure.atom_coordinate_cart(),
        supercell.coordinate_cart(),
    )

    # read structure back in as a Configuration
    in_configuration = casmconfig.Configuration.from_structure(
        prim=prim,
        structure=structure,
    )
    assert default_configuration.supercell == in_configuration.supercell
    assert default_configuration == in_configuration

    ### A ###
    # set occupation and convert to Structure
    configuration_A1 = make_configuration(supercell, [0])
    structure = configuration_A1.to_structure()
    assert structure.atom_type() == ["A"]

    # read structure back in as a Configuration
    in_configuration = casmconfig.Configuration.from_structure(
        prim=prim,
        structure=structure,
    )
    assert in_configuration.supercell == supercell
    assert in_configuration.occupation.tolist() == [0]

    ### B ###
    # set occupation and convert to Structure
    configuration_B1 = make_configuration(supercell, [1])
    structure = configuration_B1.to_structure()
    assert structure.atom_type() == ["B"]

    # read structure back in as a Configuration
    in_configuration = casmconfig.Configuration.from_structure(
        prim=prim,
        structure=structure,
    )
    assert in_configuration.supercell == supercell
    assert in_configuration.occupation.tolist() == [1]

    ### Conventional FCC unit cell ###
    supercell = casmconfig.make_canonical_supercell(
        casmconfig.Supercell(prim, T_conventional)
    )
    casmconfig.SupercellRecord(supercell)

    # convert Configuration to Structure
    default_configuration = casmconfig.Configuration(supercell)

    # convert Configuration to Structure
    structure = default_configuration.to_structure()
    assert structure.lattice() == supercell.superlattice
    assert np.allclose(
        structure.atom_coordinate_cart(),
        supercell.coordinate_cart(),
    )

    # read structure back in as a Configuration
    in_configuration = casmconfig.Configuration.from_structure(
        prim=prim,
        structure=structure,
    )
    assert default_configuration.supercell == in_configuration.supercell
    assert default_configuration == in_configuration

    ### A3B1 ###
    # set occupation and convert to Structure
    configuration_A3B1 = make_configuration(supercell, [1, 0, 0, 0])
    structure = configuration_A3B1.to_structure()
    assert structure.atom_type() == ["B", "A", "A", "A"]

    # read structure back in as a Configuration
    in_configuration = casmconfig.Configuration.from_structure(
        prim=prim,
        structure=structure,
    )
    assert in_configuration.supercell == supercell
    assert in_configuration.occupation.tolist() == [1, 0, 0, 0]

    ### A1B3 ###
    # set occupation and convert to Structure
    configuration_A1B3 = make_configuration(supercell, [1, 1, 1, 0])
    structure = configuration_A1B3.to_structure()
    assert structure.atom_type() == ["B", "B", "B", "A"]

    # read structure back in as a Configuration
    in_configuration = casmconfig.Configuration.from_structure(
        prim=prim,
        structure=structure,
    )
    assert in_configuration.supercell == supercell
    assert in_configuration.occupation.tolist() == [1, 1, 1, 0]


def test_excluded_species_1():
    xtal_prim = xtal_prims.BCC(r=1.0, occ_dof=["A", "B", "Va"])
    prim = casmconfig.Prim(xtal_prim)

    T_conventional = np.array(
        [  # conventional BCC cubic cell
            [0, 1, 1],
            [1, 0, 1],
            [1, 1, 0],
        ],
        dtype=int,
    )
    T_supercell = T_conventional * 2

    ### FCC primitive cell ###

    supercell = casmconfig.make_canonical_supercell(
        casmconfig.Supercell(prim, T_supercell)
    )
    default_configuration = casmconfig.Configuration(supercell)
    default_configuration.set_occ(1, 2)

    # convert Configuration to Structure
    structure = default_configuration.to_structure()
    assert len(structure.atom_type()) == 15

    # convert Configuration to Structure
    structure = default_configuration.to_structure(excluded_species=[])
    assert len(structure.atom_type()) == 16


def test_bcc_hcp_mapping_conversions_1():
    prim = xtal_prims.BCC(r=1.0, occ_dof=["A", "B", "Va"])
    prim_factor_group = xtal.make_factor_group(prim)

    hcp_structure = xtal_structures.HCP(r=1.0, atom_type="A")

    ### Find structure mappings
    # - this generates multiple equivalent mappings
    structure_mappings = mapmethods.map_structures(
        prim,
        hcp_structure,
        prim_factor_group=prim_factor_group,
        max_vol=4,
        max_cost=1e20,
        min_cost=0.0,
        k_best=1,
    )

    ### Find unique mapped structures
    unique_mapped_structures = make_unique_mapped_structures(
        unmapped_structure=hcp_structure,
        structure_mappings=structure_mappings,
        prim_factor_group=prim_factor_group,
    )
    assert len(unique_mapped_structures) == 1

    ### Construct equivalent mapped structures
    prototype_structure = unique_mapped_structures[0]
    equivalent_mapped_structures = []
    for S in prim_factor_group:
        test = S * prototype_structure
        found = False
        for existing in equivalent_mapped_structures:
            if test.is_equivalent_to(existing):
                found = True
                break
        if not found:
            equivalent_mapped_structures.append(test)

    # print("### Equivalent mapped structures ###")
    # for i, equiv in enumerate(equivalent_mapped_structures):
    #     print(f"~~~ {i} ~~~")
    #     print(xtal.pretty_json(equiv.to_dict()))
    #     print()
    assert len(equivalent_mapped_structures) == 12

    ### Construct equivalent configuration in a Prim with strain and displacement DoF
    BCC_strain_disp_prim = casmconfig.Prim(
        xtal_prims.BCC(
            r=1.0,
            occ_dof=["A", "B", "Va"],
            local_dof=[xtal.DoFSetBasis("disp")],
            global_dof=[xtal.DoFSetBasis("Hstrain")],
        )
    )
    mapped_configuration = casmconfig.make_canonical_configuration(
        configuration=casmconfig.Configuration.from_structure(
            prim=BCC_strain_disp_prim,
            structure=prototype_structure,
        ),
        in_canonical_supercell=True,
    )

    ### Make the fully commensurate superdupercell
    superduperlattice = xtal.make_superduperlattice(
        lattices=[mapped_configuration.supercell.superlattice],
        mode="fully_commensurate",
        point_group=BCC_strain_disp_prim.crystal_point_group.elements,
    )
    T = xtal.make_transformation_matrix_to_super(
        superlattice=superduperlattice,
        unit_lattice=BCC_strain_disp_prim.xtal_prim.lattice(),
    )
    superdupercell = casmconfig.make_canonical_supercell(
        casmconfig.Supercell(
            prim=BCC_strain_disp_prim,
            transformation_matrix_to_super=T,
        ),
    )

    ### Put the mapped configuration into the fully commensurate superdupercell
    prototype_configuration = casmconfig.copy_configuration(
        motif=mapped_configuration,
        supercell=superdupercell,
    )

    ### Generate equivalent configurations
    equivalent_configurations = casmconfig.make_equivalent_configurations(
        configuration=prototype_configuration,
    )
    # print("### Equivalent configurations ###")
    # for i, equiv in enumerate(equivalent_configurations):
    #     print(f"~~~ {i} ~~~")
    #     print(xtal.pretty_json(equiv.to_dict()))
    #     structure = equiv.to_structure()
    #     print(xtal.pretty_json(structure.to_dict()))
    #     print()
    assert len(equivalent_configurations) == 12


def test_bcc_hcp_mapping_conversions_2():
    prim = xtal_prims.BCC(r=1.0, occ_dof=["A", "B", "Va"])
    prim_factor_group = xtal.make_factor_group(prim)

    _hcp_structure = xtal_structures.HCP(
        r=1.0,
        atom_type="A",
        global_properties={"energy": np.array([[0.1]]).transpose()},
    )
    hcp_structure = xtal.Structure(
        lattice=_hcp_structure.lattice(),
        atom_coordinate_frac=_hcp_structure.atom_coordinate_frac(),
        atom_type=["A", "B"],
        global_properties=_hcp_structure.global_properties(),
    )

    ### Find structure mappings
    # - this generates multiple equivalent mappings
    structure_mappings = mapmethods.map_structures(
        prim,
        hcp_structure,
        prim_factor_group=prim_factor_group,
        max_vol=4,
        max_cost=1e20,
        min_cost=0.0,
        k_best=1,
    )

    ### Find unique mapped structures
    unique_mapped_structures = make_unique_mapped_structures(
        unmapped_structure=hcp_structure,
        structure_mappings=structure_mappings,
        prim_factor_group=prim_factor_group,
    )
    assert len(unique_mapped_structures) == 1

    ### Construct equivalent mapped structures
    prototype_structure = unique_mapped_structures[0]
    equivalent_mapped_structures = []
    for S in prim_factor_group:
        test = S * prototype_structure
        found = False
        for existing in equivalent_mapped_structures:
            if test.is_equivalent_to(existing):
                found = True
                break
        if not found:
            equivalent_mapped_structures.append(test)

    # print("### Equivalent mapped structures ###")
    # for i, equiv in enumerate(equivalent_mapped_structures):
    #     print(f"~~~ {i} ~~~")
    #     print(xtal.pretty_json(equiv.to_dict()))
    #     print()
    assert len(equivalent_mapped_structures) == 12

    ### Construct configuration with properties,
    # including strain, displacement, and energy
    BCC_prim = casmconfig.Prim(prim)

    config_w_props = casmconfig.ConfigurationWithProperties.from_structure(
        prim=BCC_prim,
        structure=prototype_structure,
    )
    assert len(config_w_props.local_properties) == 1
    assert "disp" in config_w_props.local_properties
    assert len(config_w_props.global_properties) == 2
    assert "energy" in config_w_props.global_properties
    assert "Ustrain" in config_w_props.global_properties

    # print("### Mapped configuration w/ properties ###")
    # print(xtal.pretty_json(config_w_props.to_dict()))

    ### Make canonical mapped configuration w/ properties
    S = casmconfig.to_canonical_configuration(
        configuration=config_w_props.configuration
    )
    S * config_w_props

    # print("### Mapped canonical configuration w/ properties ###")
    # print(xtal.pretty_json(canonical_config_w_props.to_dict()))

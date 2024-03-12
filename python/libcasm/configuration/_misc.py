import copy
import functools
import math
from typing import Optional

import numpy as np

import libcasm.casmglobal


@functools.total_ordering
class ApproximateFloatArray:
    def __init__(
        self,
        arr: np.ndarray,
        abs_tol: float = libcasm.casmglobal.TOL,
    ):
        """Store an array that will be compared lexicographically up to a given
        absolute tolerance using math.isclose

        Parameters
        ----------
        arr: np.ndarray
            The array to be compared

        abs_tol: float = :data:`~libcasm.casmglobal.TOL`
            The absolute tolerance
        """
        if not isinstance(arr, np.ndarray):
            raise TypeError(
                "Error in ApproximateFloatArray: arr must be a numpy.ndarray"
            )
        self.arr = arr
        self.abs_tol = abs_tol

    def __eq__(self, other):
        if len(self.arr) != len(other.arr):
            return False
        for i in range(len(self.arr)):
            if not math.isclose(self.arr[i], other.arr[i], abs_tol=self.abs_tol):
                return False
        return True

    def __lt__(self, other):
        if len(self.arr) != len(other.arr):
            return len(self.arr) < len(other.arr)
        for i in range(len(self.arr)):
            if not math.isclose(self.arr[i], other.arr[i], abs_tol=self.abs_tol):
                return self.arr[i] < other.arr[i]
        return False


def is_canonical_order_parameters(
    eta: np.ndarray,
    dof_space_rep: list[np.ndarray[np.float64]],
    abs_tol: float = libcasm.casmglobal.TOL,
) -> bool:
    _eta = ApproximateFloatArray(eta, abs_tol=abs_tol)
    for M in dof_space_rep:
        _test = ApproximateFloatArray(M * eta, abs_tol=abs_tol)
        if _test > _eta:
            return False
    return True


def make_canonical_order_parameters(
    eta: np.ndarray,
    dof_space_rep: list[np.ndarray[np.float64]],
    abs_tol: float = libcasm.casmglobal.TOL,
) -> np.ndarray:
    _most_canonical = ApproximateFloatArray(eta, abs_tol=abs_tol)
    for M in dof_space_rep:
        _test = ApproximateFloatArray(M @ eta, abs_tol=abs_tol)
        if _test > _most_canonical:
            _most_canonical = copy.deepcopy(_test)
    return _most_canonical.arr


def equivalent_order_parameters_index(
    canonical_eta_list: list[np.ndarray],
    eta: np.ndarray,
    dof_space_rep: list[np.ndarray[np.float64]],
    abs_tol: float = libcasm.casmglobal.TOL,
) -> Optional[int]:
    """Find index of first equivalent order parameters in a list, if any

    Parameters
    ----------
    canonical_eta_list: list[np.ndarray]
        A list of points, represented in the DoFSpace basis, in canonical form
    eta: np.ndarray
        A point, represented in the DoFSpace basis.
    dof_space_rep: list[numpy.ndarray[numpy.float64[subspace_dim, subspace_dim]]]
        Elements, `M`, of `dof_space_rep` transform subspace vectors according to
        ``x_subspace_after = M @ x_subspace_before``.
    abs_tol: float = :data:`~libcasm.casmglobal.TOL`
        The absolute tolerance

    Returns
    -------
    index : Optional[int]
        Returns the index of the first point in `eta_list` that is symmetrically
        equivalent to `eta`. Returns None if none are equivalent.
    """
    _test = ApproximateFloatArray(
        make_canonical_order_parameters(
            eta=eta,
            dof_space_rep=dof_space_rep,
            abs_tol=abs_tol,
        ),
        abs_tol=abs_tol,
    )
    for i, x in enumerate(canonical_eta_list):
        _x = ApproximateFloatArray(x, abs_tol=abs_tol)
        if _x == _test:
            return i
    return None

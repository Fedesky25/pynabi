from typing import Literal as _L, Union as _U

from pynabi._common import StampCollection
from ._common import Stampable as _S, _pos_int, OneLineStamp as _OLS, IndexedWithDefault as _IWD
from .units import Energy as _En
from enum import Enum as _E


__all__ = ["SCFDirectMinimization", "SCFMixing", "NonSelfConsistentCalc", "ToleranceOn", "EnergyCutoff", "MaxSteps"]


class SCFDirectMinimization(_S):
    def __init__(self) -> None:
        super().__init__()
    
    def stamp(self, index: int):
        return f"iscf{index or ''} 0"
    
    def compatible(self, coll: StampCollection):
        tol = coll.get(Tolerance)
        assert tol is not None, "SCF direct minimization requires one tolerance to be spcified"


class SCFMixing(_IWD, default="Pulay", prop="iscf"):
    """Usual ground state (GS) calculations or for structural relaxations, where the potential has to be determined self-consistently"""

    def __init__(self, density: bool = False) -> None:
        super().__init__()
        self._den = density
    
    def Simple(self):
        self._index = 2 + 10*int(self._den)
        return self

    def Anderson(self, onPrevious: bool = False):
        """Anderson mixing of the potential/density.\n
        If `onPrevious` is True, the mixing is based on the two previous iterations"""
        self._index = (4 if onPrevious else 3) + 10*int(self._den)
        return self

    def CGBased(self, alt: bool = False):
        """CG based on the minimum of the energy with respect to the potential/density.\n
        If `alt` is True, it uses an alternative implementation (still in development)"""
        self._index = (6 if alt else 5) + 10*int(self._den)
        return self

    def Pulay(self, iterations: int = 7):
        """Pulay mixing of the potential/density based on the `iterations` previous iterations"""
        assert _pos_int(iterations), "Number of iterations used in Pulay mixing must be a positive integer"
        self._index = 7 + 10*int(self._den)
        self._extra = { "npulayit": iterations }
        return self
    
    def compatible(self, coll: StampCollection):
        tol = coll.get(Tolerance)
        assert tol is not None, "SCF mixing requires one tolerance to be specified"


class NonSelfConsistentCalc(_S):
    def __init__(self, i: _L[-1,-2,-3] = -2) -> None:
        assert type(i) is int and -3 <= i <= -1, "Non self consistend calculation index must be -1, -2, or -3"
        self._i = i
    
    def stamp(self, index: int):
        return f"iscf{index or ''} {self._i}"
    

class ToleranceOn(_E):
    EnergyDifference = "dfe"
    ForceDifference = "dff"
    ForceRelativeDifference = "drff"
    PotentialResidual = "vrs"
    WavefunctionSquaredResidual = "wfr"

    def __call__(self, value: float):
        return Tolerance(value, self.value)


class Tolerance(_OLS):
    """Do not use this class directly: prefer ToleranceOn"""
    name = "tol"


class EnergyCutoff(_OLS):
    """Used to define the kinetic energy cutoff which controls the number of planewaves at given k point. The allowed plane waves are those with kinetic energy lower than ecut, which translates to the following constraint on the planewave vector G in reciprocal space"""
    name = "ecut"
    def __init__(self, value: _U[float,_En]) -> None:
        e = _En.sanitize(value);
        assert e._v > 0, "cutoff energy must be positive"
        super().__init__(str(e))


class MaxSteps(_OLS):
    """The maximum number of cycles (or “iterations”) in a SCF or non-SCF run. 
    
    Full convergence from random numbers is usually achieved in 12-20 SCF iterations. Each can take from minutes to hours. In certain difficult cases, usually related to a small or zero band gap or magnetism, convergence performance may be much worse
    
    For non-self-consistent runs, it governs the number of cycles of convergence for the wavefunctions for a fixed density and Hamiltonian.

    NOTE that a choice of nstep = 0 is permitted; this will either read wavefunctions from disk and compute the density, the total energy and stop, or else will initialize randomly the wavefunctions and compute the resulting density and total energy. This is provided for testing purposes.
    """

    name = "nstep"
    def __init__(self, value: int) -> None:
        assert type(value) is int and value >= 0, "Number of steps must be an integer greater than or equal to 0"
        super().__init__(value)


_exclusives = set((SCFMixing, SCFDirectMinimization, NonSelfConsistentCalc))
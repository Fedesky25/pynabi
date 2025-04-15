"""
WARNING: do not import this file directly!
"""

from pynabi._common import Stampable, StampCollection, _pos_int, OneLineStamp, IndexedWithDefault 
from pynabi.units.internal import Energy
from typing import Literal, Union as Union
from enum import Enum


class SCFDirectMinimization(Stampable):
    def __init__(self) -> None:
        super().__init__()
    
    def stamp(self, index: int):
        return f"iscf{index or ''} 0"
    
    def compatible(self, coll: StampCollection):
        tol = coll.get(Tolerance)
        assert tol is not None, "SCF direct minimization requires one tolerance to be spcified"


class SCFMixing(IndexedWithDefault, default="Pulay", prop="iscf"):
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


class NonSelfConsistentCalc(Stampable):
    def __init__(self, i: Literal[-1,-2,-3] = -2) -> None:
        assert type(i) is int and -3 <= i <= -1, "Non self consistend calculation index must be -1, -2, or -3"
        self._i = i
    
    def stamp(self, index: int):
        return f"iscf{index or ''} {self._i}"
    

class ToleranceOn(Enum):
    EnergyDifference = "dfe"
    ForceDifference = "dff"
    ForceRelativeDifference = "drff"
    PotentialResidual = "vrs"
    WavefunctionSquaredResidual = "wfr"

    def __call__(self, value: float):
        return Tolerance(value, self.value)


class Tolerance(OneLineStamp):
    """Do not use this class directly: prefer ToleranceOn"""
    name = "tol"


class EnergyCutoff(OneLineStamp):
    """Used to define the kinetic energy cutoff which controls the number of planewaves at given k point. The allowed plane waves are those with kinetic energy lower than ecut, which translates to the following constraint on the planewave vector G in reciprocal space"""
    name = "ecut"
    def __init__(self, value: Union[float,Energy]) -> None:
        e = Energy.sanitize(value);
        assert e._v > 0, "cutoff energy must be positive"
        super().__init__(str(e))


class MaxSteps(OneLineStamp):
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
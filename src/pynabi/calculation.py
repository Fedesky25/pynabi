from typing import Optional as _O, Literal as _L, Union as _U
from ._common import Stampable as _S, _pos_int, OneLineStamp as _OLS, IndexedWithDefault as _IWD
from .units import Energy as _En
from enum import Enum as _E


__all__ = ["SCFDirectMinimization", "SCFMixing", "NonSelfConsistentCalc", "ToleranceOn", "EnergyCutoff", "MaxSteps"]


class SCFDirectMinimization(_S):
    def __init__(self) -> None:
        super().__init__()
    
    def stamp(self, index: int):
        return f"iscf{index or ''} 0"


class SCFMixing(_IWD, default="Pulay", prop="iscf"):
    """Usual ground state (GS) calculations or for structural relaxations, where the potential has to be determined self-consistently"""

    def __init__(self, density: bool = False) -> None:
        super().__init__()
        self._den = density
    
    def Simple(self):
        self._index = 2 + 10*int(self._den)
        return self

    def Anderson(self, onPrevious: bool = False):
        self._index = (4 if onPrevious else 3) + 10*int(self._den)
        return self

    def CGBased(self, alt: bool = False):
        self._index = (6 if alt else 5) + 10*int(self._den)
        return self

    def Pulay(self, iterations: int = 7):
        assert _pos_int(iterations), "Number of iterations used in Pulay mixing must be a positive integer"
        self._index = 7 + 10*int(self._den)
        self._extra = { "npulayit": iterations }
        return self


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
    name = "ecut"
    def __init__(self, value: _U[float,_En]) -> None:
        e = _En.sanitize(value);
        assert e._v > 0, "cutoff energy must be positive"
        super().__init__(str(e))


class MaxSteps(_OLS):
    name = "nstep"
    def __init__(self, value: int) -> None:
        assert type(value) is int and value >= 0, "Number of steps must be an integer greater than or equal to 0"
        super().__init__(value)


_exclusives = set((SCFMixing, SCFDirectMinimization, NonSelfConsistentCalc))
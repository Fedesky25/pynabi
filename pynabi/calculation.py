from typing import Optional as _O, Literal as _L
from ._common import Stampable as _S, _pos_int, OneLineStamp as _OLS
from .units import EUnit as _U
from enum import Enum as _E


__all__ = ["SCFDirectMinimization", "SCFMixing", "NonSelfConsistentCalc", "ToleranceOn", "EnergyCutoff", "MaxSteps"]


class SCFDirectMinimization(_S):
    def __init__(self, steps: _O[int] = None) -> None:
        self._s = steps
    
    def stamp(self, index: int):
        s = index or ''
        r = f"iscf{s} 0"
        if self._s is not None:
            r = r + f"\nnsteps{s} {self._s}"
        return r


class SCFMixing(_S):
    def __init__(self, density: bool = False):
        self._i: _O[int] = None
        self._d = density
        self._p: _O[int] = None

    def Simple(self):
        self._i = 2
        return self

    def Anderson(self, onPrevious: bool = False):
        self._i = 4 if onPrevious else 3
        return self

    def CGBased(self, alt: bool = False):
        self._i = 6 if alt else 5
        return self

    def Pulay(self, iterations: int = 7):
        assert _pos_int(iterations), "Number of iterations used in Pulay mixing must be a positive integer"
        self._i = 7
        self._p = iterations
        return self
    
    def stamp(self, index: int):
        s = index or ''
        res: list[str] = []
        if self._i is not None:
            res.append(f"iscf{s} {self._i+10 if self._d else self._i}")
        if self._p is not None:
            res.append(f"npulayit{s} {self._p}")
        if len(res) == 0:
            print("Warning: empty SCFMixing")
        return '\n'.join(res)


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
    def __init__(self, value: float, unit: _U = _U.getDefault()) -> None:
        assert value > 0, "cutoff energy must be positive"
        super().__init__(f"{value} {unit.name}")


class MaxSteps(_OLS):
    name = "nstep"
    def __init__(self, value: int) -> None:
        assert type(value) is int and value >= 0, "Number of steps must be an integer greater than or equal to 0"
        super().__init__(value)
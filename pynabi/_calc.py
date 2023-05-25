from enum import Enum
from ._units import Unit
from ._common import OneLineStamp

__all__ = ["ToleranceOn", "UseWaveletBasis", "SCFProcedure", "StepNumber", "EnergyCutoff"]

class ToleranceOn(Enum):
    EnergyDifference = "dfe"
    ForceDifference = "dff"
    ForceRelativeDifference = "drff"
    PotentialResidual = "vrs"
    WavefunctionSquaredResidual = "wfr"

    def __call__(self, value: float):
        return Tolerance(value, self.value)


class Tolerance(OneLineStamp):
    name = "tol"
    longname = "tolerance"
    

class UseWaveletBasis(OneLineStamp):
    name = "usewvl"
    longname = "use wavelet basis"
    def __init__(self, use: bool) -> None:
        super().__init__(int(use))


class SCFProcedure(OneLineStamp):
    name = "iscf"
    longname = "SCF Procedure"
    def __init__(self, algorithm: int) -> None:
        super().__init__(algorithm)
        assert -3 <= algorithm <= 17, "SCF algorithm must go from -3 to 17" 


class StepNumber(OneLineStamp):
    name = "nstep"
    longname = "Number of steps"
    def __init__(self, value: int) -> None:
        super().__init__(value)
        assert type(value) is int and value >= 0, "Number of steps must be an integer greater than or equal to 0"


class EnergyCutoff(OneLineStamp):
    name = "ecut"
    longname = "Energy cutoff"
    def __init__(self, value: float, unit: Unit = Unit.Hartree) -> None:
        super().__init__(f"{value} {unit.name}")
        assert unit.value == 1, "unit of cutoff energy must be energy"
        assert value > 0, "cutoff energy must be positive"


# class CalculationCycle(Stampable):
#     def __init__(self, algorith: int, steps: int) -> None:
#         super().__init__()
#         self.algorithm = algorith
#         self.steps = steps

#     def stamp(self, index: int):
#         s = index or ''
#         return f"{sectionTitle(index, 'SCF procedure')}\niscf{s} {self.algorithm}\nnstep{s} {self.steps}"
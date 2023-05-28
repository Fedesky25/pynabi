from enum import Enum as _E


__all__ = ["LUnit", "Pos3D", "Length", "EUnit", "EUnit"]

class LUnit(_E):
    Bohr = 0
    nm = 1
    Angstrom = 2

    def setAsDefault(self):
        global _def_LU
        _def_LU = self
    
    @staticmethod
    def getDefault():
        return _def_LU

_def_LU = LUnit.Bohr


class Pos3D:
    def __init__(self, x: float, y: float, z: float, unit: LUnit = _def_LU) -> None:
        assert type(unit) is LUnit, f"{unit} is not valid length unit"
        self.x = x
        self.y = y
        self.z = z
        self.u = unit
    
    def __str__(self) -> str:
        return f"{self.x} {self.y} {self.z} {self.u.name}"
    
    @staticmethod
    def uniform(l: 'Length'):
        return Pos3D(l._v,l._v,l._v,l._u) # type: ignore


class Length:
    def __init__(self, value: float, unit: LUnit = _def_LU) -> None:
        assert type(unit) is LUnit, f"{unit} is not  avlid length unit"
        self._v = value
        self._u = unit

    def __str__(self) -> str:
        return f"{self._v} {self._u.name}"


class EUnit(_E):
    Ha = 0
    eV = 1
    Ry = 3
    K = 4

    def setAsDefault(self):
        global _def_EU
        _def_EU = self
    
    @staticmethod
    def getDefault():
        return _def_EU

_def_EU = EUnit.Ha


class Energy:
    def __init__(self, value: float, unit: EUnit = _def_EU) -> None:
        assert type(unit) is EUnit, f"{unit} is not valid energy unit"
        self._v = value
        self._u = unit

    def __str__(self) -> str:
        return f"{self._v} {self._u.name}"
from enum import Enum as Enum


__all__ = ["Unit", "Pos3D", "Length"]


class Unit(Enum):
    Hartree = 1
    eV = 1
    meV = 1
    Ry = 1
    K = 1
    Bohr = 2
    nm = 2
    angstrom = 2


class Pos3D:
    def __init__(self, x: float, y: float, z: float, unit: Unit = Unit.Bohr) -> None:
        self.x = x
        self.y = y
        self.z = z
        self.u = unit
        assert unit.value == 2, "Position in 3D must have length dimension"
    
    def __str__(self) -> str:
        return f"{self.x} {self.y} {self.z} {self.u.name}"
    
    @staticmethod
    def uniform(v: float, unit: Unit = Unit.Bohr):
        return Pos3D(v,v,v,unit)
    

class Length:
    def __init__(self, value: float, unit: Unit) -> None:
        self.v = value
        self.u = unit
        assert unit.value == 2, "Length must have length unit"

    def __str__(self) -> str:
        return f"{self.v} {self.u.name}"

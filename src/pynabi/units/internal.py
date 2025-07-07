"""
WARNING: do not import this file directly!
"""

from pynabi._common import Vec3D
from typing import Self, Union, Tuple, Optional, Any


class AbMeasure:
    _U: Tuple[Tuple[float, str],...] = ()
    _R = (1.0, 0)

    def __init__(self, value: float, unit: int) -> None:
        """Do not use this constructor"""
        self._v = value
        self._u = unit

    def __mul__(self, other: Union[int,float]) -> Self:
        if type(other) is int or type(other) is float:
            return type(self)(self._v*other, self._u) # type: ignore
        else:
            raise NotImplementedError(f"Quantity of type {type(self).__name__} can only be multiplied by a scalar number")

    def __rmul__(self, other: Union[int,float]) -> Self:
        return type(self).__mul__(self, other)

    def __add__(self, other: Self):
        if type(other) is not type(self):
            raise NotImplementedError(f"Cannot add quantities of different type ({type(self).__name__} and {type(other).__name__})")
        U = type(self)._U
        delta = other._v * U[other._u][0] / U[self._u][0]
        return type(self)(self._v + delta, self._u);

    def __sub__(self, other: Self) -> Self:
        if type(other) is not type(self):
            raise NotImplementedError(f"Cannot subtract quantities of different type ({type(self).__name__} and {type(other).__name__})")
        U = type(self)._U
        delta = other._v * U[other._u][0] / U[self._u][0]
        return type(self)(self._v - delta, self._u);

    def __truediv__(self, other: Union[float,int]):
        if type(other) is float or type(other) is int:
            return type(self)(self._v/other, self._u)
        else:
            n = type(self).__name__
            raise NotImplementedError(f"Quantity of type {n} can only be divided by a scalar number or a {n}")
        
    def __rtruediv__(self, other: Self):
        if type(other) is type(self):
            U = type(self)._U
            return (other._v / self._v) * (U[other._u][0]/U[self._u][0])
        else:
            n = type(self).__name__
            raise NotImplementedError(f"Quantity of type {n} can only be divided by a scalar number or a {n}")

    def __neg__(self):
        return type(self)(-self._v, self._u)

    def __str__(self):
        return f"{self._v} {type(self)._U[self._u][1]}"
    
    def setAsReference(self):
        type(self)._R = (self._v, self._u)
    
    @classmethod
    def fromReference(cls, value: float, ref: Optional[Self]) -> Self:
        assert type(value) is float or type(value) is int, f"Value of {cls.__name__} must be a number"
        if ref is None:
            return cls(value*cls._R[0], cls._R[1])
        else:
            assert type(ref) is cls, f"Reference for a {cls.__name__} mus be a {cls.__name__} itself"
            return value*ref

    @classmethod
    def sanitize(cls, value: Any) -> Self:
        t = type(value)
        if t is float or t is int:
            return cls(value*cls._R[0], cls._R[1])
        elif t is cls:
            return value
        else:
            raise TypeError(f"Cannot construct {cls.__name__} from {value} (of type {t})")
        
    @classmethod
    def getDefaultReference(cls):
        return cls(cls._R[0], cls._R[1])


class Length(AbMeasure):
    _U = (
        (0.529177249, "Bohr"), 
        (1.0, "Angstrom"), 
        (10.0, "nm")
    )


Bohr = Length(1.0, 0)
Ang = Length(1.0, 1)
nm = Length(1.0, 2)


class Energy(AbMeasure):
    _U = (
        (27.2114, "Ha"),
        (1.0, "eV"),
        (13.6057039763, "Ry"),
        (8.617328149741e-5, "K")
    )


Ha = Energy(1.0, 0)
eV = Energy(1.0, 1)
Ry = Energy(1.0, 2)
Kelvin = Energy(1.0, 3)

def _ratio(x: Union[Length,float,int], unit: int):
    if type(x) is Length:
        return x._v*Length._U[unit][0]/Length._U[x._u][0]
    elif type(x) is float or type(x) is int:
        return x
    else:
        raise TypeError("Position component must be a scalar or length")

class Pos3D:
    def __init__(self, x: Union[Length,float], y: Union[Length,float], z: Union[Length,float], unit: Optional[Length] = None) -> None:
        m = 1.0
        if unit is None:
            m = Length._R[0]
            self.u = Length._R[1]
        elif type(unit) is Length:
            m = unit._v
            self.u = unit._u
        else:
            raise TypeError("3D Position reference must be a Length itself")

        self.x = _ratio(x, self.u)*m
        self.y = _ratio(y, self.u)*m
        self.z = _ratio(z, self.u)*m
    
    def __str__(self) -> str:
        return f"{self.x} {self.y} {self.z} {Length._U[self.u][1]}"
    
    @staticmethod
    def uniform(l: Union[float,'Length']):
        return Pos3D(1,1,1,Length.sanitize(l))
    
    @staticmethod
    def sanitize(v: Any) -> 'Pos3D':
        t = type(v)
        if t is Pos3D:
            return v
        elif t is Vec3D:
            return Pos3D(v.x,v.y,v.z)
        else:
            raise TypeError(f"Cannot construct Pos3D from {v} (of type {t})")
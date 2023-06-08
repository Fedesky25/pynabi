from typing import TypeVar, Type, Tuple, Callable


__all__ = ["Vec3D", "SKO"]

class Vec3D:
    def __init__(self, x: float, y: float, z: float) -> None:
        self.x = x
        self.y = y
        self.z = z
    
    def __str__(self) -> str:
        return "{} {} {}".format(self.x, self.y, self.z)
    
    @staticmethod
    def uniform(v: float = 1.0):
        return Vec3D(v,v,v)
    
    @staticmethod
    def zero():
        return Vec3D(0,0,0)

 
S = TypeVar("S", bound="Stampable")

class StampCollection:
    def __init__(self, siblings: 'dict[type,Stampable]', base: 'dict[type,Stampable]') -> None:
        self._s = siblings
        self._b = base
    
    def get(self, t: Type[S]) -> S | None:
        s = self._s.get(t);
        if s is not None:
            s = self._b.get(t)
        return s # type: ignore
    
    def nextto(self, t: Type['Stampable']):
        return t in self._s


class Stampable:
    def stamp(self, index: int):
        raise NotImplementedError(f"{type(self).__name__} has not implemented stamp")
    
    def compatible(self, coll: 'StampCollection'):
        pass
    
    def __str__(self):
        return self.stamp(0)


class OneLineStamp(Stampable):
    name = ''
    
    def __init__(self, value, suffix: str = '') -> None:
        super().__init__()
        self.value = value
        self.suffix = suffix

    def stamp(self, index: int):
        t = type(self)
        return f"{t.name}{self.suffix}{index or ''} {self.value}"


class Later(object):
    _instance = None
    def __new__(cls):
        if cls._instance is None:
            cls._instance = object.__new__(cls)
        return cls._instance


class CanDelay(Stampable): 
    _delayables: Tuple[Tuple[str, str, Callable],...] = ()
    def __init__(self, *values):
        self._dv = values;
        for i,v in enumerate(values):
            if v is not Later._instance:
                self._delayables[i][2](v)

    def _doesDelay(self, i: int):
        return self._dv[i] is Later._instance
    
    def stamp(self, index: int):
        res: list[str] = []
        s = index or ''
        for i,v in self._dv:
            if v is not Later._instance:
                res.append(f"{self._delayables[i][0]}{s} {v}")
        return '\n'.join(res)


class Delayed(Stampable):
    def __init__(self, cls: Type[CanDelay], index: int, value) -> None:
        super().__init__() 
        cls._delayables[index][2](value)
        self.c = cls
        self.i = index
        self.v = value
    
    def compatible(self, coll: StampCollection):
        v = coll.get(self.c)
        if v is None:
            raise TypeError(f"{self.c._delayables[self.i][1]} definition requires {self.c.__name__} definition before")
        if not v._doesDelay(self.i):
            raise ValueError(f"{self.c.__name__} already defines {self.c._delayables[self.i][1]}")
    
    def stamp(self, index: int):
        return f"{self.c._delayables[self.i][0]}{index or ''} {self.v}"


class SKO:
    """Single K Occupation"""
    def __init__(self, k: Vec3D, *occupations: float) -> None:
        self.k = k
        self.occ = occupations


def sectionTitle(index: int, name: str):
    if index > 0:
        return f"\n# DS{index} - {name}"
    else:
        return f"\n# {name}"
    
    
def _pos_int(v):
    return type(v) is int and v > 0
from typing import TypeVar, Type, Tuple, Callable, Any, TypeGuard

from pynabi._common import StampCollection


__all__ = ["Vec3D"]


class Singleton(object):
    _instance = None
    def __new__(cls):
        if cls._instance is None:
            cls._instance = object.__new__(cls)
        return cls._instance
    
    @classmethod
    def eq(cls, value):
        return type(value) is cls or value is cls


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
        if s is None:
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


class Later(Singleton):
    pass


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
        for i,v in enumerate(self._dv):
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
            raise TypeError(f"{self.getName()} definition requires {self.c.__name__} definition before")
        if not v._doesDelay(self.i):
            raise ValueError(f"{self.c.__name__} already defines {self.getName()}")
    
    def stamp(self, index: int):
        return f"{self.getProp()}{index or ''} {self.v}"

    def getProp(self):
        return self.c._delayables[self.i][0]

    def getName(self):
        return self.c._delayables[self.i][1]


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
    

def decorate_IWD_method(fn):
    def method(self: 'IndexedWithDefault', *args):
        assert self._index is None, f"Multiple definitions for {type(self).__name__}"
        return fn(self, *args)
    return method


class IndexedWithDefault(Stampable):
    _default: str = ""
    _prop: str = ""
    
    def __init__(self) -> None:
        super().__init__()
        self._index: int|None = None
        self._extra: dict[str, Any] = {}
    
    def __init_subclass__(cls, default: str, prop: str = "") -> None:
        super().__init_subclass__()
        cls._prop = prop
        if not callable(getattr(cls, default, 0)): # 0 is not callable
            raise AttributeError(f"Property {default} of {cls.__name__} is not callable")
        for (name,method) in cls.__dict__.items():
            if callable(method) and not hasattr(IndexedWithDefault, name):
                setattr(cls, name, decorate_IWD_method(method))
    
    def stamp(self, index: int):
        s = index or ''
        if self._index is None:
            getattr(self,self._default)()
        res = f"{self._prop}{s} {self._index}"
        extra = '\n'.join(f"{k}{s} {v}" for (k,v) in self._extra.items())
        if len(extra) > 0:
            res += '\n' + extra
        return res


def _pos_int(v) -> TypeGuard[int]:
    return type(v) is int and v > 0


def _pos0_int(v) -> TypeGuard[int]:
    return type(v) is int and v >= 0


def _pos_num(v) -> TypeGuard[float|int]:
    return type(v) is float or type(v) is int and v > 0


def _pos0_num(v) -> TypeGuard[float|int]:
    return type(v) is float or type(v) is int and v >= 0
from typing import Union, Tuple, Literal, Callable, Type, Optional
from enum import Enum

from pynabi._common import StampCollection
from ._common import Stampable


__all__ = ["SpinType", "Metal", "Smearing", "Semiconductor", "TwoQuasiFermilevels", "OccupationPerBand"]


class SpinPolarization(Stampable):
    def __init__(self, pol: Literal[1,2], inor: Literal[1,2], den: Literal[1,2,4]) -> None:
        """DO NOT USE: use precomputed spin polazitation instead"""
        super().__init__()
        self.polarizationNumber = pol
        self.spinorNumber = inor
        self.density = den

    def stamp(self, index: int):
        s = index or ''
        return f"""nsppol{s} {self.polarizationNumber}
nspinor{s} {self.spinorNumber}
nspden{s} {self.density}"""
    

class SpinType:
    """Possible spin polarization options\n\nDO NOT MODIFY"""
    Unpolarized = SpinPolarization(1,1,1)
    Antiferromagnetic = SpinPolarization(1,1,2)
    SpinOrbitCoupling = SpinPolarization(1,2,1)
    NonCollinearMagnetism = SpinPolarization(1,2,4)
    Polarized = SpinPolarization(2,1,2)


class Smearing(Enum):
    FermiDirac = 3
    Marzari5634 = 4
    Marzari8165 = 5
    MethfesselPaxton = 6
    Gaussian = 7
    Uniform = 8


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


def _checkbands(v):
    assert type(v) is int and v > 0, "Number of bands must be positive integer"

def _checktsmear(v):
    assert type(v) is float and 0 <= v <= 1, "Smearing broadening must be a float between 0 and 1"


class Metal(CanDelay):
    _delayables = ( 
        ("tsmear", "Smearing broadening", _checktsmear),
        ("nbands", "Number of bands", _checkbands)
    )

    def __init__(self, smearing: Smearing, broadening: Union[float,Later] = Later(), bands: Union[int,Later] = Later()) -> None:
        super().__init__(broadening, bands)
        self._opt = smearing.value

    def stamp(self, index: int):
        return f"occopt{index or ''} {self._opt}\n{super().stamp(index)}"

    @classmethod
    def setBroadening(cls, value: float):
        return Delayed(cls, 0, value)

    @classmethod
    def setBands(cls, value: int):
        return Delayed(cls, 1, value)
    

class Semiconductor(CanDelay):
    _delayables = (("nbands", "Number of bands", _checkbands),)

    def __init__(self, spinMagnetizationTarget: Optional[float], bands: Union[int,Later] = Later()):
        super().__init__(bands)
        self._smt = spinMagnetizationTarget
    
    def compatible(self, coll: StampCollection):
        spin = coll.get(SpinPolarization)
        if spin is not None and spin.polarizationNumber == 2:
            assert self._smt is not None, "Semiconducotor with polarized spin must specify the spin magnetization target"

    def stamp(self, index: int):
        res = [f"occopt{index or ''} 1"]
        if self._smt is not None:
            res.append(f"spinmagntarget{index or ''} {self._smt}")
        s = super().stamp(index)
        if len(s) > 0:
            res.append(s)
        return '\n'.join(res)

    @classmethod
    def setBands(cls, value: int):
        return Delayed(cls, 0, value)
    

class TwoQuasiFermilevels(CanDelay):
    _delayables = (("nbands", "Number of bands", _checkbands),)

    def __init__(self, carriers: int, valenceBands: int, bands: Union[int,Later] = Later()):
        super().__init__(bands)
        assert type(carriers) is int and carriers > 0, "Number of carriers must be a positive integer"
        assert type(valenceBands) is int and carriers >= 0, "Number of valenceBands must be a positive (or null) integer"
        self._c = carriers
        self._v = valenceBands

    def stamp(self, index: int):
        s = index or ''
        r = f"occopt{s} 9\nivalence{s} {self._v}\nnqfd{s} {self._c}"
        d = super().stamp(index)
        if len(d) > 0:
            r += d
        return r

    @classmethod
    def setBands(cls, value: int):
        return Delayed(cls, 0, value)
    

class OccupationPerBand(Stampable):
    def __init__(self, *occupations: Union[float, Tuple[float, float]], repeat: Optional[int] = None) -> None:
        super().__init__()
        assert all(
            type(v) is float or (type(v) is tuple and len(v) == 2 and type(v[0]) is float and type(v[1]) is float)
            for v in occupations
        ), "Occupations must be float or tuple of two floats"

        self._d = False
        self._o = occupations
        if repeat is None:
            self._r = 0
        elif type(repeat) is int:
            assert len(occupations) == 1, "If repeat is given, only one occupation must be given"
            self._r = repeat
        else:
            raise TypeError("OccupationPerBands 'repeat' argument must be None or positive int")
        
    def compatible(self, coll: StampCollection):
        spin = coll.get(SpinPolarization)
        if spin is None or spin.polarizationNumber == 1:
            self._d = False
            assert all(type(v) is float for v in self._o), "Since the spin is not polarized, ony one occupation per band must be given"
        else:
            self._d = True
    
    def stamp(self, index: int):
        o = ''
        b = 1
        if self._r == 0:
            b = len(self._o)
            if self._d:
                s1 = [0.0]*b
                s2 = [0.0]*b
                for i,o in enumerate(self._o):
                    if type(o) is tuple:
                        s1[i] = o[0]
                        s2[i] = o[1]
                    else:
                        s1[i] = s2[i] = o # type: ignore
                o = ' '.join(str(v) for v in s1) + ' ' + ' '.join(str(v) for v in s2)
            else:
                o = ' '.join(str(o) for o in self._o)
        else:
            b = self._r
            if self._d:
                if type(self._o[0]) is tuple:
                    o = f"{self._r}*{self._o[0][0]} {self._r}*{self._o[0][1]}"
                else:
                    o = f"{self._r*2}*{self._o[0]}"
            else:
                o = f"{self._r}*{self._o[0]}"
        s = index or ''
        return f"occopt{s} 0\nocc{s} {o}\nbands{s} {b}"
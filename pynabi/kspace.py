from ._common import Vec3D as _V, Stampable as _Stmp, _pos_int
from enum import Enum as _E
from typing import Dict as _dict, Any as _any, Union as _union, Tuple as _tuple, Iterable as _iter


__all__ = ["BZ", "CriticalPointsOf", "manual", "SymmetricGrid", "UsualKShifts", "path"]


class BZ(_E):
    Irreducible = 1
    Half = 2
    Full = 3
    NoTimeReversal = 4


class CriticalPointsOf(_E):
    """
    Critical points of (some) Brillouin zones\n
    Taken from http://lampx.tugraz.at/~hadley/ss1/bzones/
    """
    CUB = {
        'R': _V(0.5, 0.5, 0.5),
        'X': _V(0.0, 0.5, 0.0),
        'M': _V(0.5, 0.5, 0.0)
    }
    BCC = {
        'H': _V(-0.5, 0.5, 0.5),
        'P': _V.uniform(0.25),
        'N': _V(0.0, 0.5, 0.0)
    }
    FCC = {
        'X': _V(0.0, 0.5, 0.5),
        'L': _V.uniform(0.5),
        'W': _V(0.25, 0.75, 0.5),
        'U': _V(0.25, 0.625, 0.625),
        'K': _V(0.375, 0.75, 0.375)
    }
    HEX = {
        'A': _V(0,0,1/2),
        'K': _V(2/3,1/3,0),
        'H': _V(2/3,1/3,1/2),
        'M': _V(1/2,0,0),
        'L': _V(1/2,0,1/2)
    }
    TET = {
        'X': (0.5,0.0,0.0),
        'M': (0.5,0.5,0.0),
        'Z': (0.0,0.0,0.5),
        'R': (0.5,0.0,0.5),
        'A': _V.uniform(0.5)
    }
    BCT = {
        'X': _V(0.5,0.0,0.0),
        'Z': _V(0.5,0.5,-0.5),
        'N': _V(0.0,0.5,0.0),
        'P': _V.uniform(0.25)
    }
    ORC = {
        'X': _V(0.5,0.0,0.0),
        'Y': _V(0.0,0.5,0.0),
        'Z': _V(0.0,0.0,0.5),
        'T': _V(0.0,0.5,0.5),
        'U': _V(0.5,0.0,0.5),
        'S': _V(0.5,0.5,0.5),
        'R': _V.uniform(0.5)
    }
    ORCC = {
        'Y': _V(0.5,0.5,0.0),
        'Y': _V(-0.5,0.5,0.0),
        'Z': _V(0.0,0.0,0.5),
        'T': _V.uniform(0.5),
        'T': _V(-0.5,0.5,0.5),
        'S': _V(0.0,0.5,0.0),
        'R': _V(0.0,0.5,0.5),
    }


class UsualKShifts(_E):
    Unshifted = (_V.zero(),)
    Default = (_V.uniform(0.5),)
    BCC = (_V.uniform(0.25), _V.uniform(-0.25))
    FCC = (_V.uniform(0.5), _V(0.5,0.0,0.0), _V(0.0,0.5,0.0), _V(0.0,0.0,0.5))
    HEX = (_V(1.0,0.0,0.0), _V(-0.5,0.8660254037844386,0.0), _V(0.0,0.0,1.0))


class KSpaceDefinition(_Stmp):
    def __init__(self, option: int, props: _dict[str, _any]) -> None:
        super().__init__()
        self.kptop = option
        self.props = props
    
    def stamp(self, index: int):
        s = index or ''
        head = '' if self.kptop == 100 else f"kptopt{s} {self.kptop}\n"
        body = '\n'.join(f"{k}{s} {v}" for k,v in self.props.items())
        return f"{head}{body}"


def manual(*points: _V, normalize: float = 1.0):
    assert normalize >= 1, "k-points normalization faction cannot be lower than 1"
    return KSpaceDefinition(0, {
        "nkpt": len(points),
        "kpt": "   ".join(str(v) for v in points),
        "kptnrm": normalize
    })


class SymmetricGrid:
    def __init__(self, symmetry: BZ, shifts: _union[_tuple[_V,...], UsualKShifts] = ()):
        self._sy = symmetry
        if type(shifts) is UsualKShifts:
            self._sh = shifts.value
        elif type(shifts) is tuple:
            assert all(type(v) is _V for v in shifts), "K Shifts must be vectors"
            self._sh = shifts
        else:
            raise TypeError("Invalid type of k shifts")
    
    def _u(self, n: str, v: str):
        assert len(self._sh) > 0, "There must be at least one k shift"
        return KSpaceDefinition(self._sy.value, {
            "nshiftk": len(self._sh),
            "shiftk": "   ".join(str(s) for s in self._sh),
            n: v
        })
    
    def ofMonkhorstPack(self, a: int, b: _union[int,None] = None, c: _union[int,None] = None):
        b = a if b is None else b
        c = a if c is None else c
        assert _pos_int(a) and _pos_int(b) and _pos_int(c), "Number of k points in Monkhorst-Pack grid must be a positive integer"
        return self._u("ngkpt", f"{a} {b} {c}")
    
    def fromSuperLattice(self, a: _V, b: _V, c: _V):
        return self._u("kptrlatt", f"{a} {b} {c}")
    

    def automatic(self, length: float = 30.0):
        assert length > 0, "Real space length used for automatic k grid must be positive"
        return KSpaceDefinition(self._sy.value, { "kptrlen": length })


def setMPGridPointsNumber(number: _union[int, _tuple[int, int, int]]):
    if type(number) is tuple:
        assert len(number) == 3, "numbers of k points must be three"
        assert all(_pos_int(v) for v in number), "numbers of k points must be positive integers"
        return KSpaceDefinition(100, { "ngkpt": f"{number[0]} {number[1]} {number[2]}" })
    else:
        assert _pos_int(number), "Number of k points (per primitive) must be a positive integer"
        return KSpaceDefinition(100, { "ngkpt": f"{number} {number} {number}" })


def path(divisions: _union[int,_tuple[int,...]], points: _union[str,_iter[_union[str,_V]]], pointSet: _union[CriticalPointsOf,_dict[str,_V]] = {}):
    s: dict[str,_V] = pointSet.value if type(pointSet) is CriticalPointsOf else pointSet  # type: ignore
    b: list[_V] = []
    for p in points:
        if type(p) is str:
            for c in p:
                if c == 'G':
                    b.append(_V.zero())
                else:
                    assert c in s, f"(critical) point '{c}' is not defined"
                    b.append(s[c])
        elif type(p) is _V:
            b.append(p)
        else:
            raise TypeError(f"Invalid type of k-path point (got {type(p)})")
    assert len(b) > 1, "Number of boundaries must be at least 2 (i.e. one segment)"
    
    if type(divisions) is tuple:
        assert len(divisions) == len(b)-1, f"lenght of division must be equal to number of segments: got {len(divisions)} instead of {len(b)-1}"
        assert all(_pos_int(v) for v in divisions), "Number of division per segment must be all positive"
        return KSpaceDefinition(1-len(b), { 
            "kptbounds": "   ".join(str(v) for v in b), 
            "ndivk": ' '.join(str(d) for d in divisions)
        })
    else:
        assert _pos_int(divisions), "number of division (for smallest segment) must be positive"
        return KSpaceDefinition(1-len(b), { 
            "kptbounds": "   ".join(str(v) for v in b), 
            "ndivsm": divisions
        })
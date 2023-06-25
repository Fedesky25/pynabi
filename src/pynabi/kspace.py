from ._common import Vec3D as _V, Stampable as _Stmp, _pos_int, _pos_num, CanDelay as _CD, Delayed as _D, Later as _L
from enum import Enum as _E
from typing import Dict as _dict, Any as _any, Union as _union, Tuple as _tuple, Iterable as _iter


__all__ = ["BZ", "CriticalPointsOf", "ManualGrid", "SymmetricGrid", "UsualKShifts", "Path"]


class BZ(_E):
    """Brillouin Zone symmetries required to setup a symmetric grid"""
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
        'X': _V(0.5,0.0,0.0),
        'M': _V(0.5,0.5,0.0),
        'Z': _V(0.0,0.0,0.5),
        'R': _V(0.5,0.0,0.5),
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


class ManualGrid(_Stmp):
    def __init__(self, *points: _V, normalize: float = 1.0) -> None:
        assert all(type(v) is _V for v in points), "Points of the manual grid must be Vec3D"
        assert normalize >= 1, "k-points normalization faction cannot be lower than 1"
        self.p = points
        self.n = normalize
    
    def stamp(self, index: int):
        s = index or ''
        kpt = "   ".join(str(v) for v in self.p)
        return f"kptopt{s} 0\nnkpt{s} {len(self.p)}\nkpt{s} {kpt}\nkptnrm{s} {self.n}"


def _parse_shifts(value: _tuple[_V,...]|UsualKShifts) -> _tuple[_V,...]:
    if type(value) is UsualKShifts:
        return value.value
    elif type(value) is tuple:
        assert all(type(v) is _V for v in value), "K Shifts must be vectors"
        return value
    else:
        raise TypeError("Invalid type of k shifts")


class _D_grid_points(_CD.info):
    def sanitize(self, value):
        t = type(value)
        if t is int:
            assert value > 0, f"{self.name} must be positive (if integer)"
            return (value,value,value)
        elif t is tuple:
            assert len(value) == 3 and all(_pos_int(v) for v in value), f"{self.name} must be a tuple of three positive integers"
            return value
        else:
            raise TypeError(f"{self.name} must be either a positive integer or a tuple of 3 positive integers")

    def stamp(self, suffix, value):
        return f"{self.prop}{suffix} {value[0]} {value[1]} {value[2]}"


class _D_super_lattice(_CD.info):
    def sanitize(self, value):
        try:
            assert type(value) is tuple and len(value) == 3
            assert all(type(v) is _V for v in value) 
            return value
        except:
            raise TypeError(f"{self.name} must be three Vec3D")
    
    def stamp(self, suffix, value):
        s = '  '.join(str(v) for v in value)
        return f"{self.prop}{suffix} {s}"


class SymmetricGrid(_CD):
    """Constructs a grid of k-points in the reciprocal space leveraging a symmetry of the Brillouin Zone"""

    _delayables = (
        _D_grid_points("ngkpt", "Number of grid points"),
        _D_super_lattice("kptrlatt", "Super lattice vectors")
    )

    def __init__(self, symmetry: BZ, shifts: _union[_tuple[_V,...], UsualKShifts] = ()):
        super().__init__(_L(),_L())
        assert type(symmetry) is BZ, "Symmetry must one of the entries of BZ"
        self.sym = symmetry
        self.shi = _parse_shifts(shifts)
        self.type = -1
    
    def _doesDelay(self, i: int):
        return self.type == i and super()._doesDelay(i)
    
    def ofMonkhorstPack(self, gridPointsNumber: int|tuple[int,int,int]|_L = _L()):
        self.type = 0
        self._dv = (self._delayables[0].laterOrSanitized(gridPointsNumber), _L())
        return self
    
    def fromSuperLattice(self, a: _V, b: _V, c: _V):
        self.type = 1
        self._dv = (_L(), self._delayables[1].laterOrSanitized((a,b,c)))
        return self
    
    @classmethod
    def setMPgridPointNumber(cls, num: int|tuple[int,int,int]):
        return _D(cls, 0, num)
    
    @classmethod
    def setSuperLatticeVectors(cls, a: _V, b: _V, c: _V):
        return _D(cls, 1, (a,b,c))


class AutomaticGrid(_Stmp):
    """ABINIT will automatically generate a large set of possible k point grids, and select among this set, the grids that give a length of smallest vector larger than the provided lenght.
    
    Note that this procedure can be time-consuming. It is worth doing it once for a given unit cell and set of symmetries, but not use this procedure by default. The best is then to use `AbOut().KPointsSets()`, in order to get a detailed analysis of the set of grids."""

    def __init__(self, symmetry: BZ, length: float = 30.0):
        self.sym = symmetry
        self.len = length
    
    def stamp(self, index: int):
        s = index or ''
        return f"kptopt{s} {self.sym.value}\nkptrlen{s} {self.len}"


def _parse_crit_point(c: str, s: dict[str,_V]):
    if c == 'G':
        return _V.zero()
    else:
        assert c in s, f"(critical) point '{c}' is not defined"
        assert type(s[c]) is _V, f"type of critial point '{c}' is not Vec3D"
        return s[c]


class Path(_Stmp):
    """A path though points in the reciprocal space"""

    def __init__(self, points: list[_V]|tuple[_V], prop: str, val: str) -> None:
        """DO NOT USE this constructor"""
        super().__init__()
        self.points = points
        self.prop = prop
        self.val = val
    
    def stamp(self, index: int):
        s = index or ''
        bounds = '  '.join(str(v) for v in self.points)
        return f"kptopt{s} {1-len(self.points)}\nkptbounds{s} {bounds}\n{self.prop}{s} {self.val}"

    @staticmethod
    def auto(minDivisions: int, points: str|_iter[_union[str,_V]], pointSet: CriticalPointsOf|_dict[str,_V] = {}):
        """Path through k points. The number of division for each segment is scaled based on the length and `minDivisions`, the number of division of the smallest segment.
        
        ## Example
        ```python
        # this path
        p1 = Path.auto(10, (
            V.zero(),
            V(0.0, 0.5, 0.5),
            V.uniform(0.5),
            V(0.25, 0.75, 0.5)
        ))
        # is equivalent to
        p2 = Path.auto(10, "GXLW", CriticalPointsOf.FCC)

        # You don't have to sacrifice the comfort of strings 
        # to pass through non critical points
        p3 = Path.auto(10, ("GX", V(0.25,0.5,0.4), "LW"), CriticalPointsOf.FCC)

        # For some application, it could be useful to define
        # custom critical points and use them
        ccp = {
            'A': V(0.25,0.5,0.75),
            'B': V(0.75,0.5,0.25),
            'C': V(-0.5,0.25,0.4)
        }
        p4 = Path.auto(10, "GABGC", ccp) # note that 'G' is always (0,0,0)
        ```"""
        assert _pos_int(minDivisions), "Smallest division must be a positive integer"
        s: dict[str,_V] = pointSet.value if type(pointSet) is CriticalPointsOf else pointSet  # type: ignore
        b: list[_V] = []
        for p in points:
            if type(p) is str:
                for c in p:
                    b.append(_parse_crit_point(c,s))
            elif type(p) is _V:
                b.append(p)
            else:
                raise TypeError(f"Invalid type of k-path point (got {type(p)})")
        assert len(b) > 1, "Number of boundaries must be at least 2 (i.e. one segment)"
        return Path(b, "ndivsm", str(minDivisions))
    
    @staticmethod
    def manual(*args: int|_V|str, pointSet: CriticalPointsOf|_dict[str,_V] = {}):
        """Path though k points where each segment has its own number of divisions. To specify points and divisions, the arguments must be a sequence of alternating positive integers and k points, starting and eending in a k point.

        ## Example
        ```python
        # this path
        p1 = Path.manual( V.zero(), 10, V(0.0, 0.5, 0.5), 15, V.uniform(0.5), 20, V(0.25, 0.75, 0.5))
        # is equivalent to
        p2 = Path.auto('G',10,'X',15,'L',20,'W', pointSet=CriticalPointsOf.FCC)

        # You don't have to sacrifice the comfort of strings 
        # to pass through non critical points
        p3 = Path.auto('G',10,'X',18,V(0.25,0.5,0.4),12,'L',15,'W', pointSet=CriticalPointsOf.FCC)

        # For some application, it could be useful to define
        # custom critical points and use them
        ccp = {
            'A': V(0.25,0.5,0.75),
            'B': V(0.75,0.5,0.25),
            'C': V(-0.5,0.25,0.4)
        }
        p4 = Path.auto('G',7,'A',8,'B',16'G',10,'C', pointSet=ccp) # note that 'G' is always (0,0,0)
        ```"""
        s: dict[str,_V] = pointSet.value if type(pointSet) is CriticalPointsOf else pointSet  # type: ignore
        p: list[_V] = []
        d: list[int] = []
        assert len(args) & 1, "Invalid path sequence"
        for i in range(0, len(args)-1, 2):
            t = type(args[i])
            if t is str:
                p.append(_parse_crit_point(args[i],s)) # type: ignore
            elif t is _V:
                p.append(args[i]) # type: ignore
            else:
                raise TypeError(f"Element number {i+1} must be a vector or a critical point name")
            assert _pos_int(args[i+1]), f"Number of divisions must be a positive integer (at position {i+2})"
            d.append(args[i+1]) # type: ignore
        last = args[len(args)-1]
        t = type(last)
        if t is str:
            p.append(_parse_crit_point(last,s)) # type: ignore
        elif t is _V:
            p.append(last) # type: ignore
        else:
            raise TypeError(f"Last element must be a vector or a critical point name")
        return Path(p, "ndivk", ' '.join(str(v) for v in d))


_ex = (ManualGrid, SymmetricGrid, AutomaticGrid, Path)
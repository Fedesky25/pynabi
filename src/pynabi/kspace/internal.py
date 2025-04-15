"""
WARNING: do not import this file directly!
"""

from pynabi._common import Vec3D as Vec3D, Stampable as Stampable, _pos_int, CanDelay as CanDelay, Delayed as Delayed, Later as Later
from typing import Dict as Dict, Union as Union, Tuple as Tuple, Iterable as Iterable
from enum import Enum as Enum


class BrillouinZone(Enum):
    """Brillouin Zone symmetries required to setup a symmetric grid"""
    Irreducible = 1
    Half = 2
    Full = 3
    NoTimeReversal = 4


class CriticalPointsOf(Enum):
    """
    Critical points of (some) Brillouin zones\n
    Taken from http://lampx.tugraz.at/~hadley/ss1/bzones/
    """
    CUB = {
        'R': Vec3D(0.5, 0.5, 0.5),
        'X': Vec3D(0.0, 0.5, 0.0),
        'M': Vec3D(0.5, 0.5, 0.0)
    }
    BCC = {
        'H': Vec3D(-0.5, 0.5, 0.5),
        'P': Vec3D.uniform(0.25),
        'N': Vec3D(0.0, 0.5, 0.0)
    }
    FCC = {
        'X': Vec3D(0.0, 0.5, 0.5),
        'L': Vec3D.uniform(0.5),
        'W': Vec3D(0.25, 0.75, 0.5),
        'U': Vec3D(0.25, 0.625, 0.625),
        'K': Vec3D(0.375, 0.75, 0.375)
    }
    HEX = {
        'A': Vec3D(0,0,1/2),
        'K': Vec3D(2/3,1/3,0),
        'H': Vec3D(2/3,1/3,1/2),
        'M': Vec3D(1/2,0,0),
        'L': Vec3D(1/2,0,1/2)
    }
    TET = {
        'X': Vec3D(0.5,0.0,0.0),
        'M': Vec3D(0.5,0.5,0.0),
        'Z': Vec3D(0.0,0.0,0.5),
        'R': Vec3D(0.5,0.0,0.5),
        'A': Vec3D.uniform(0.5)
    }
    BCT = {
        'X': Vec3D(0.5,0.0,0.0),
        'Z': Vec3D(0.5,0.5,-0.5),
        'N': Vec3D(0.0,0.5,0.0),
        'P': Vec3D.uniform(0.25)
    }
    ORC = {
        'X': Vec3D(0.5,0.0,0.0),
        'Y': Vec3D(0.0,0.5,0.0),
        'Z': Vec3D(0.0,0.0,0.5),
        'T': Vec3D(0.0,0.5,0.5),
        'U': Vec3D(0.5,0.0,0.5),
        'S': Vec3D(0.5,0.5,0.5),
        'R': Vec3D.uniform(0.5)
    }
    ORCC = {
        'Y': Vec3D(0.5,0.5,0.0),
        'Y': Vec3D(-0.5,0.5,0.0),
        'Z': Vec3D(0.0,0.0,0.5),
        'T': Vec3D.uniform(0.5),
        'T': Vec3D(-0.5,0.5,0.5),
        'S': Vec3D(0.0,0.5,0.0),
        'R': Vec3D(0.0,0.5,0.5),
    }


class UsualKShifts(Enum):
    Unshifted = (Vec3D.zero(),)
    Default = (Vec3D.uniform(0.5),)
    BCC = (Vec3D.uniform(0.25), Vec3D.uniform(-0.25))
    FCC = (Vec3D.uniform(0.5), Vec3D(0.5,0.0,0.0), Vec3D(0.0,0.5,0.0), Vec3D(0.0,0.0,0.5))
    HEX = (Vec3D(1.0,0.0,0.0), Vec3D(-0.5,0.8660254037844386,0.0), Vec3D(0.0,0.0,1.0))


class ManualGrid(Stampable):
    def __init__(self, *points: Vec3D, normalize: float = 1.0) -> None:
        assert all(type(v) is Vec3D for v in points), "Points of the manual grid must be Vec3D"
        assert normalize >= 1, "k-points normalization faction cannot be lower than 1"
        self.p = points
        self.n = normalize
    
    def stamp(self, index: int):
        s = index or ''
        kpt = "   ".join(str(v) for v in self.p)
        return f"kptopt{s} 0\nnkpt{s} {len(self.p)}\nkpt{s} {kpt}\nkptnrm{s} {self.n}"


def _parse_shifts(value: Tuple[Vec3D,...]|UsualKShifts) -> Tuple[Vec3D,...]:
    if type(value) is UsualKShifts:
        return value.value
    elif type(value) is tuple:
        assert all(type(v) is Vec3D for v in value), "K Shifts must be vectors"
        return value
    else:
        raise TypeError("Invalid type of k shifts")


class _D_grid_points(CanDelay.info):
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


class _D_super_lattice(CanDelay.info):
    def sanitize(self, value):
        try:
            assert type(value) is tuple and len(value) == 3
            assert all(type(v) is Vec3D for v in value) 
            return value
        except:
            raise TypeError(f"{self.name} must be three Vec3D")
    
    def stamp(self, suffix, value):
        s = '  '.join(str(v) for v in value)
        return f"{self.prop}{suffix} {s}"


class SymmetricGrid(CanDelay):
    """Constructs a grid of k-points in the reciprocal space leveraging a symmetry of the Brillouin Zone"""

    _delayables = (
        _D_grid_points("ngkpt", "Number of grid points"),
        _D_super_lattice("kptrlatt", "Super lattice vectors")
    )

    def __init__(self, symmetry: BrillouinZone, shifts: Union[Tuple[Vec3D,...], UsualKShifts] = ()):
        super().__init__(Later(),Later())
        assert type(symmetry) is BrillouinZone, "Symmetry must one of the entries of BZ"
        self.sym = symmetry
        self.shi = _parse_shifts(shifts)
        self.type = -1
    
    def _doesDelay(self, i: int):
        return self.type == i and super()._doesDelay(i)
    
    def ofMonkhorstPack(self, gridPointsNumber: int|tuple[int,int,int]|Later = Later()):
        assert self.type == -1, "Symmetric grid type redefined"
        self.type = 0
        self._dv = (self._delayables[0].laterOrSanitized(gridPointsNumber), Later())
        return self
    
    def fromSuperLattice(self, a: Vec3D, b: Vec3D, c: Vec3D):
        assert self.type == -1, "Symmetric grid type redefined"
        self.type = 1
        self._dv = (Later(), self._delayables[1].laterOrSanitized((a,b,c)))
        return self
    
    def stamp(self, index: int):
        s = index or ''
        shift = " ".join(str(t) for t in self.shi)
        return f"kptopt{s} {self.sym.value}\nnshiftk{s} {len(self.shi)}\nshiftk{s} {shift}\n{super().stamp(index)}"
    
    @classmethod
    def setMPgridPointNumber(cls, num: int|tuple[int,int,int]):
        """Sets the number of k points in the Monkhorst-Pack grid"""
        return Delayed(cls, 0, num)
    
    @classmethod
    def setSuperLatticeVectors(cls, a: Vec3D, b: Vec3D, c: Vec3D):
        """Sets the vectors of the super lattice"""
        return Delayed(cls, 1, (a,b,c))


class AutomaticGrid(Stampable):
    """ABINIT will automatically generate a large set of possible k point grids, and select among this set, the grids that give a length of smallest vector larger than the provided lenght.
    
    Note that this procedure can be time-consuming. It is worth doing it once for a given unit cell and set of symmetries, but not use this procedure by default. The best is then to use `AbOut().KPointsSets()`, in order to get a detailed analysis of the set of grids."""

    def __init__(self, symmetry: BrillouinZone, length: float = 30.0):
        self.sym = symmetry
        self.len = length
    
    def stamp(self, index: int):
        s = index or ''
        return f"kptopt{s} {self.sym.value}\nkptrlen{s} {self.len}"


def _parse_crit_point(c: str, s: dict[str,Vec3D]):
    if c == 'G':
        return Vec3D.zero()
    else:
        assert c in s, f"(critical) point '{c}' is not defined"
        assert type(s[c]) is Vec3D, f"type of critial point '{c}' is not Vec3D"
        return s[c]


class Path(Stampable):
    """A path though points in the reciprocal space"""

    def __init__(self, points: list[Vec3D]|tuple[Vec3D], prop: str, val: str) -> None:
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
    def auto(minDivisions: int, points: str|Iterable[Union[str,Vec3D]], pointSet: CriticalPointsOf|Dict[str,Vec3D] = {}):
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
        s: dict[str,Vec3D] = pointSet.value if type(pointSet) is CriticalPointsOf else pointSet  # type: ignore
        b: list[Vec3D] = []
        for p in points:
            if type(p) is str:
                for c in p:
                    b.append(_parse_crit_point(c,s))
            elif type(p) is Vec3D:
                b.append(p)
            else:
                raise TypeError(f"Invalid type of k-path point (got {type(p)})")
        assert len(b) > 1, "Number of boundaries must be at least 2 (i.e. one segment)"
        return Path(b, "ndivsm", str(minDivisions))
    
    @staticmethod
    def manual(*args: int|Vec3D|str, pointSet: CriticalPointsOf|Dict[str,Vec3D] = {}):
        """Path though k points where each segment has its own number of divisions. To specify points and divisions, the arguments must be a sequence of alternating positive integers and k points, starting and ending in a k point.

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
        p4 = Path.auto('G',7,'A',8,'B',16,'G',10,'C', pointSet=ccp) # note that 'G' is always (0,0,0)
        ```"""
        s: dict[str,Vec3D] = pointSet.value if type(pointSet) is CriticalPointsOf else pointSet  # type: ignore
        p: list[Vec3D] = []
        d: list[int] = []
        assert len(args) & 1, "Invalid path sequence"
        for i in range(0, len(args)-1, 2):
            t = type(args[i])
            if t is str:
                p.append(_parse_crit_point(args[i],s)) # type: ignore
            elif t is Vec3D:
                p.append(args[i]) # type: ignore
            else:
                raise TypeError(f"Element number {i+1} must be a vector or a critical point name")
            assert _pos_int(args[i+1]), f"Number of divisions must be a positive integer (at position {i+2})"
            d.append(args[i+1]) # type: ignore
        last = args[len(args)-1]
        t = type(last)
        if t is str:
            p.append(_parse_crit_point(last,s)) # type: ignore
        elif t is Vec3D:
            p.append(last) # type: ignore
        else:
            raise TypeError(f"Last element must be a vector or a critical point name")
        return Path(p, "ndivk", ' '.join(str(v) for v in d))


_exclusives = (ManualGrid, SymmetricGrid, AutomaticGrid, Path)
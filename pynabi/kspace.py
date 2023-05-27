from ._common import Vec3D as _V, Stampable as _Stmp
from enum import Enum as _E
from typing import Dict as _dict, Any as _any, Union as _union, Tuple as _tuple, Iterable as _iter


__all__ = ["BZ", "CriticalPointsOf", "manual", "MonkhorstPackGrid", "fromSuperLattice"]


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


class KSpaceDefinition(_Stmp):
    def __init__(self, option: int, props: _dict[str, _any]) -> None:
        super().__init__()
        self.kptop = option
        self.props = props
    
    def stamp(self, index: int):
        s = index or ''
        body = '\n'.join(f"{k}{s} {v}" for k,v in self.props.items())
        return f"kptopt{s} {self.kptop}\n{body}"


def manual(*points: _V, normalize: float = 1.0):
    assert normalize >= 1, "k-points normalization faction cannot be lower than 1"
    return KSpaceDefinition(0, {
        "nkpt": len(points),
        "kpt": "   ".join(str(v) for v in points),
        "kptnrm": normalize
    })


def _pos_int(v):
    return type(v) is int and v > 0


def MonkhorstPackGrid(symmetry: BZ, number: _union[int, _tuple[int, int, int]], *shifts: _V):
    """Creates a symmetric k grid
    
    | name | | description |
    | ---: | - | :---------- |
    | symmetry | - | Brillouin zone symmetry to use |
    | number | - | Number of k points of Monkhorst-Pack grids |
    | shiftk | - | Shifts of the homogeneous grid of k points |
    """
    assert len(shifts) > 0, "There must be at least one k shift"
    if type(number) is tuple:
        assert len(number) == 3, "numbers of k points must be three"
        assert all(_pos_int(v) for v in number), "numbers of k points must be positive integers"
    else:
        assert _pos_int(number), "Number of k points (per primitive) must be a positive integer"
        number = (number, number, number) # type: ignore
    
    return KSpaceDefinition(symmetry.value, {
        "ngkpt": ' '.join(number), # type: ignore
        "nshiftk": len(shifts),
        "shiftk": "   ".join(str(s) for s in shifts)
    })


def fromSuperLattice(symmetry: BZ, coordinates: _tuple[_V, _V, _V], *shifts: _V):
    assert len(shifts) > 0, "There must be at least one k shift"
    assert len(coordinates) == 3 and all(type(v) is _V for v in coordinates), "There must be three super lattice vectors"
    return KSpaceDefinition(symmetry.value, {
        "kptrlatt": "   ".join(str(v) for v in coordinates),
        "nshiftk": len(shifts),
        "shiftk": "   ".join(str(s) for s in shifts)
    })


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
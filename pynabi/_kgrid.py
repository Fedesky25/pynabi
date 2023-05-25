from ._common import Vec3D, Stampable
from typing import Union, Tuple
from enum import Enum

__all__ = ["BZ", "ManualKGrid", "SymmetricKGrid", "KPath"]

class KGrid(Stampable):
    pass
    
class ManualKGrid(KGrid):
    def __init__(self, *points: Vec3D, normalize: float = 1) -> None:
        super().__init__()
        self.points = points
        self.norm = normalize
        assert normalize >= 1, "k-points normalization faction cannot be lower than 1"
    
    def stamp(self, index: int):
        suf = index or ''
        kpt_name = f"kpt{suf}"
        kpt_value = ('\n ' + ' '*len(kpt_name)).join(str(v) for v in self.points)
        return f"""kptopt{suf} 0
nkpt{suf} {len(self.points)}
{kpt_name} {kpt_value}
kptnrm{suf} {self.norm}"""


class BZ(Enum):
    Irreducible = 1
    Half = 2
    Full = 3
    NoTimeReversal = 4


class SymmetricKGrid(KGrid):
    def __init__(self, symmetry: BZ, number: Union[int, Tuple[int,int,int]], *shifts: Vec3D) -> None:
        super().__init__()
        self.symmetry = symmetry.value
        self.shifts = shifts
        n: Tuple[int,int,int] = (number, number, number) if type(number) is int else number # type: ignore
        assert n[0] > 0 and n[1] > 0 and n[2] > 0, "number of k points must be positive"
        self.number = n
    
    def stamp(self, index: int):
        s = index or ''
        shift_name = f"kpt{s}"
        shift_value = ('\n ' + ' '*len(shift_name)).join(str(v) for v in self.shifts)
        return f"""kptopt{s} {self.symmetry}
ngkpt{s} {self.number[0]} {self.number[1]} {self.number[2]}
nshiftk{s} {len(self.shifts)}
{shift_name} {shift_value}"""


class KPath(KGrid):
    def __init__(self, divisions: Union[int,Tuple[int,...]], *boundaries: Vec3D) -> None:
        super().__init__()
        self.divisions = divisions
        self.boundaries = boundaries
        assert len(boundaries) >= 2, "k-path must have at least 2 segment boundaries"
        if type(divisions) is int:
            assert divisions > 0, "number of division (for smallest segment must be positive)"
        elif type(divisions) is tuple:
            assert len(divisions) == len(boundaries) - 1, f"lenght of division must be equal to number of segments: got {len(divisions)} instead of {len(boundaries)-1}"
            assert all(n > 0 for n in divisions), "Number of division per segment must be all positive"
        else:
            raise TypeError("Divisions for KPath must be either integer of tuple of integers")
    
    def stamp(self, index: int):
        s = index or ''
        bn = f"kptbounds{s}"
        bv = ('\n '+' '*len(bn)).join(str(b) for b in self.boundaries)
        dn, dv = ("ndivsm", self.divisions) if type(self.divisions) is int else ("ndivk", ' '.join(str(d) for d in self.divisions)) # type: ignore
        return f"""kptopt{s} -{len(self.boundaries)-1}
{bn} {bv}
{dn}{s} {dv}"""

from .occupation import Occupation as _OCC, SpinPolarization as _SP, SpinType as _ST
from .kspace import manual as _MKG
from ._common import Vec3D as _Vec3D, SKO as _SKO

__all__ = ["uniformBandNumber", "spinUnpolarizedVariableBandNumber"]


def uniformBandNumber(spin: _SP, occupation: float, bands: int, *points: _Vec3D, normalize: float = 1):
    bands = len(points)
    occ = _OCC.EqualBandNumber(occupation, bands, spin.polarizationNumber == 2)
    kgrid = _MKG(*points, normalize=normalize)
    return (kgrid, occ, spin)


def spinUnpolarizedVariableBandNumber(spin: _SP, *points: _SKO, normalize: float = 1):
    assert spin.polarizationNumber != 2, "Spin must be unpolarized in spinUnpolarizedVariableBandNumber"
    kgrid = _MKG(*(p.k for p in points), normalize=normalize)
    bands = tuple(len(p.occ) for p in points)
    occ_str = " ".join(" ".join(str(o) for o in p.occ) for p in points)
    return (spin, kgrid, _OCC(2, bands, occ=occ_str))

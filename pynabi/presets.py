from ._crystal import Lattice as _Lattice, AtomBasis as _AtomBasis, Atom as _Atom
from ._occupation import Occupation as _OCC, SpinPolarization as _SP, SpinType as _ST
from .kspace import manual as _MKG
from ._common import Vec3D as _Vec3D, SKO as _SKO
from typing import Union

__all__ = ["cubic", "FCC", "BCC", "HCP", "uniformBandNumber", "spinUnpolarizedVariableBandNumber"]

def cubic(a: float, atom: _Atom, pseudosAt: str = ''):
    return (
        _AtomBasis(
            (atom, _Vec3D.zero()),
            pseudosAt=pseudosAt
        ),
        _Lattice.fromAngles(
            _Vec3D.uniform(90),
            _Vec3D.uniform(a)
        )
    )

def BCC(a: float, atom: _Atom, center: Union[_Atom,None] = None, pseudosAt: str = ''):
    r = cubic(a, atom, pseudosAt)
    r[0].add(atom if center is None else center, _Vec3D.uniform(0.5))
    return r

def FCC(a: float, atom: _Atom, face: Union[_Atom,None] = None, pseudosAt: str = ''):
    r = cubic(a, atom, pseudosAt)
    f = atom if face is None else face
    r[0].add(f, _Vec3D(0.5,0.5,0.0))
    r[0].add(f, _Vec3D(0.0,0.5,0.5))
    r[0].add(f, _Vec3D(0.5,0.0,0.5))
    return r

def ZincBlende(a: float, atomA: _Atom, atomB: _Atom, pseudosAt: str = ''):
    return (
        _AtomBasis(
            (atomA, _Vec3D.zero()),
            (atomA, _Vec3D(0.5,0.5,0.0)),
            (atomA, _Vec3D(0.0,0.5,0.5)),
            (atomA, _Vec3D(0.5,0.0,0.5)),
            (atomB, _Vec3D.uniform(0.25)),
            (atomB, _Vec3D(0.75, 0.75, 0.25)),
            (atomB, _Vec3D(0.75, 0.25, 0.75)),
            (atomB, _Vec3D(0.25, 0.75, 0.75)),
            pseudosAt=pseudosAt
        ),
        _Lattice.fromAngles(
            _Vec3D.uniform(90),
            _Vec3D.uniform(a)
        )
    )

def HCP(a: float, c: float, atomA: _Atom, atomB: _Atom, pseudosAt: str = ''):
    return (
        _AtomBasis(
            (atomA, _Vec3D.zero()),
            (atomB, _Vec3D(2/3, 1/3, 0.5)),
            pseudosAt=pseudosAt
        ),
        _Lattice.fromAngles(
            _Vec3D(90,120,90),
            _Vec3D(a,a,c)
        )
    )


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

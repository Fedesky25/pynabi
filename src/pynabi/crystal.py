from ._common import Vec3D as _V, Stampable as _S, _pos_int
from typing import Optional as _O, Union as _U
from .units import Length as _L, Pos3D as _P

__all__ = ["Atom", "AtomBasis", "Lattice", "HCP", "CsClLike", "RockSaltLike", "FluoriteLike", "ZincBlendeLike", "HCP", "WurtziteLike", "NiAsLike"]

_atom_symbols = ["H","He","Li","Be","B","C","N","O","F","Ne","Na","Mg","Al","Si","P","S","Cl","Ar","K","Ca","Sc","Ti","V","Cr","Mn","Fe","Co","Ni","Cu","Zn","Ga","Ge","As","Se","Br","Kr","Rb","Sr","Y","Zr","Nb","Mo","Tc","Ru","Rh","Pd","Ag","Cd","In","Sn","Sb","Te","I","Xe","Cs","Ba","La","Ce","Pr","Nd","Pm","Sm","Eu","Gd","Tb","Dy","Ho","Er","Tm","Yb","Lu","Hf","Ta","W","Re","Os","Ir","Pt","Au","Hg","Tl","Pb","Bi","Po","At","Rn","Fr","Ra","Ac","Th","Pa","U ","Np","Pu","Am","Cm","Bk","Cf","Es","F","Md","No","Lr","Rf","Db","Sg","Bh","Hs","Mt","Ds","Rg","Cn","Nh","Fl","Mc","Lv","Ts","Og"]

class Atom:
    def __init__(self, Z: int, file: str):
        assert _pos_int(Z), "Atomic number must be a positive integer"
        self.num = Z
        self.file = file
        assert self.file.find(",") == -1, "File name cannot contain ','"
    
    def __str__(self) -> str:
        return f"{_atom_symbols[self.num-1]} atom (located at {self.file})"
    
    # def __eq__(self, __value: object) -> bool:
    #     return (self is __value) or (type(__value) is Atom and __value.num == self.num and __value.file == self.file)

    @staticmethod
    def of(name: str):
        """Constructs an atom using default filenaming
        
        returns Atom(Z of element, element name + '.psp8')
        """
        try:
            Z = _atom_symbols.index(name) + 1;
            return Atom(Z, name + '.psp8')
        except:
            raise ValueError(name + " is not a valid atom type")
    
    @staticmethod
    def poolstr(atoms: 'list[Atom]'):
        return f"""# Atoms definition
ntypat {len(atoms)}
znucl {' '.join(str(a.num) for a in atoms)}
pseudos \"{', '.join(a.file for a in atoms)}\""""


class AtomBasis(_S):
    def __init__(self, *atoms: 'tuple[Atom,_V]', cartesian: bool = False) -> None:
        """Construct an atom basis given a sequence of tuples containing an atom and a position (Vec3D)\n
        If `cartesian` is True, the coordinates of the atoms' position are cartesian instead of reduced."""
        assert len(atoms) > 0, "There must be at least one atom in basis"
        self.atoms = atoms
        self.cartesian = cartesian
    
    def add(self, atom: Atom, where: _V):
        if type(self.atoms) is tuple:
            self.atoms = list(self.atoms)
        self.atoms.append((atom, where)) # type: ignore
    
    def getAtoms(self):
        return (a[0] for a in self.atoms)
    
    def stamp(self, index: int, pool: 'list[Atom]'):
        indexes = [pool.index(a[0]) for a in self.atoms]
        suffix = str(index or '');
        x_type = "xcart" if self.cartesian else "xred";
        return f"""natom{suffix} {len(indexes)}
typeat{suffix} {' '.join(str(i+1) for i in indexes)}
{x_type}{suffix} {'   '.join(str(a[1]) for a in self.atoms)}"""

    @staticmethod
    def ofOne(atom: Atom):
        """Construct an atom basis with only one atom (placed at Vec3D(0,0,0))"""
        return AtomBasis((atom, _V.zero()))


def _2uni(v: _U[float,_L]):
    t = type(v)
    if t is _L:
        return _P.uniform(v) # type: ignore
    elif t is float:
        return _V.uniform(v) # type: ignore
    else:
        raise TypeError(f"Lattice constant is of wrong value (got {t} instead of float or Length)")


class Lattice(_S):
    def __init__(self, **props):
        """Do not use directly: prefer static methods like fromAngle, fromPrimitives"""
        self._p = props

    def stamp(self, index: int):
        suffix = index if index > 0 else ''
        return '\n'.join(f"{k}{suffix} {v}" for k,v in self._p.items())
    
    @staticmethod
    def fromAngles(angles: _V, scaling: _U[_V,_P]):
        """
        Constructs lattice from angles [&alpha;, &beta;, &gamma;]\n
        Scaling is defined per primitive vector (e.g. a is multiplied by scaling[0])

        Angles are defined between vectors as follows:
        | angle | from | to |
        | :---: | :----------: | :-----------: |
        | &alpha;(1) | b(2) | c(3) |
        | &beta;(2) | a(1) | c(3) |
        | &gamma;(3) | a(1) | b(2) |
        """
        assert type(angles) is _V, "Angles must be of type Vec3D"
        return Lattice(acell=_P.sanitize(scaling), angdeg=angles)
    
    @staticmethod
    def fromPrimitives(a: _V, b: _V, c: _V, scaling: _U[_V,_P]):
        """
        Construct lattice from primitives [a, b, c]\n
        Scaling is applied to the cartesian axis
        """
        assert type(a) is _V and type(b) is _V and type(c) is _V, "Primitive vectors must be of type Vec3D"
        return Lattice(scalecart=_P.sanitize(scaling), rprim=f"{a}   {b}   {c}")
    
    @staticmethod
    def CUB(a: _U[float,_L]):
        return Lattice.fromAngles(_V.uniform(90),_2uni(a))

    @staticmethod
    def BCC(a: _U[float,_L]):
        return Lattice.fromPrimitives(
            _V(-0.5,0.5,0.5),
            _V(0.5,-0.5,0.5),
            _V(0.5,0.5,-0.5),
            _2uni(a)
        )

    @staticmethod
    def FCC(a: _U[float,_L]):
        """Rhombohedral primitive cell of the FCC lattice"""
        return Lattice.fromPrimitives(
            _V(0.5, 0.5, 0.0),
            _V(0.0, 0.5, 0.5),
            _V(0.5, 0.0, 0.5),
            _2uni(a)
        )

    @staticmethod
    def HEX(a: float, c: float, unit: _O[_L] = None):
        return Lattice.fromPrimitives(
            _V(-1.0,0.0,0.0),
            _V(-0.5,0.8660254037844386,0.0),
            _V(0.0,0.0,1.0),
            _P(a,a,c,unit)
        )

    @staticmethod
    def TET(a: float, c: float, unit: _O[_L] = None):
        """Simple tetragonal"""
        return Lattice.fromAngles(
            _V.uniform(90),
            _P(a,a,c,unit)
        )
    
    @staticmethod
    def BCT(a: float, c: float, unit: _O[_L] = None):
        """Body-centered tetragonal"""
        return Lattice.fromPrimitives(
            _V(-0.5,0.5,0.5),
            _V(0.5,-0.5,0.5),
            _V(0.5,0.5,-0.5),
            _P(a,a,c,unit)
        )

    @staticmethod
    def ORC(dimensions: _U[_V,_P]):
        """Simple orthorhombic"""
        return Lattice.fromAngles(_V.uniform(90),dimensions)
    
    @staticmethod
    def ORCC(dimensions: _U[_V,_P]):
        """Base-centered orthorhombic"""
        return Lattice.fromPrimitives(
            _V(0.5,-0.5,0.0),
            _V(0.5,0.5,0.0),
            _V(0.0,0.0,1.0),
            dimensions
        )


def CsClLike(atomA: Atom, atomB: Atom, a: _U[float,_L]):
    """Two interpenetrating primitive cubic structure"""
    l = Lattice.CUB(a)
    b = AtomBasis(
        (atomA, _V.zero()),
        (atomB, _V.uniform(0.5))
    )
    return (b,l)


def RockSaltLike(atomA: Atom, atomB: Atom, a: _U[float,_L]):
    """
    Two interpenetrating FCC that form a chessboard like crystal\n
    Examples: NaCl, PbS
    """
    l = Lattice.FCC(a)
    b = AtomBasis(
        (atomA, _V.zero()),
        (atomB, _V.uniform(0.5))
    )
    return (b,l)


def FluoriteLike(Ca: Atom, F: Atom, a: _U[float,_L]):
    """
    FCC crystal with stechiometric ration of 1:2 between atoms\n
    Examples: BaF2, β-PbF2, PuO2, SrF2, UO2, CaF2, ZrO2, K2O , K2S , Li2O, Na2O, Na2S, Rb2O, Mg2Si
    """
    l = Lattice.FCC(a)
    b = AtomBasis(
        (Ca, _V.zero()),
        (F, _V.uniform(1/3)),
        (F, _V.uniform(2/3)),
    )
    return (b,l)


def ZincBlendeLike(atomA: Atom, atomB: Atom, a: _U[float,_L]):
    """
    Two interpenetrating FCC crystal\n
    Examples:
     - Diamond: both atoms are C
     - Zincblende: basis of Zn & S
    """
    l = Lattice.FCC(a)
    b = AtomBasis(
        (atomA, _V.zero()),
        (atomB, _V.uniform(0.25))
    )
    return (b,l)


def HCP(atomA: Atom, atomB: Atom, a: float, c: float, unit: _O[_L] = None):
    """Hexagonal Close Packet crystal"""
    l = Lattice.HEX(a,c,unit)
    b = AtomBasis(
        (atomA, _V.zero()),
        (atomB, _V(2/3, 1/3, 0.5)),
    )
    return (b,l)


def WurtziteLike(atomA: Atom, atomB: Atom, a: float, c: float, unit: _O[_L] = None):
    """
    Two interpenetrating HCP crystal\n
    Examples: wurtzite (ZnS), silver iodide (AgI), zinc oxide (ZnO), cadmium sulfide (CdS), cadmium selenide (CdSe), silicon carbide (α-SiC), gallium nitride (GaN), aluminium nitride (AlN), boron nitride (w-BN)
    """
    l = Lattice.HEX(a,c,unit)
    b = AtomBasis(
        (atomA, _V.zero()),
        (atomB, _V(2/3, 1/3, 0.25)),
        (atomA, _V(2/3, 1/3, 0.5)),
        (atomB, _V(0.0,0.0,0.75))
    )
    return (b,l)


def NiAsLike(Ni: Atom, As: Atom, a: float, c: float, unit: _O[_L] = None):
    """
    Iterpenetrating HEX and HCP\n
    Examples: 
     - Achavalite: FeSe
     - Breithauptite: NiSb
     - Freboldite: CoSe
     - Kotulskite: Pd(Te,Bi)
     - Langistite: (Co,Ni)As
     - Nickeline: NiAs
     - Sobolevskite: Pd(Bi,Te)
     - Sudburyite: (Pd,Ni)Sb
    """
    l = Lattice.HEX(a,c,unit)
    b = AtomBasis(
        (Ni, _V.zero()),
        (As, _V(2/3, 1/3, 0.25)),
        (Ni, _V(0.0,0.0,0.5)),
        (As, _V(1/3, 2/3, 0.75))
    )
    return (b,l)
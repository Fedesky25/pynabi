"""
WARNING: do not import this file directly!
"""

from pynabi._common import Vec3D, Stampable
from typing import Optional, Union, Tuple
from pynabi.units.internal import Length, Pos3D


atom_symbols = ["H","He","Li","Be","B","C","N","O","F","Ne","Na","Mg","Al","Si","P","S","Cl","Ar","K","Ca","Sc","Ti","V","Cr","Mn","Fe","Co","Ni","Cu","Zn","Ga","Ge","As","Se","Br","Kr","Rb","Sr","Y","Zr","Nb","Mo","Tc","Ru","Rh","Pd","Ag","Cd","In","Sn","Sb","Te","I","Xe","Cs","Ba","La","Ce","Pr","Nd","Pm","Sm","Eu","Gd","Tb","Dy","Ho","Er","Tm","Yb","Lu","Hf","Ta","W","Re","Os","Ir","Pt","Au","Hg","Tl","Pb","Bi","Po","At","Rn","Fr","Ra","Ac","Th","Pa","U ","Np","Pu","Am","Cm","Bk","Cf","Es","F","Md","No","Lr","Rf","Db","Sg","Bh","Hs","Mt","Ds","Rg","Cn","Nh","Fl","Mc","Lv","Ts","Og"]

class Atom:
    def __init__(self, id: int|str, file: Optional[str] = None):
        """Constructs an atom

        `id` can be either the atomic number or symbol 

        `file` is the path the pseudopotential file. If None, it defaults to `f"{atom_symbol}.psp8"`
        """
        if type(id) is int:
            assert 0 < id <= len(atom_symbols), f"Atomic number must be an integer between 1 and {len(atom_symbols)}"
            self.num = id
        elif type(id) is str:
            self.num = atom_symbols.index(id) + 1
        else:
            raise ValueError(f"Invalid atom identifier ({id})")
        if file is None:
            self.file = f"{atom_symbols[self.num-1]}.psp8"
        else:
            self.file = file

    def __hash__(self) -> int:
        return hash(self.file)
    
    def __str__(self) -> str:
        return f"{atom_symbols[self.num-1]} atom (pseudopotential at {self.file})"
    
    def __eq__(self, other: object) -> bool:
        return type(other) is Atom and other.num == self.num and other.file == self.file
    
    @staticmethod
    def poolstr(atoms: 'list[Atom]'):
        return f"""# Atoms definition
ntypat {len(atoms)}
znucl {' '.join(str(a.num) for a in atoms)}
pseudos \"{', '.join(a.file for a in atoms)}\""""


class AtomBasis(Stampable):
    def __init__(self, *atoms: Tuple[Atom,Vec3D], cartesian: bool = False) -> None:
        """Construct an atom basis given a sequence of tuples containing an atom and a position (Vec3D)\n
        If `cartesian` is True, the coordinates of the atoms' position are cartesian instead of reduced."""
        assert len(atoms) > 0, "There must be at least one atom in basis"
        self.atoms = atoms
        self.cartesian = cartesian
    
    def add(self, atom: Atom, where: Vec3D):
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
typat{suffix} {' '.join(str(i+1) for i in indexes)}
{x_type}{suffix} {'   '.join(str(a[1]) for a in self.atoms)}"""

    @staticmethod
    def ofOne(atom: Atom):
        """Construct an atom basis with only one atom (placed at Vec3D(0,0,0))"""
        return AtomBasis((atom, Vec3D.zero()))


def _2uni(v: Union[float,Length]):
    t = type(v)
    if t is Length:
        return Pos3D.uniform(v) # type: ignore
    elif t is float:
        return Vec3D.uniform(v) # type: ignore
    else:
        raise TypeError(f"Lattice constant is of wrong value (got {t} instead of float or Length)")


class Lattice(Stampable):
    def __init__(self, **props):
        """Do not use directly: prefer static methods like fromAngle, fromPrimitives"""
        self._p = props

    def stamp(self, index: int):
        suffix = index if index > 0 else ''
        return '\n'.join(f"{k}{suffix} {v}" for k,v in self._p.items())
    
    @staticmethod
    def fromAngles(angles: Vec3D, scaling: Union[Vec3D,Pos3D]):
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
        assert type(angles) is Vec3D, "Angles must be of type Vec3D"
        return Lattice(acell=Pos3D.sanitize(scaling), angdeg=angles)
    
    @staticmethod
    def fromPrimitives(a: Vec3D, b: Vec3D, c: Vec3D, scaling: Union[Vec3D,Pos3D]):
        """
        Construct lattice from dimensionless primitives [a, b, c]. Each primitive gets scaled by the corresponfing component of `scaling`.
        """
        assert type(a) is Vec3D and type(b) is Vec3D and type(c) is Vec3D, "Primitive vectors must be of type Vec3D"
        return Lattice(acell=Pos3D.sanitize(scaling), rprim=f"{a}   {b}   {c}")
    
    @staticmethod
    def CUB(a: Union[float,Length]):
        return Lattice.fromAngles(Vec3D.uniform(90),_2uni(a))

    @staticmethod
    def BCC(a: Union[float,Length]):
        return Lattice.fromPrimitives(
            Vec3D(-0.5,0.5,0.5),
            Vec3D(0.5,-0.5,0.5),
            Vec3D(0.5,0.5,-0.5),
            Pos3D.uniform(a)
        )

    @staticmethod
    def FCC(a: Union[float,Length]):
        """Rhombohedral primitive cell of the FCC lattice"""
        return Lattice.fromPrimitives(
            Vec3D(0.5, 0.5, 0.0),
            Vec3D(0.0, 0.5, 0.5),
            Vec3D(0.5, 0.0, 0.5),
            Pos3D.uniform(a)
        )

    @staticmethod
    def HEX(a: Union[Length,float], c: Union[Length,float]):
        return Lattice.fromPrimitives(
            Vec3D(-1.0,0.0,0.0),
            Vec3D(-0.5,0.8660254037844386,0.0),
            Vec3D(0.0,0.0,1.0),
            Pos3D(a,a,c)
        )

    @staticmethod
    def TET(a: Union[Length,float], c: Union[Length,float]):
        """Simple tetragonal"""
        return Lattice.fromAngles(
            Vec3D.uniform(90),
            Pos3D(a,a,c)
        )
    
    @staticmethod
    def BCT(a: Union[Length,float], c: Union[Length,float]):
        """Body-centered tetragonal"""
        a = Length.sanitize(a)
        c = Length.sanitize(c)
        r = c / a
        return Lattice.fromPrimitives(
            Vec3D(-0.5,0.5,0.5*r),
            Vec3D(0.5,-0.5,0.5*r),
            Vec3D(0.5,0.5,-0.5*r),
            Pos3D.uniform(a)
        )

    @staticmethod
    def ORC(dimensions: Union[Vec3D,Pos3D]):
        """Simple orthorhombic"""
        return Lattice.fromAngles(Vec3D.uniform(90),dimensions)
    
    @staticmethod
    def ORCC(dimensions: Union[Vec3D,Pos3D]):
        """Base-centered orthorhombic"""
        return Lattice.fromPrimitives(
            Vec3D(0.5*dimensions.x,-0.5*dimensions.y, 0.0*dimensions.z),
            Vec3D(0.5*dimensions.x, 0.5*dimensions.y, 0.0*dimensions.z),
            Vec3D(0.0*dimensions.x, 0.0*dimensions.y, 1.0*dimensions.z),
            Pos3D(1,1,1,Length(1, dimensions.u)) # type: ignore
        )


def CsClLike(atomA: Atom, atomB: Atom, a: Union[float,Length]):
    """Two interpenetrating primitive cubic structure"""
    l = Lattice.CUB(a)
    b = AtomBasis(
        (atomA, Vec3D.zero()),
        (atomB, Vec3D.uniform(0.5))
    )
    return (b,l)


def RockSaltLike(atomA: Atom, atomB: Atom, a: Union[float,Length]):
    """
    Two interpenetrating FCC that form a chessboard like crystal\n
    Examples: NaCl, PbS
    """
    l = Lattice.FCC(a)
    b = AtomBasis(
        (atomA, Vec3D.zero()),
        (atomB, Vec3D.uniform(0.5))
    )
    return (b,l)


def FluoriteLike(Ca: Atom, F: Atom, a: Union[float,Length]):
    """
    FCC crystal with stechiometric ration of 1:2 between atoms\n
    Examples: BaF2, β-PbF2, PuO2, SrF2, UO2, CaF2, ZrO2, K2O , K2S , Li2O, Na2O, Na2S, Rb2O, Mg2Si
    """
    l = Lattice.FCC(a)
    b = AtomBasis(
        (Ca, Vec3D.zero()),
        (F, Vec3D.uniform(1/3)),
        (F, Vec3D.uniform(2/3)),
    )
    return (b,l)


def ZincBlendeLike(atomA: Atom, atomB: Atom, a: Union[float,Length]):
    """
    Two interpenetrating FCC crystal\n
    Examples:
     - Diamond: both atoms are C
     - Zincblende: basis of Zn & S
    """
    l = Lattice.FCC(a)
    b = AtomBasis(
        (atomA, Vec3D.zero()),
        (atomB, Vec3D.uniform(0.25))
    )
    return (b,l)


def HCP(atomA: Atom, atomB: Atom, a: Union[Length,float], c: Union[Length,float]):
    """Hexagonal Close Packet crystal"""
    l = Lattice.HEX(a,c)
    b = AtomBasis(
        (atomA, Vec3D.zero()),
        (atomB, Vec3D(2/3, 1/3, 0.5)),
    )
    return (b,l)


def WurtziteLike(atomA: Atom, atomB: Atom, a: Union[Length,float], c: Union[Length,float]):
    """
    Two interpenetrating HCP crystal\n
    Examples: wurtzite (ZnS), silver iodide (AgI), zinc oxide (ZnO), cadmium sulfide (CdS), cadmium selenide (CdSe), silicon carbide (α-SiC), gallium nitride (GaN), aluminium nitride (AlN), boron nitride (w-BN)
    """
    l = Lattice.HEX(a,c)
    b = AtomBasis(
        (atomA, Vec3D.zero()),
        (atomB, Vec3D(2/3, 1/3, 0.25)),
        (atomA, Vec3D(2/3, 1/3, 0.5)),
        (atomB, Vec3D(0.0,0.0,0.75))
    )
    return (b,l)


def NiAsLike(Ni: Atom, As: Atom, a: Union[Length,float], c: Union[Length,float]):
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
    l = Lattice.HEX(a,c)
    b = AtomBasis(
        (Ni, Vec3D.zero()),
        (As, Vec3D(2/3, 1/3, 0.25)),
        (Ni, Vec3D(0.0,0.0,0.5)),
        (As, Vec3D(1/3, 2/3, 0.75))
    )
    return (b,l)
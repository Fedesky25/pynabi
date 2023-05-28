from ._common import Vec3D as _V, Stampable as _S
from .units import Pos3D as _P, Length as _L, LUnit as _LU, _def_LU
from typing import Optional as _O, Union as _U

__all__ = ["Atom", "AtomBasis", "Lattice", "cubic", "BCC", "FCC", "interpenetratingFCC", "HCP"]

_atom_symbols = ["H","He","Li","Be","B","C","N","O","F","Ne","Na","Mg","Al","Si","P","S","Cl","Ar","K","Ca","Sc","Ti","V","Cr","Mn","Fe","Co","Ni","Cu","Zn","Ga","Ge","As","Se","Br","Kr","Rb","Sr","Y","Zr","Nb","Mo","Tc","Ru","Rh","Pd","Ag","Cd","In","Sn","Sb","Te","I","Xe","Cs","Ba","La","Ce","Pr","Nd","Pm","Sm","Eu","Gd","Tb","Dy","Ho","Er","Tm","Yb","Lu","Hf","Ta","W","Re","Os","Ir","Pt","Au","Hg","Tl","Pb","Bi","Po","At","Rn","Fr","Ra","Ac","Th","Pa","U ","Np","Pu","Am","Cm","Bk","Cf","Es","F","Md","No","Lr","Rf","Db","Sg","Bh","Hs","Mt","Ds","Rg","Cn","Nh","Fl","Mc","Lv","Ts","Og"]

class Atom:
    def __init__(self, Z: int, file: str):
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
        
        returns Atom(Z of element, element name + 'psp8')
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
        return f"""natoms{suffix} {len(indexes)}
typeat{suffix} {' '.join(str(i+1) for i in indexes)}
{x_type}{suffix} {'   '.join(str(a[1]) for a in self.atoms)}"""
        

class Lattice(_S):
    def __init__(self, scaling: _U[_V,_P,None], prop: _O['tuple[str,str]']):
        """Do not use directly: prefer fromAngle, fromPrimitives, setScaling"""
        self.scaling = scaling
        self.prop = prop

    def stamp(self, index: int):
        suffix = index if index > 0 else ''
        res: list[str] = []
        if self.scaling is not None:
            res.append(f"acell{suffix} {self.scaling}")
        if self.prop is not None:
            res.append(f"{self.prop[0]}{suffix} {self.prop[1]}")
        return '\n'.join(res)

    def isCompleteWith(self, other: _O['Lattice']):
        s = self.scaling is not None
        p = self.prop is not None
        if other is None:
            return s and p
        else:
            return (s or other.scaling is not None) and (p or other.prop is not None)
    
    @staticmethod
    def fromAngles(angles: _V, scaling: _U[_V,_P,None] = None):
        """Constructs lattice from angles [&alpha;, &beta;, &gamma;]

        Angles are defined between vectors as follows:
        | angle | from | to |
        | :---: | :----------: | :-----------: |
        | &alpha;(1) | b(2) | c(3) |
        | &beta;(2) | a(1) | c(3) |
        | &gamma;(3) | a(1) | b(2) |
        """
        return Lattice(scaling, ("angdeg", str(angles)))
    
    @staticmethod
    def fromPrimitives(a: _V, b: _V, c: _V, scaling: _U[_V,_P,None] = None):
        """Construct lattice from primitives [a, b, c]"""
        return Lattice(scaling, ("rprim", f"{a}\n  {b}\n  {c}"))

    @staticmethod
    def setScaling(scaling: _V):
        return Lattice(scaling, None)


def _2uni(v: _U[float,_L]):
    t = type(v)
    if t is _L:
        return _P.uniform(v) # type: ignore
    elif t is float:
        return _V.uniform(v) # type: ignore
    else:
        raise TypeError(f"Lattice constant is of wrong value (got {t} instead of float or Length)")


def cubic(a: _U[float,_L], atom: Atom):
    return (
        AtomBasis((atom, _V.zero())),
        Lattice.fromAngles(_V.uniform(90),_2uni(a))
    )

def BCC(a: _U[float,_L], atom: Atom, center: _O[Atom] = None):
    r = cubic(a, atom)
    r[0].add(atom if center is None else center, _V.uniform(0.5))
    return r

def FCC(a: _U[float,_L], atom: Atom, face: _O[Atom] = None):
    r = cubic(a, atom)
    f = atom if face is None else face
    r[0].add(f, _V(0.5,0.5,0.0))
    r[0].add(f, _V(0.0,0.5,0.5))
    r[0].add(f, _V(0.5,0.0,0.5))
    return r


def interpenetratingFCC(a: _U[float,_L], atomA: Atom, atomB: Atom):
    """
    Two interpenetrating FCC crystal (Rhombohedral primitive cell)\n
    Examples: ZincBlende
    """
    l = Lattice.fromPrimitives(
        _V(0.5, 0.5, 0.0),
        _V(0.0, 0.5, 0.5),
        _V(0.5, 0.0, 0.5),
        _2uni(a)
    )
    b = AtomBasis(
        (atomA, _V.zero()),
        (atomB, _V.uniform(0.25))
    )
    return (b,l)


def HCP(atomA: Atom, atomB: Atom, a: float, c: float, unit: _LU = _def_LU):
    return (
        AtomBasis(
            (atomA, _V.zero()),
            (atomB, _V(2/3, 1/3, 0.5)),
        ),
        Lattice.fromAngles(
            _V(90,120,90),
            _P(a,a,c,unit)
        )
    )
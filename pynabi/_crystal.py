from ._common import Vec3D, Stampable
from typing import Union, List, Tuple

__all__ = ["Atom", "AtomBasis", "Lattice"]

atoms = ["H","He","Li","Be","B","C","N","O","F","Ne","Na","Mg","Al","Si","P","S","Cl","Ar","K","Ca","Sc","Ti","V","Cr","Mn","Fe","Co","Ni","Cu","Zn","Ga","Ge","As","Se","Br","Kr","Rb","Sr","Y","Zr","Nb","Mo","Tc","Ru","Rh","Pd","Ag","Cd","In","Sn","Sb","Te","I","Xe","Cs","Ba","La","Ce","Pr","Nd","Pm","Sm","Eu","Gd","Tb","Dy","Ho","Er","Tm","Yb","Lu","Hf","Ta","W","Re","Os","Ir","Pt","Au","Hg","Tl","Pb","Bi","Po","At","Rn","Fr","Ra","Ac","Th","Pa","U ","Np","Pu","Am","Cm","Bk","Cf","Es","F","Md","No","Lr","Rf","Db","Sg","Bh","Hs","Mt","Ds","Rg","Cn","Nh","Fl","Mc","Lv","Ts","Og"]

class Atom:
    def __init__(self, Z: int, file: str):
        self.num = Z
        self.file = file
        assert self.file.find(",") == -1, "File name cannot contain ','"
    
    def __str__(self) -> str:
        return f"{atoms[self.num-1]} atom (located at {self.file})"
    
    # def __eq__(self, __value: object) -> bool:
    #     return (self is __value) or (type(__value) is Atom and __value.num == self.num and __value.file == self.file)

    @staticmethod
    def of(name: str):
        """Constructs an atom using default filenaming
        
        returns Atom(Z of element, element name + 'psp8')
        """
        try:
            Z = atoms.index(name) + 1;
            return Atom(Z, name + '.psp8')
        except:
            raise ValueError(name + " is not a valid atom type")
    
    @staticmethod
    def poolstr(atoms: 'list[Atom]'):
        return f"""# Atoms definition
ntypat {len(atoms)}
znucl {' '.join(str(a.num) for a in atoms)}
pseudos \"{', '.join(a.file for a in atoms)}\""""


class AtomBasis(Stampable):
    def __init__(self, *atoms: 'tuple[Atom,Vec3D]', cartesian: bool = False) -> None:
        assert len(atoms) > 0, "There must be at least one atom in basis"
        self.atoms = atoms
        self.cartesian = cartesian
    
    def add(self, atom: Atom, where: Vec3D):
        if type(self.atoms) is tuple:
            self.atoms = list(self.atoms)
        self.atoms.append((atom, where)) # type: ignore
    
    def getAtoms(self):
        return (a[0] for a in self.atoms)
    
    def stamp(self, index: int, pool: List[Atom]):
        indexes = [pool.index(a[0]) for a in self.atoms]
        suffix = str(index or '');
        x_type = "xcart" if self.cartesian else "xred";
        return f"""natoms{suffix} {len(indexes)}
typeat{suffix} {' '.join(str(i+1) for i in indexes)}
{x_type}{suffix} {'   '.join(str(a[1] for a in self.atoms))}"""
        

class Lattice(Stampable):
    def __init__(self, scaling: Union[None,Vec3D], prop: Union[None,Tuple[str,str]]):
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

    def isCompleteWith(self, other: Union['Lattice',None]):
        s = self.scaling is not None
        p = self.prop is not None
        if other is None:
            return s and p
        else:
            return (s or other.scaling is not None) and (p or other.prop is not None)
    
    @staticmethod
    def fromAngles(angles: Vec3D, scaling: Union[None,Vec3D] = None):
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
    def fromPrimitives(a: Vec3D, b: Vec3D, c: Vec3D, scaling: Union[None,Vec3D] = None):
        """Construct lattice from primitives [a, b, c]"""
        return Lattice(scaling, ("rprim", f"{a}\n  {b}\n  {c}"))

    @staticmethod
    def setScaling(scaling: Vec3D):
        return Lattice(scaling, None)
from ._common import Vec3D, sectionTitle, Stampable
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
    def __init__(self, *atoms: 'tuple[Atom,Vec3D]', cartesian: bool = False, pseudosAt: str = '') -> None:
        self.atoms = atoms
        self.cartesian = cartesian
        self.dir = pseudosAt
    
    def add(self, atom: Atom, where: Vec3D):
        if type(self.atoms) is tuple:
            self.atoms = list(self.atoms)
        self.atoms.append((atom, where)) # type: ignore
    
    def getAtoms(self):
        return (a[0] for a in self.atoms)
    
    def stamp(self, index: int, pool: List[Atom]):
        indexes = [pool.index(a[0]) for a in self.atoms]
        suffix = str(index or '');
        res: list[str] = [sectionTitle(index, 'Atom basis')]
        l = len(indexes)
        if l > 0:
            res.append(f"natoms{suffix} {l}")
            res.append(f"typeat{suffix} {' '.join(str(i+1) for i in indexes)}")
            x_type = "xcart" if self.cartesian else "xred";
            x_space = '\n  ' + (' '*(len(x_type)+len(suffix)))
            res.append(f"{x_type}{suffix} {x_space.join([str(a[1]) for a in self.atoms])}")
        if len(self.dir) > 0:
            res.append(f"pp_dirpath{suffix} \"{self.dir}\"")
        return '\n'.join(res)
        

class Lattice(Stampable):
    def __init__(self, scaling: Union[None,Vec3D], prop: Union[None,Tuple[str,str]]):
        self.scaling = scaling
        self.prop = prop

    def stamp(self, index: int):
        suffix = index if index > 0 else ''
        res: list[str] = [sectionTitle(index, "Lattice")]
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
        """Constructs lattice from angles

        Angles are defined between vectors as follows:
        | index | from | to |
        | :----: | :----------: | :-----------: |
        | 1 | 2 | 3 |
        | 2 | 1 | 2 |
        | 3 | 1 | 3 |
        """
        return Lattice(scaling, ("angdeg", str(angles)))
    
    @staticmethod
    def fromPrimitives(x: Vec3D, y: Vec3D, z: Vec3D, scaling: Union[None,Vec3D] = None):
        """Construct lattice from primitives"""
        return Lattice(scaling, ("rprim", f"{x}\n  {y}\n  {z}"))

    @staticmethod
    def setScaling(scaling: Vec3D):
        return Lattice(scaling, None)
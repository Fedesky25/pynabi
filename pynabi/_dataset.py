from typing import Union, Tuple, List, Literal, Iterable
from ._common import Stampable
from ._crystal import AtomBasis, Atom
from enum import Enum

__all__ = ["DataSet", "PreviousRun", "I", "O", "AbIO", "createAbi"]

def splat(i: Iterable[Union[Stampable,Iterable[Stampable]]]) -> Iterable[Stampable]:
    for v in i:
        if isinstance(v,Stampable):
            yield v
        else:
            try:
                yield from splat(v)
            except:
                raise TypeError("Arguments provided to dataset must be DataSet stampables of iterables of them")

class DataSet:
    def __init__(self, *stampables: Union[Stampable,Iterable[Stampable]]) -> None:
        self.index = 0
        self.atoms: Union[AtomBasis,None] = None
        self.stamps: list[Stampable] = []
        types = set()
        for s in splat(stampables):
            t = type(s)
            if t in types:
                raise ValueError(f"Multiple {s.__class__.__name__} given")
            types.add(t)
            if t is AtomBasis:
                self.atoms = s # type: ignore
            else:
                self.stamps.append(s) # type: ignore
    
    def stamp(self, atompool: List[Atom]):
        res: list[str] = []
        if self.atoms is not None:
            res.append(self.atoms.stamp(self.index, atompool))
        for s in self.stamps:
            res.append(s.stamp(self.index))
        return '\n'.join(res)


class PreviousRun(object):
    _instance = None
    def __new__(cls):
        if cls._instance is None:
            cls._instance = object.__new__(cls)
        return cls._instance


class I(Enum):
    FirstOrderDensity = "1den", 0
    FirstOrderWavefunction = "1wf", 0
    BetheSalpeterCouplingBlock = "bscoup", 0
    BetheSalpeterEigenstates = "bseig", 0
    BetheSalpeterResonantBlock = "bsreso", 0
    DDB = "ddb", 1
    DDKWavefunctions = "ddk", 0
    dE = "delfd", 0
    ElectronDensity = "den", 1
    dkdE = "dkde", 0
    dkdk = "dkdk", 0
    PotentialDerivative = "dvdb", 1
    EffectiveMasses = "efmas", 0
    HaydockRestartFile = "haydock", 0
    OccupationNumbers = "occ", 0
    KSPotential = "pot", 2
    QuasiParticleStructure = "qps", 0
    Screening = "scr", 1
    SIGEPH = "sigeph", 2
    Susceptibility = "suscep", 0
    WavefunctionsK = "wfk", 1
    FineWavefunctionsK = "wfkfine", 2
    WavefunctionsQ = "wfq", 1
    #chkprdm

    def ok(self, value):
        t = type(value)
        if self.value[1] == 0 and t is str:
            raise ValueError(f"{self.name} cannot be read from a file")
        if self.value[1] == 2 and t is not str:
            raise ValueError(f"{self.name} can only be read from a file")

    def stamp(self, value: Union[PreviousRun, DataSet, str], index: Union[int,Literal['']]):
        t = type(value)
        if t is PreviousRun:
            return f"ird{self.value[0]}{index or ''} 1"
        elif t is str:
            return f"get{self.value[0]}_filepath{index or ''} \"{value}\""
        else:
            if value.index >= (index or 0): # type: ignore
                raise IndexError(f"Cannot read {self.name} for {index}-th dataset from the {value.index}-th dataset") # type: ignore
            return f"get{self.value[0]}{index or ''} {value.index}" # type: ignore


class O(Enum):
    PotentialAndDensity1D = "1dm"
    CheckPoint = "chkprdm"
    ElectronDensity = "den"
    DensityOfStates = "dos"
    MResolvedPartialDOS = "dosm"
    EigenEnergies = "eig"
    ElectronLocalizedFunction = "elf"
    FermiSurface = "surf"
    ElectronDensityGradient = "gden"
    GeometryAnalysis = "geo"
    MatrixGKK = "gkk"
    GSR = "gsr"
    KleynmanBylanderFormFactors = "kbff"
    KineticEnergyDensity = "kden"
    KPointsSets = "kpt"
    ElectronDensityLaplacian = "lden"
    Potential = "pot"
    PSPS = "psps"
    SpinCurrentDensity = "spcur"
    STMDensity = "stm"
    Susceptibility = "suscep"
    CoulombPotential = "vclmb"
    HartreePotential = "vha"
    HartreeANdExchangeCorrelationPotential = "vhxc"
    Volume = "vol"
    VolumeForImages = "volimg"
    LocalPseudoPotential = "vpsp"
    ExchangeCorrelationPotential = "vxc"
    WanT = "want"
    Wavefunction = "wf"
    FullMeshWavefunction = "wf_full"
    XML = "xml"


class AbIO(Stampable):
    def __init__(self) -> None:
        self._in: dict[I, Union[PreviousRun, DataSet, str]] = {}
        self._out: dict[O, int] = {}
        self._prefixes = { "out": '', "in": '', "tmp": '' }

    def prefix(self, input: str = '', output: str = '', temporary: str = ''):
        self._prefixes["in"] = input
        self._prefixes["out"] = output
        self._prefixes["tmp"] = temporary
        return self
    
    def print(self, *what: Union[O, Tuple[O,int]]):
        for o in what:
            if type(o) is O:
                self._out[o] = 1
            elif type(o) is tuple:
                assert type(o[0]) is O and type(o[1]) is int, f"Print tuple must be (name, int): got ({o[0]}, {o[1]})"
                self._out[o[0]] = o[1]
            else:
                raise ValueError(f"{o} is not a valid print option")
        return self

    def unprint(self, *what: O):
        for o in what:
            assert type(o) is O, f"{o} is not a valid print option"
            self._out[o] = 0
        return self
    
    def get(self, what: I, _from: Union[PreviousRun, DataSet, str]):
        what.ok(_from)
        self._in[what] = _from;
        return self
    
    def stamp(self, index: int):
        suffix = index or ''
        res: list[str] = []
        tmp: list[str] = []
        for (key, path) in self._prefixes.items():
            if len(path) > 0:
                tmp.append(f"{key}data_prefix{suffix} \"{path}\"")
        if len(tmp) > 0:
            res.extend(tmp)
        tmp = []
        if index < 2:
            for (i, val) in self._in.items():
                if type(val) is not str:
                    raise ValueError(f"First dataset can only read variables from files (at {i.name})")
                tmp.append(i.stamp(val, suffix))
        else:
            for (i,val) in self._in.items():
                tmp.append(i.stamp(val, suffix))
        if len(tmp) > 0:
            res.extend(tmp)
        tmp = []
        for (o,n) in self._out.items():
            tmp.append(f"prt{o.value}{suffix} {n}")
        if len(tmp) > 0:
            res.extend(tmp)
        return '\n'.join(res)


def createAbi(setup: Union[DataSet,None], *datasets: DataSet):
    n = len(datasets)
    if setup is None:
        if n == 0:
            return ''
        elif n == 1:
            return createAbi(datasets[0])
        else:
            setup = DataSet()
    elif len(datasets) == 1:
        raise ValueError("Cannot use a single dataset")
    
    res: list[str] = [f"ndtset {len(datasets)}\n"]

    atomSet = set() if setup.atoms is None else set(setup.atoms.getAtoms())
    initialAtomCount = len(atomSet)
    # isBaseLatticeComplete = setup.lattice is not None and setup.lattice.isCompleteWith(None) 
    for (i,s) in enumerate(datasets):
        s.index = i+1
        # assert isBaseLatticeComplete or (s.lattice is not None and s.lattice.isCompleteWith(setup.lattice)), f"{i}-th lattice is not completely defined"
        if s.atoms is not None:
            atomSet = atomSet.union(s.atoms.getAtoms())
        elif initialAtomCount == 0:
            raise ValueError(f"All datasets (in particular the {i+1}-th one) must define the atom basis since no common one was defined")
    
    atomPool = list(atomSet)
    res.append(Atom.poolstr(atomPool))
    res.append("\n# Common DataSet")
    res.append(setup.stamp(atomPool))
    for s in datasets:
        res.append(f"\n# DataSet {s.index}")
        res.append(s.stamp(atomPool))
    return '\n'.join(res)
    

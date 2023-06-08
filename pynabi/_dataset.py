from typing import Union, List, Iterable, Literal, Callable, Optional, Type
from ._common import Stampable, Singleton, Delayed, StampCollection
from .crystal import AtomBasis, Atom
from inspect import stack


__all__ = ["DataSet", "PreviousRun", "AbIn", "AbOut", "createAbi"]


_RS = Union[Stampable,Iterable['_RS']]


def splat(i: Iterable[_RS]) -> Iterable[Stampable]:
    for v in i:
        if isinstance(v,Stampable):
            yield v
        else:
            try:
                yield from splat(v)
            except:
                raise TypeError("Arguments provided to dataset must be DataSet stampables of iterables of them")


class DataSet:
    def __init__(self, *stampables: _RS) -> None:
        self.index = 0
        self.atoms: Union[AtomBasis,None] = None
        self.stamps: list[Stampable] = []
        self.map: dict[Type[Stampable], Stampable] = {}
        delayed_props = set()
        for s in splat(stampables):
            if not isinstance(s, Stampable):
                raise TypeError(f"{s.__class__.__name__} is not a valid type for a dataset")
            if isinstance(s, Delayed):
                p = s.getProp()
                if p in delayed_props:
                    raise ValueError(f"Multiple {s.getName()} definition are present")
                delayed_props.add(p)
            else:
                t = type(s)
                if t in self.map:
                    raise ValueError(f"Multiple {s.__class__.__name__} given")
                self.map[t] = s
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


class PreviousRun(Singleton):
    pass


def _AbInMethod(prop: str, sel: Literal[0,1,2] = 0):
    def method(self: 'AbIn', value: Union['DataSet', str, PreviousRun]): 
        if type(value) is str:
            assert sel != 0, f"{method._name} cannot be read from a file"
        else:
            assert sel < 2, f"{method._name} can only be read from a file"
        self._d[method] = value
        return self
    method._name = ""
    method._prop = prop
    return method


class AbIn(Stampable):
    def __init__(self, prefix: Optional[str] = None):
        self._p = prefix
        self._d: dict[Callable, Union['DataSet', str, PreviousRun]] = dict()
        self._ppd: Optional[str] = None
    
    def stamp(self, index: int):
        res: list[str] = []
        if self._ppd is not None:
            res.append(f"pp_dirpath{index or ''} \"{self._ppd}\"")
        if self._p is not None:
            res.append(f"indata_prefix{index or ''} \"{self._p}\"")
        body = '\n'.join(AbIn._print_helper(k,v,index) for k, v in self._d.items())
        if len(body) > 0:
            res.append(body)
        return '\n'.join(res);
    
    @staticmethod
    def _print_helper(m: Callable, v: Union['DataSet', str, PreviousRun], i: int):
        if type(v) is DataSet:
            assert v.index < i, f"Cannot read {m._name} for {i}-th dataset from the {v.index}-th dataset"
            return f"get{m._prop}{i or ''} {v.index}"
        elif type(v) is PreviousRun:
            return f"ird{m._prop}{i or ''} 1"
        else:
            return f"get{m._prop}_filepath{i or ''} \"{v}\""
    
    def PseudoPotentials(self, path: str):
        self._ppd = path;
        return self
            
    FirstOrderDensity = _AbInMethod("1den")
    FirstOrderWavefunction = _AbInMethod("1wf")
    BetheSalpeterCouplingBlock = _AbInMethod("bscoup")
    BetheSalpeterEigenstates = _AbInMethod("bseig")
    BetheSalpeterResonantBlock = _AbInMethod("bsreso")
    DDB = _AbInMethod("ddb",1)
    DDKWavefunctions = _AbInMethod("ddk")
    dE = _AbInMethod("delfd")
    ElectronDensity = _AbInMethod("den",1)
    dkdE = _AbInMethod("dkde")
    dkdk = _AbInMethod("dkdk")
    PotentialDerivative = _AbInMethod("dvdb",1)
    EffectiveMasses = _AbInMethod("efmas")
    HaydockRestartFile = _AbInMethod("haydock")
    OccupationNumbers = _AbInMethod("occ")
    KSPotential = _AbInMethod("pot",2)
    QuasiParticleStructure = _AbInMethod("qps")
    Screening = _AbInMethod("scr",1)
    SIGEPH = _AbInMethod("sigeph",2)
    Susceptibility = _AbInMethod("suscep")
    WavefunctionsK = _AbInMethod("wfk",1)
    FineWavefunctionsK = _AbInMethod("wfkfine",2)
    WavefunctionsQ = _AbInMethod("wfq",1)


for k,v in AbIn.__dict__.items():
    if callable(v) and hasattr(v,"_name"):
        v._name = k


def _os(name: str):
    def method(self: 'AbOut', flag: bool = True):
        self._d[name] = int(flag)
        return self
    return method


def _om(name: str, max: int, min: int = 0):
    def method(self: 'AbOut', value: int):
        assert type(value) is int and min <= value <= max, f"Value of {stack()[0][3]} must be integer between {min} and {max}"
        self._d[name] = value
        return self
    return method


class AbOut(Stampable):
    def __init__(self, prefix: Union[str,None] = None):
        self._p = prefix
        self._d: dict[str, int] = {}

    def stamp(self, index: int):
        s = index or ''
        body = '\n'.join(f"prt{k}{s} {v}" for k,v in self._d.items())
        if self._p is None:
            return body
        elif len(body) == 0:
            return f"outdata_prefix \"{self._p}\""
        else:
            return f"outdata_prefix \"{self._p}\"\n{body}"

    PotentialAndDensity1D = _os("1dm")
    CheckPoint = _os("chkprdm")
    ElectronDensity = _os("den")
    DensityOfStates = _om("dos", 5)
    MResolvedPartialDOS = _os("dosm")
    EigenEnergies = _os("eig")
    ElectronLocalizedFunction = _os("elf")
    FermiSurface = _os("surf")
    ElectronDensityGradient = _os("gden")
    GeometryAnalysis = _os("geo")
    MatrixGKK = _os("gkk")
    GSR = _os("gsr")
    KleynmanBylanderFormFactors = _os("kbff")
    KineticEnergyDensity = _os("kden")
    KPointsSets = _os("kpt")
    ElectronDensityLaplacian = _os("lden")
    Potential = _os("pot")
    PSPS = _os("psps")
    SpinCurrentDensity = _os("spcur")
    STMDensity = _os("stm")
    Susceptibility = _os("suscep")
    CoulombPotential = _os("vclmb")
    HartreePotential = _os("vha")
    HartreeANdExchangeCorrelationPotential = _os("vhxc")
    Volume = _om("vol", 11, -10)
    ImagesVolume = _os("volimg")
    LocalPseudoPotential = _os("vpsp")
    ExchangeCorrelationPotential = _os("vxc")
    UseWanTInterface = _om("want", 3)
    Wavefunction = _om("wf", 2, -1)
    FullMeshWavefunction = _os("wf_full")
    XML = _os("xml")


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
    base_coll = {} if setup is None else setup.map
    for (i,d) in enumerate(datasets):
        d.index = i+1
        coll = StampCollection(d.map, base_coll)
        for s in d.stamps:
            s.compatible(coll)
        if d.atoms is not None:
            atomSet = atomSet.union(d.atoms.getAtoms())
        elif initialAtomCount == 0:
            raise ValueError(f"All datasets (in particular the {i+1}-th one) must define the atom basis since no common one was defined")
    
    atomPool = list(atomSet)
    res.append(Atom.poolstr(atomPool))
    res.append("\n# Common DataSet")
    res.append(setup.stamp(atomPool))
    for d in datasets:
        res.append(f"\n# DataSet {d.index}")
        res.append(d.stamp(atomPool))
    return '\n'.join(res)
    

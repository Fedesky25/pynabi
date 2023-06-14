from ._common import Stampable as _S, StampCollection as _SC
from .units import Energy as _En
from typing import Any as _A, Literal as _L, Tuple as _T


# TODO: missing ionmov 28


def _pos_num(v):
    return (type(v) is int or type(v) is float) and v >= 0


# class FixedAtoms(_S):
#     pass


def _check_co(coll: _SC, name: str):
    coov = coll.get(CellOptimizationOnVector)
    if coov is not None:
        raise ValueError(f"{name} is not compatible with cell optimization on vector")
    co = coll.get(CellOptimization)
    if co is not None:
        raise ValueError(f"{name} is not compatible with cell optimization")


def _mdtemp(v):
    assert type(v) is tuple and len(v) == 2 and _pos_num(v[0]) and _pos_num(v[1]), "Temperatures must be a tuple of to positive numbers"
    return f"[{v[0]},{v[1]}]"


# class _MD_SO_def:
#     def __init__(self, fn) -> None:
#         self.fn = fn
    
#     def __call__(self, *args: _A, **kwds: _A) -> _A:
#         self.fn(*args, **kwds)
    
#     def __set_name__(self, owner: type, name: str):
#         setattr(owner, "_def", self.fn)
#         setattr(owner, name, self.fn) # replace ourselves with the original method


class _MD_SO(_S):
    def __init__(self, timeStep: float = 100, maxSteps: int = 1000) -> None:
        super().__init__()
        assert _pos_num(timeStep), "time step (for molecular dynamics or structural optimization) must be a positive number"
        assert type(maxSteps) is int and maxSteps > 0, "number of maximum steps (for molecular dynamics or structural optimization) must be a positive integer"
        self._dt = timeStep
        self._ms = maxSteps
        self._n = 0
        self._d: dict[str,_A] = {}
    
    def stamp(self, index: int):
        s = index or ''
        r = f"ionmov{s} {self._n}"
        extra = '\n'.join(f"{k}{s} {v}" for (k,v) in self._d.items())
        if len(extra) > 0:
            r += "\n" + extra
        return r
    
    def compatible(self, coll: _SC):
        assert self._n != 0, f"Incomplete definition of {type(self).__name__}"


class MolecularDynamics(_MD_SO):
    """Move atoms using molecular dynamics"""

    def compatible(self, coll: _SC):
        if self._n != 13: # isoenthaplic is the only one that accept cell optimization
            _check_co(coll, "Molecular dynamics")
    
    def basic(self):
        """Move atoms using undumped molecular dynamics\n
        For viscous damping use `StructuralOptimization(...).damping(...)`"""
        self._n = 1
        self._d = {"vis": 0}
        return self

    def Verlet(self):
        """Molecular dynamics using the Verlet algorithm, see [[Allen1987a]](https://docs.abinit.org/theory/bibliography#allen1987a) p 81"""
        self._n = 6
        self._d = {}
        return self

    def quenchedVerlet(self):
        """Quenched Molecular dynamics using the Verlet algorithm, and stopping each atom for which the scalar product of velocity and force is negative.\n
        The goal is not to produce a realistic dynamics, but to go as fast as possible to the minimum. For this purpose, it is advised to set all the masses to the same value"""
        self._n = 7
        self._d = {}
        return self
    
    def NoseHoover(self, inertia: float = 100, temperatures: _T[float, float] = (200,300)):
        """
        Molecular dynamics with Nose-Hoover thermostat, using the Verlet algorithm
        
        `inertia` is the inertia factor WT of the thermostat in atomic units (see [the docs](https://docs.abinit.org/variables/rlx/#noseinert))

        `temperatures` is a tuple of floats who represent the initial and final temperatures of the thermostat
        """
        assert _pos_num(inertia), "thermostat intertia must a positive number"
        self._n = 8
        self._d = { "mdtemp": _mdtemp(temperatures), "noseinert": inertia }
        return self

    def Langevin(self, friction: float, temperatures: _T[float, float] = (200,300)):
        """
        Langevin molecular dynamics

        `friction` is friction coefficient for Langevin dynamics with units Hartree*Electronic mass*(atomic unit of Time)/Bohr^2

        `temperatures` is a tuple of floats who represent the initial and final temperatures of the dynamics
        """
        assert _pos_num(friction), "thermostat intertia must a positive number"
        self._n = 9
        self._d = { "mdtemp": _mdtemp(temperatures), "friction": friction }
        return self
    
    def isokinetic(self, temperature: float = 200):
        """
        Isokinetic ensemble molecular dynamics\n
        The velocity is initialized from the given temperature (in Kelvin)
        """
        assert _pos_num(temperature), "The temperature must be a positive number"
        self._n = 12
        self._d = { "mdtemp": f"[{temperature},{temperature+1}]" }
        return self

    def SRKNa14(self):
        """Simple molecular dynamics with a symplectic algorithm proposed in [[Blanes2002]](https://docs.abinit.org/theory/bibliography#blanes2002) (called SRKNa14) of the kind first published in [[Yoshida1990]](https://docs.abinit.org/theory/bibliography#bitzek2006). This algorithm requires at least 14 evaluation of the forces (actually 15 are done within Abinit) per time step. At this cost it usually gives much better energy conservation than the verlet algorithm for a 30 times bigger value of time step.\n
        NOTE: the potential energy of the initial atomic configuration is never evaluated using this algorithm."""
        self._n = 14
        self._d = {}
        return self

    def learnOnTheFly(self, iterations: int = 10):
        """Use of Learn on The Fly method (LOTF) for Molecular Dynamics. In the framework of isokinetic MD, the atomic forces and positions are computed by using LOTF interpolation. A SCF computation is performed only _A `iterations` steps. The results of the SCF are used to compute the parameters of a short range classical potential (for the moment only the glue potential for gold is implemented). Then these parameters are continuously tuned to compute atomic trajectories.\n
        The LOTF cycle is divided in the following steps: 
         1. Initialization (SFC at t=0) and computation of potential parameters. 
         2. Extrapolation of the atomic forces and positions for `iterations` time step. To perform this extrapolation, the potential computed in 1 is used (Verlet algorithm). 
         3. SFC at t=`iterations`. Computation of the potential parameters. 
         4. LOTF interpolation, linear interpolation of the potential parameters and computation of the atomic forces and positions between t=0 and t=`iterations`.
        
        NOTE: LOTF has to be enabled at configure time. If LOTF is not enabled, abinit will set automatically isokinetic MD.\n
         """
        assert type(iterations) is int and iterations > 0, "Number of LOTF iterations must be a positive integer"
        self._n = 23
        self._d = {"lotf_nitex": iterations }
        return self

    def constantEnergy(self):
        """Simple constant energy molecular dynamics using the velocity Verlet symplectic algorithm (second order), see [[Hairer2003]](https://docs.abinit.org/theory/bibliography#hairer2003)."""
        self._n = 24
        self._d = {}
        return self


class StructuralOptimization(_MD_SO):
    def compatible(self, coll: _SC):
        if self._n not in (2, 3, 22):
            _check_co(coll, "Structural optimization")

    def damping(self, viscosity: float = 100):
        """Move atoms using molecular dynamics with viscous damping (friction linearly proportional to velocity). \n
        The implemented algorithm is the generalisation of the Numerov technique (6th order), but is NOT invariant upon time-reversal, so that the energy is not conserved."""
        assert _pos_num(viscosity), "Viscosity must be a positive number"
        self._n = 1
        self._d = {"vis": viscosity}
        return self

    def BFGS(self, delocalizedCoordinates: bool = False, withEnergy: bool = False):
        """Conduct structural optimization using the Broyden-Fletcher-Goldfarb-Shanno minimization (BFGS)\n
         - `withEnergy` determines whether to take into account the total energy or not (can be very unstable)
         - `delocalizedCoordinates` detemines whether to use delocalized coordinates or not (if true, no cell optimization can be done) 
        """
        assert type(delocalizedCoordinates) is bool, "delocalizedCoordinates must be a bool"
        assert type(withEnergy) is bool, "withEnergy must be a bool"
        self._n = 2 + 8*int(delocalizedCoordinates) + int(withEnergy)
        self._d = {}
        return self
    
    def conjugateGradient(self):
        """Conjugate gradient algorithm for simultaneous optimization of potential and ionic degrees of freedom.\n
        WARNING: this is under development, and does not work very well in m_A cases"""
        self._n = 4
        self._d = {}
        return self

    def simpleRelaxation(self, preconditioning: _L[-1, 0, 1, 2] = 0):
        """Simple relaxation of ionic positions according to (converged) forces\n
        `preconditioning` is an integer between -1 and 2 and describes the way a change of force is derived from a change of atomic position. 
        In particular: hessian is 2^(-preconditioning) times the identity matrix
        """
        assert type(preconditioning) is int and -1 <= preconditioning <= 2, "preconditioning must a integer between -1 and 2"
        self._n = 5
        self._d = { "iprcfc": preconditioning }
        return self

    # def directInversion(self):
    #     self._n = 20
    #     self._d = {}
    #     return self

    def LBFGS(self):
        """Conduct structural optimization using the Limited-memory Broyden-Fletcher-Goldfarb-Shanno minimization (L-BFGS) [[Nocedal1980]](https://docs.abinit.org/theory/bibliography#nocedal1980). 
        The routines are based on the original implementation by J. Nocedal available on netlib.org. 
        This algorithm can be much better than the native implementation of BFGS in ABINIT"""
        self._n = 22
        self._d = {}
        return self


class FIRE(_S):
    """Fast inertial relaxation engine (FIRE) algorithm proposed by Erik Bitzek, Pekka Koskinen, Franz Gähler, Michael Moseler, and Peter Gumbsch in [[Bitzek2006]](https://docs.abinit.org/theory/bibliography#bitzek2006). According to the authors, the efficiency of this method is nearly the same as L-BFGS. It is based on conventional molecular dynamics with additional velocity modifications and adaptive time steps. The purpose of this algorithm is relaxation, not molecular dynamics. \n
    The initial time step is set with `initialTimeStep`: it governs the ion position changes, but the cell parameter changes as well. The suggested first guess is 0.03.\n
    The positions are in reduced coordinates instead of in cartesian coordinates."""
    def __init__(self, initialTimeStep: float = 0.03) -> None:
        super().__init__()
        assert _pos_num(initialTimeStep), "The initial time step must be a positive number"
        self.t = initialTimeStep
    
    def stamp(self, index: int):
        s = index or ''
        return f"ionmov{s} 15\ndtion{s} {self.t}"


class MonteCarloSampling(_S):
    """Hybrid Monte Carlo sampling of the ionic positions at fixed temperature and unit cell geometry (NVT ensemble). 
    The underlying molecular dynamics corresponds to the constant energy one. \n
    Within the HMC algorithm [[Duane1987]](https://docs.abinit.org/theory/bibliography#duane1987), the trial states are generated via short trajectories (ten steps in current implementation). The initial momenta for each trial are randomly sampled from Boltzmann distribution, and the final trajectory state is either accepted or rejected based on the Metropolis criterion. Such strategy allows to simultaneously update all reduced coordinates, achieve higher acceptance ratio than classical Metropolis Monte Carlo and better sampling efficiency for shallow energy landscapes [[Prokhorenko2018]](https://docs.abinit.org/theory/bibliography#prokhorenko2018)"""

    def __init__(self, temperatures: _T[float, float] = (200,300), timeStep: float = 100) -> None:
        assert _pos_num(timeStep)
        self.t = _mdtemp(temperatures)
        self.s = timeStep
    
    def stamp(self, index: int):
        s = index or ''
        return f"ionmov{s} 24\ndtion{s} {self.s}\nmdtemp{s} {self.t}"

    def compatible(self, coll: _SC):
        _check_co(coll, "Monte Carlo Sampling")


def _CO_method(n: int):
    def fn(self: 'CellOptimization'):
        assert self.n == 0, "Cell optimization defined multiple times"
        self.n = n
        return self
    return fn


class CellOptimization(_S):
    def __init__(self, energyCutoffSmearing: _En|float) -> None:
        super().__init__()
        self.s = _En.sanitize(energyCutoffSmearing)
        self.n = 0
    
    def stamp(self, index: int):
        s = index or ''
        return f"optcell{s} {self.n}\necutsm{s} {self.s}"
    
    def compatible(self, coll: _SC):
        assert self.n != 0, "Cell optimization not fully defined"
    
    OnVolumeOnly = _CO_method(1)
    Full = _CO_method(2)
    ConstantVolume = _CO_method(3)


class CellOptimizationOnVector(_S):
    def __init__(self, energyCutoffSmearing: _En|float, axis: _L[1,2,3,'x','y','z','X','Y','Z'], keep: bool) -> None:
        super().__init__()
        self.s = _En.sanitize(energyCutoffSmearing)
        if type(axis) is int:
            assert 1 <= axis <= 3, "Axis index must be between 1 and 3"
            self._a = axis
        elif type(axis) is str:
            n = ('x','y','z')
            s = axis.lower()
            assert s in n, "Axis name can only be one of the following charachters: xyxXYZ"
            self._a = n.index(n)
        else:
            raise TypeError("Invalid type for axis argument")
        assert type(keep) is bool, "Argument keeep must be a bool"
        self._k = keep
    
    def stamp(self, index: int):
        s = index or ''
        return f"optcell{s} {4 + self._a + 3*int(self._k)}\necutsm{s} {self.s}"


class MaxLatticeDilatation(_S):
    def __init__(self, additional: float = 0, exceed: bool = False) -> None:
        super().__init__()
        assert _pos_num(additional), "Additional percentage of memory must be a positive number"
        if not exceed:
            assert additional <= 15, "Additional percentage of memory must be less than or equal to 15 (when exceed is False)"
        self.a = additional
        self.x = exceed
    
    def stamp(self, index: int):
        s = index or ''
        return f"dilatmx{s} {1 + self.a*0.01}\nchkdilatmx{s} {int(not self.x)}"


_ex1 = set((MolecularDynamics, StructuralOptimization, FIRE, MonteCarloSampling))
_ex2 = set((CellOptimization, CellOptimizationOnVector))


# check se non SCF che bande siano specificate
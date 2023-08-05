from ._common import Stampable, StampCollection, CanDelay, Delayed, Later, _pos_int, _pos0_int, DelayedInfo
from typing import Union, Tuple, Literal, Optional
from enum import Enum
from .units import Energy


__all__ = ["SpinType", "Metal", "Smearing", "Semiconductor", "TwoQuasiFermilevels", "OccupationPerBand"]


class SpinPolarization(Stampable):
    def __init__(self, pol: Literal[1,2], inor: Literal[1,2], den: Literal[1,2,4]) -> None:
        """DO NOT USE: use precomputed spin polazitation instead"""
        super().__init__()
        self.polarizationNumber = pol
        self.spinorNumber = inor
        self.density = den

    def stamp(self, index: int):
        s = index or ''
        return f"""nsppol{s} {self.polarizationNumber}
nspinor{s} {self.spinorNumber}
nspden{s} {self.density}"""
    

class SpinType:
    """Possible spin polarization options\n\nDO NOT MODIFY"""
    Unpolarized = SpinPolarization(1,1,1)
    Antiferromagnetic = SpinPolarization(1,1,2)
    SpinOrbitCoupling = SpinPolarization(1,2,1)
    NonCollinearMagnetism = SpinPolarization(1,2,4)
    Polarized = SpinPolarization(2,1,2)


class Smearing(Enum):
    FermiDirac = 3
    """
    Fermi-Dirac smearing for a finite-temperature metal.
    
    Smeared delta function: 0.5*(cosh(x/2))^(-2)
     
    For usual calculations, at zero temperature, do not use this, but prefer Gaussian smearing instead. If you want to do a calculation at finite temperature, please refer to the electron physical temperature.
    """

    Marzari5634 = 4
    """“Cold smearing” of N. Marzari (see his thesis work), with a=-0.5634 -> minimization of the bump"""

    Marzari8165 = 5
    """“Cold smearing” of N. Marzari (see his thesis work), with a=-0.8165 -> monotonic function in the tail """

    MethfesselPaxton = 6
    """Smearing of Methfessel and Paxton [Methfessel1989](https://docs.abinit.org/theory/bibliography#methfessel1989) with Hermite polynomial of degree 2, corresponding to “Cold smearing” of N. Marzari with a=0"""

    Gaussian = 7
    """
    Gaussian smearing, corresponding to the 0-order Hermite polynomial of Methfessel and Paxton. 
    
    Smeared delta function: exp(-x^2)/sqrt(π)
    
    Robust and quite efficient.
    """
    
    Uniform = 8
    """Uniform smearing: the delta function is replaced by a constant function of value one over ]-0.5,+0.5[ (with one-half value at the boundaries). Used for testing purposes only."""


_D_pos_int = DelayedInfo.basic(_pos_int, "must be a positive integer")
_D_bands = _D_pos_int("nband", "Number of bands")


class _D_pos_energy(DelayedInfo):
    def sanitize(self, value):
        e = Energy.sanitize(value)
        assert e._v > 0, f"{self.name} must be a positive energy (or float)"
        return e


class Metal(CanDelay):
    """Metallic occupation of levels, using different occupation schemes. All k points have the same number of bands. The combination of a broadening and a physical temperature can be obtained by using both tsmear and tphysel."""

    _delayables = ( _D_bands, _D_pos_energy("tsmear", "Smearing brodening temperature"), _D_pos_energy("tphysel", "Electrons' physical temperature"))

    def __init__(self, 
                 smearing: Smearing, 
                 bands: Union[int,Later] = Later(), 
                 broadening: Union[float,Energy,Later] = Later(), 
                 ePhysicalTemp: float|Energy|Later = Later()
                 ) -> None:
        """
         * `smearing`: the smerign function to use.
         * `bands`: number of bands, occupied plus possibly unoccupied, for which wavefunctions are being computed along with eigenvalues; must be doubled in case of 2 wavefunction spinorial component (`SpinType.SpinOrbitCoupling` and `Spintype.NonCollinearMagnetism`)
         * `broadening`: temperature of smearing which gives the broadening of occupation numbers. Abinit default value is 0.01 Hartree, which should be OK using gaussian like smearings (Marzari, Methfessel-Paxton, Gaussian) for a free-electron metal, like Al. For d-band metals, you may need to use less.
            * For the Fermi-Dirac smearing, it is the physical temperature, as the broadening is based on Fermi-Dirac statistics
            * For Gaussian-like smearing (Marzari, Methfessel-Paxton, Gaussian), it is only a convergence parameter, while the pysical temperature can be set using `ePhysicalTemp`.
         * `ePhysicalTemp`: the physical temperature of the system; note that the signification of the entropy is modified with respect to the usual entropy: the choice has been made to use `broadening` as a prefactor of the entropy, to define the entropy contribution to the free energy.
        """
        super().__init__(bands, broadening, ePhysicalTemp)
        self._opt = smearing.value

    def stamp(self, index: int):
        return f"occopt{index or ''} {self._opt}\n{super().stamp(index)}"

    @classmethod
    def setBands(cls, value: int):
        return Delayed(cls, 0, value)
    
    @classmethod
    def setBroadening(cls, value: Energy|float):
        return Delayed(cls, 1, value)
    
    @classmethod
    def setEPhysicalTemp(cls, value: Energy|float):
        return Delayed(cls, 2, value)

    

class Semiconductor(CanDelay):
    """All k points and spins have the same number of bands. All k points have the same occupancies of bands for a given spin (but these occupancies may differ for spin up and spin down - typical for ferromagnetic insulators), and they are automatically generated by the code to give a semiconductor."""

    _delayables = (_D_bands,)

    def __init__(self, spinMagnetizationTarget: Optional[float], bands: Union[int,Later] = Later()):
        """If the spin is polarized (you have used SpinType.Polarized), then `spinMagnetizationTarget` is mandatory"""
        super().__init__(bands)
        if spinMagnetizationTarget is not None and type(spinMagnetizationTarget) is not float:
            raise TypeError("spinMagnetizationTarget must be a float or None")
        self._smt = spinMagnetizationTarget
    
    def compatible(self, coll: StampCollection):
        spin = coll.get(SpinPolarization)
        if spin is not None and spin.polarizationNumber == 2:
            assert self._smt is not None, "Semiconducotor with polarized spin must specify the spin magnetization target"

    def stamp(self, index: int):
        res = [f"occopt{index or ''} 1"]
        if self._smt is not None:
            res.append(f"spinmagntarget{index or ''} {self._smt}")
        s = super().stamp(index)
        if len(s) > 0:
            res.append(s)
        return '\n'.join(res)

    @classmethod
    def setBands(cls, value: int):
        return Delayed(cls, 0, value)
    

class TwoQuasiFermilevels(CanDelay):
    """Fermi-Dirac occupation is enforced with two distinct quasi-Fermi levels:  See details in [Paillard2019]. At present, the number of holes and electrons should be the same. 
    
    Cannot be used with fixed magnetization calculation."""

    _delayables = (_D_bands,)

    def __init__(self, carriers: int, valenceBands: int, bands: Union[int,Later] = Later()):
        """`carriers` holes are forced in the first `valenceBands` and `carriers` electrons are forced in bands above."""
        super().__init__(bands)
        assert _pos_int(carriers), "Number of carriers must be a positive integer"
        assert _pos0_int(valenceBands), "Number of valenceBands must be a positive (or null) integer"
        self._c = carriers
        self._v = valenceBands

    def stamp(self, index: int):
        s = index or ''
        r = f"occopt{s} 9\nivalence{s} {self._v}\nnqfd{s} {self._c}"
        d = super().stamp(index)
        if len(d) > 0:
            r += d
        return r

    @classmethod
    def setBands(cls, value: int):
        return Delayed(cls, 0, value)
    

class OccupationPerBand(Stampable):
    def __init__(self, *occupations: Union[float, Tuple[float, float]], repeat: Optional[int] = None) -> None:
        super().__init__()
        assert all(
            type(v) is float or (type(v) is tuple and len(v) == 2 and type(v[0]) is float and type(v[1]) is float)
            for v in occupations
        ), "Occupations must be float or tuple of two floats"

        self._d = False
        self._o = occupations
        if repeat is None:
            self._r = 0
        elif type(repeat) is int:
            assert len(occupations) == 1, "If repeat is given, only one occupation must be given"
            self._r = repeat
        else:
            raise TypeError("OccupationPerBands 'repeat' argument must be None or positive int")
        
    def compatible(self, coll: StampCollection):
        spin = coll.get(SpinPolarization)
        if spin is None or spin.polarizationNumber == 1:
            self._d = False
            assert all(type(v) is float for v in self._o), "Since the spin is not polarized, ony one occupation per band must be given"
        else:
            self._d = True
    
    def stamp(self, index: int):
        o = ''
        b = 1
        if self._r == 0:
            b = len(self._o)
            if self._d:
                s1 = [0.0]*b
                s2 = [0.0]*b
                for i,o in enumerate(self._o):
                    if type(o) is tuple:
                        s1[i] = o[0]
                        s2[i] = o[1]
                    else:
                        s1[i] = s2[i] = o # type: ignore
                o = ' '.join(str(v) for v in s1) + ' ' + ' '.join(str(v) for v in s2)
            else:
                o = ' '.join(str(o) for o in self._o)
        else:
            b = self._r
            if self._d:
                if type(self._o[0]) is tuple:
                    o = f"{self._r}*{self._o[0][0]} {self._r}*{self._o[0][1]}"
                else:
                    o = f"{self._r*2}*{self._o[0]}"
            else:
                o = f"{self._r}*{self._o[0]}"
        s = index or ''
        return f"occopt{s} 0\nocc{s} {o}\nband{s} {b}"
    

_exclusives = set((Metal, Semiconductor, TwoQuasiFermilevels, OccupationPerBand))
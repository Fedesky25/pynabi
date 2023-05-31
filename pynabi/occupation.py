from typing import Union, Tuple, Literal
from enum import Enum
from ._common import Stampable


__all__ = ["SpinType", "Occupation", "Smearing"]


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
    Marzari5634 = 4
    Marzari8165 = 5
    MethfesselPaxton = 6
    Gaussian = 7
    Uniform = 8


class Occupation(Stampable):
    def __init__(self, option: int, bands: Union[None,int,Tuple[int]], **props) -> None:
        super().__init__()
        self.option = option
        self.bands = bands
        self.props = props

    def stamp(self, index: int):
        s = index or ''
        extra = "".join(f"\n{n}{s} {v}" for (n,v) in self.props.items())
        return f"occopt{s} {self.option}{_band2str(self.bands, s)}{extra}"
    
    @staticmethod
    def EqualBandNumber(occupations: float, bands: int, polarized: bool = False):
        """DO NOT USE THIS METHOD unless you know what you are doing\n\nPlease use the related preset function"""
        m = bands * (1 + int(polarized))
        return Occupation(0, bands, occ=f"{m}*{occupations}")

    @staticmethod
    def Semiconductor(bands: Union[int,None] = None, spinMagnetizationTarget: Union[float,None] = None):
        if spinMagnetizationTarget is None:
            return Occupation(1, bands)
        else:
            return Occupation(1, bands, spinmagntarget=spinMagnetizationTarget)
    
    # @staticmethod
    # def Manual(*occupations: Tuple[float]):
    #     """DO NOT USE THIS METHOD unless you know what you are doing\n\nPlease use the related preset function"""
    #     bands = tuple(len(o) for o in occupations)
    #     return Occupation(2, bands, occ=occupations)

    @staticmethod
    def Metallic(smearing: Smearing, broadening: float, bands: Union[int,None] = None):
        return Occupation(smearing.value, bands, tsmear=broadening)
    
    @staticmethod
    def TwoQuasiFermiLevels(carriers: int, valenceIndex: int, bands: Union[int,None] = None):
        return Occupation(9, bands, nqfd=carriers, ivalence=valenceIndex)
    

def _band2str(b: Union[None,int,Tuple[int]], s):
    if b is None:
        return ''
    s = f"\nnbands{s} "
    if type(b) is tuple:
        s = s + (' '.join(str(v) for v in b))
    else:
        s = s + str(b)
    return s
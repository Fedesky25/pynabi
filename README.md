# PynAbi

Python package to easily create [Abinit](https://www.abinit.org/) input files.

```cmd
pip install pynabi==0.0.3
```

## Example

```python
from pynabi import createAbi, DataSet, AbIn, AbOut
from pynabi.kspace import CriticalPointsOf, BZ, SymmetricGrid, UsualKShifts, Path
from pynabi.calculation import ToleranceOn, EnergyCutoff, MaxSteps, SCFMixing, NonSelfConsistentCalc
from pynabi.crystal import Atom, FluoriteLike
from pynabi.occupation import OccupationPerBand
from pynabi.units import eV, nm

# create manually an atom -> Atom(<Z>, <pseudo potential name>)
# or using sensible defaults as follows
Zr = Atom.of("Zr")  # Z=40 and pseudos located at "Zr.psp8"
Oxy = Atom.of("O")  # Z=8 and pseudos located at "O.psp8"

# base dataset with common variables
base = DataSet(
    AbOut("./scf/scf"),                             # prefix for output files
    AbIn().PseudoPotentials("./pseudos/PBE-SR"),    # folder with pseudo potentials

    FluoriteLike(Zr, Oxy, 0.5135*nm),               # creates AtomBasis and Lattice of a crystal like fluorite
                                                    # with lattice constant 0.5135nm

    SymmetricGrid(BZ.Irreducible, UsualKShifts.FCC) # easily define kptopt, ngkpt, nshiftk, kpt
        .ofMonkhorstPack(4),

    SCFMixing(density=True).Pulay(10),              # scf cycle with Pulay mixing of the density 
                                                    # based on the last 10 iteration

    ToleranceOn.EnergyDifference(1e-6),             # expressively define the tolerance
    MaxSteps(30)                                    # nstep
)

# set the default energy unit in eV (from now on)
eV.setAsReference()

# datasets to see the convergenge as a function of energy
sets = [DataSet(EnergyCutoff(8.0 + i*0.25)) for i in range(0,17)]

# final non-self-consistent round to find bands 
bands = DataSet(
    NonSelfConsistentCalc(),
    ToleranceOn.WavefunctionSquaredResidual(1e-12),
    AbIn().ElectronDensity(sets[-1]),                   # get the electron density from the last dataset
    OccupationPerBand(2.0, repeat=8),                   # same number of bands (max 8) for each k point
    Path.auto(10, "GXWKGLUWLK", CriticalPointsOf.FCC)   # easily define a path in the k-space   
)

with open("./out.txt", 'w') as f:
    f.write(createAbi(base, *sets, bands))
```

<details>
<summary><b>Output</b></summary>

```txt
ndtset 18

# Atoms definition
ntypat 2
znucl 8 40
pseudos "O.psp8, Zr.psp8"

# Common DataSet
natoms 3
typeat 2 1 1
xred 0 0 0   0.3333333333333333 0.3333333333333333 0.3333333333333333   0.6666666666666666 0.6666666666666666 0.6666666666666666
outdata_prefix "./scf/scf"
pp_dirpath "./pseudos/PBE-SR"
scalecart 0.5135 0.5135 0.5135 nm
rprim 0.5 0.5 0.0   0.0 0.5 0.5   0.5 0.0 0.5
ngkpt 4 4 4
iscf 17
npulayit 10
toldfe 1e-06
nstep 30

# DataSet 1
ecut1 8.0 eV

# DataSet 2
ecut2 8.25 eV

# DataSet 3
ecut3 8.5 eV

# DataSet 4
ecut4 8.75 eV

# DataSet 5
ecut5 9.0 eV

# DataSet 6
ecut6 9.25 eV

# DataSet 7
ecut7 9.5 eV

# DataSet 8
ecut8 9.75 eV

# DataSet 9
ecut9 10.0 eV

# DataSet 10
ecut10 10.25 eV

# DataSet 11
ecut11 10.5 eV

# DataSet 12
ecut12 10.75 eV

# DataSet 13
ecut13 11.0 eV

# DataSet 14
ecut14 11.25 eV

# DataSet 15
ecut15 11.5 eV

# DataSet 16
ecut16 11.75 eV

# DataSet 17
ecut17 12.0 eV

# DataSet 18
iscf18 -2
tolwfr18 1e-12
getden18 17
occopt18 0
occ18 8*2.0
bands18 8
kptopt18 -9
kptbounds18 0 0 0  0.0 0.5 0.5  0.25 0.75 0.5  0.375 0.75 0.375  0 0 0  0.5 0.5 0.5  0.25 0.625 0.625  0.25 0.75 0.5  0.5 0.5 0.5  0.375 0.75 0.375
ndivsm18 10
```

</details>

## Features

 - Multi dataset support
 - Helper functions for common crystal structures (caesium chloride, rock-salt, fluorite, zincblende, wurtzite, nickeline, HCP)
 - Registered critical points of and methodsto create lattices of CUB, BCC, FCC, HEX, TET, BCT, ORC, ORCC
 - Smooth experience in defining the k-points
 - Handy management of [file handling variables](https://docs.abinit.org/variables/files/)
 - Almost full covarage of [basic input variables](https://docs.abinit.org/variables/basic/) (missing nbandhf, symrel, tnons, wvl_hgrid)
 - Partial coverage of [ground state variables](https://docs.abinit.org/variables/gstate/)
 - Partial coverage of [relaxation variables](https://docs.abinit.org/variables/rlx/)
 - _More to come..._

## Why PynAbi over pure Abinit files?

If you're a very experienced Abinit user who knows its variables and their possible values (and associated meaning), then this package probably isn't for you.
For all other users, on top of all the aforementioned features, here's a list of reasons why you could find PynAbi useful:

 1. You have all the power of coding to generate Abinit instructions, e.g. reusability, loops to generate datasets programmatically
 2. It provides some useful presets and helper functions/methods that allows you to skip to the fun part of the simulation
 3. Under the hood, it checks for the validity of variable values *before* starting Abinit
 4. It makes use of expressive declarations and definition, leading to readable and comprehensible istructions
 5. If you're using a code editor (with autocompletition), you'll get suggetions of all the possible options in a more natural language 

Some prior knowledge of Abinit (and DFT in general) is nonetheless needed to fully understand what actually is to be simulated: this package only provides a handier way to generate the required files.

## Documentation

Although currently severely incomplete, a wiki is available on the [GitHub repository Wiki](https://github.com/Fedesky25/pynabi/wiki). Ideally a wiki for this package shouldn't be needed if your editor has some autocompletition: the aim is to document every class and function with a detailed description taken from the relevant parts of ABINIT documentation. 
# PynAbi

Python package to easily create [Abinit](https://www.abinit.org/) input files.

## Example

```python
from pynabi import createAbi, DataSet, AbIn, AbOut, Occupation
from pynabi.kspace import CriticalPointsOf, BZ, SymmetricGrid, path, UsualKShifts
from pynabi.calculation import ToleranceOn, EnergyCutoff, MaxSteps, SCFMixing, NonSelfConsistentCalc
from pynabi.crystal import Atom, BCC

# create manually an atom 
# or using sensible defaults as follows
Si = Atom.of("Si") # Z=14 and pseudos located at "Si.psp8"

# base dataset with common variables
base = DataSet(
    AbOut("./scf/scf"),                             # prefixes for output files
    AbIn().PseudoPotentials("./pseudos/PBE-SR"),    # folder with pseudo potentials
    BCC(5.09, Si, Si),                              # creates AtomBasis and Lattice of a BCC crystal
    SymmetricGrid(BZ.Irreducible, UsualKShifts.BCC) # easily define kptopt, ngkpt, nshiftk, kpt
        .ofMonkhorstPack(4),
    SCFMixing(density=True).Pulay(10),              # scf cycle with Pulay mixing of the density based on the last 10 iteration
    ToleranceOn.EnergyDifference(1e-6),             # expressively define the tolerance
    MaxSteps(30)                                    # nstep
)

# datasets to see the convergenge as a function of energy
sets = [DataSet(EnergyCutoff(8.0 + i*0.25)) for i in range(0,17)]

# final non-self-consistent round to find bands 
bands = DataSet(
    NonSelfConsistentCalc(),
    ToleranceOn.WavefunctionSquaredResidual(1e-12),
    AbIn().ElectronDensity(sets[-1]),               # get the electron density from the last dataset
    Occupation.Semiconductor(bands=8),              # semiconductor occupation (occopt=1) with 8 bands
    path(10, "GHNGPH", CriticalPointsOf.BCC),       # easily define a path in the k-space   
)

with open("./out.txt", 'w') as f:
    f.write(createAbi(base, *sets, bands))
```

<details>
<summary><b>Output</b></summary>

```txt
ndtset 18

# Atoms definition
ntypat 1
znucl 14
pseudos "Si.psp8"

# Common DataSet
natoms 2
typeat 1 1
xred 0 0 0   0.5 0.5 0.5
outdata_prefix "./scf/scf"
pp_dirpath "./pseudos/PBE-SR"
acell 5.09 5.09 5.09
angdeg 90 90 90
kptopt 1
nshiftk 2
shiftk 0.25 0.25 0.25   -0.25 -0.25 -0.25
ngkpt 4 4 4
iscf 17
npulayit 10
toldfe 1e-06
nstep 30

# DataSet 1
ecut1 8.0 Hartree

# DataSet 2
ecut2 8.25 Hartree

# DataSet 3
ecut3 8.5 Hartree

# DataSet 4
ecut4 8.75 Hartree

# DataSet 5
ecut5 9.0 Hartree

# DataSet 6
ecut6 9.25 Hartree

# DataSet 7
ecut7 9.5 Hartree

# DataSet 8
ecut8 9.75 Hartree

# DataSet 9
ecut9 10.0 Hartree

# DataSet 10
ecut10 10.25 Hartree

# DataSet 11
ecut11 10.5 Hartree

# DataSet 12
ecut12 10.75 Hartree

# DataSet 13
ecut13 11.0 Hartree

# DataSet 14
ecut14 11.25 Hartree

# DataSet 15
ecut15 11.5 Hartree

# DataSet 16
ecut16 11.75 Hartree

# DataSet 17
ecut17 12.0 Hartree

# DataSet 18
iscf18 -2
tolwfr18 1e-12
getden18 17
occopt18 1
nbands18 8
kptopt18 -5
kptbounds18 0 0 0   -0.5 0.5 0.5   0.0 0.5 0.0   0 0 0   0.25 0.25 0.25   -0.5 0.5 0.5
ndivsm18 10
```

</details>

## Features

 - Multi dataset support
 - Helper funtion/methods for common crystal structures (CUB, BCC, FCC, HEX, RHL)
 - Smooth experience in defining the k-points
 - Handy management of [file handling variables](https://docs.abinit.org/variables/files/)
 - Almost full covarage of [basic input variables](https://docs.abinit.org/variables/basic/) (missing nbandhf, symrel, tnons, wvl_hgrid)
 - Partial coverage of [ground state variables](https://docs.abinit.org/variables/gstate/)
 - _More to come..._

## Why PynAbi over pure Abinit files?

If you're a very experienced Abinit user who knows its variables and their possible values (and associated meaning), then this package probably isn't for you.
For all other users, here's a list of reasons why you could find PynAbi useful:

 1. It allows for programmatic use of the datasets, e.g. reusability of instuctions, loop to generate datasets programmatically
 2. It provides some useful presets and helper functions/methods that allows you to skip to the fun part of the simulation
 3. Under the hood, it checks for the validity of variable values *before* starting Abinit
 4. It makes use of expressive declarations and definition, leading to readable and comprehensible istructions
 5. If you're using a code editor (with autocompletition), you'll get suggetions of all the possible options in a more natural language 

Some prior knowledge of Abinit (and DFT in general) is nonetheless needed to fully understand what actually is to be simulated: this package only provides a handier way to generate the required files.

## Documentation

_coming soon_
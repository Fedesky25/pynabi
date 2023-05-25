# PynAbi

Python package to easily create [Abinit](https://www.abinit.org/) input files.

## Example

```python
from pynabi import createAbi, DataSet, presets, Atom, BZ, SymmetricKGrid, KPath, Vec3D, ToleranceOn, EnergyCutoff, StepNumber, SCFProcedure, AbIO, I, Occupation

# folder with pseudo potentials
pseudo_folder = "./pseudos/PBE-SR"

# create manually an atom 
# or using sensible defaults as follows
Si = Atom.of("Si") # Z=14 and pseudos located at "Si.psp8"

# base dataset with common variables
base = DataSet(
    AbIO().prefix(output="./scf/scf"),          # prefixes for input, output, temporary files
    presets.BCC(5.09, Si, Si, pseudo_folder),   # creates AtomBasis and Lattice of a BCC
    SymmetricKGrid(BZ.Irreducible, 2,           # easily define kptopt, ngkpt, nshiftk, kpt
                   Vec3D(0.5, 0.5, 0.5),
                   Vec3D(0.5, 0.0, 0.0),
                   Vec3D(0.0, 0.5, 0.0),
                   Vec3D(0.0, 0.0, 0.5)),
    ToleranceOn.EnergyDifference(1e-6),         # expressively define the tolerance
    SCFProcedure(0),                            # iscf
    StepNumber(30)                              # nstep
)

# datasets to see the convergenge as a function of energy
sets = [DataSet(EnergyCutoff(8.0 + i*0.25)) for i in range(0,17)]

# final non-self-consistent round to find bands 
bands = DataSet(
    SCFProcedure(-2),
    ToleranceOn.WavefunctionSquaredResidual(1e-12),
    AbIO().get(I.ElectronDensity, sets[-1]),            # get the elctron density from the last dataset
    Occupation.Semiconductor(bands=8),                  # semiconductor occupation (occopt=1) with 8 bands
    KPath.through(10, CriticalPointsOf.BCC, "GHNGPH"),  # define a path in the k-space
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
xred 0 0 0
      0.5 0.5 0.5
pp_dirpath "./pseudos/PBE-SR"
outdata_prefix "./output/prefix"
acell 5.09 5.09 5.09
angdeg 90 90 90
kptopt 1
ngkpt 2 2 2
nshiftk 4
kpt 0.5 0.5 0.5
    0.5 0.0 0.0
    0.0 0.5 0.0
    0.0 0.0 0.5
toldfe 1e-06
iscf 0
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
kptbounds18 0 0 0
            -0.5 0.5 0.5
            0.0 0.5 0.0
            0 0 0
            0.25 0.25 0.25
            -0.5 0.5 0.5
ndivsm18 10
```
</details>

## Features

 - Presets for common crystal structures (CUB, BCC, FCC, HCP)
 - Helper functions to better control the k-points definition
 - Multi dataset support
 - Handy management of [file handling variables](https://docs.abinit.org/variables/files/)
 - Almost full covarage of [basic input variables](https://docs.abinit.org/variables/basic/)
 - _More to come..._

## Documentation

_coming soon_
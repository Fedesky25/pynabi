import sys, os
sys.path.append(os.getcwd()+"\\src")
# TODO: find the correct way to import pynabi for testing purposes
# ================================================================


from pynabi import createAbi, DataSet, AbOut, AbIn
from pynabi.crystal import ZincBlendeLike, Atom
from pynabi.calculation import EnergyCutoff, ToleranceOn
from pynabi.kspace import BZ, UsualKShifts, SymmetricGrid

Si = Atom.of("Si")

base = DataSet(
    AbOut("./scf/scf"),
    AbIn().PseudoPotentials("./pseudos/PBE-SR"),
    ZincBlendeLike(Si, Si, 7.31),
    SymmetricGrid(BZ.Irreducible,UsualKShifts.FCC)
        .ofMonkhorstPack(4),
    EnergyCutoff(12.0),
    ToleranceOn.EnergyDifference(1e-6)
)

# study convergenge w.r.t. energy cutoff
E_conv = [DataSet(EnergyCutoff(8.0 + i)) for i in range(0,9)]

# study convergenge w.r.t. number of k points
K_conv = [DataSet(SymmetricGrid.setMPgridPointNumber(i)) for i in range(4,10)]

with open("examples/2.txt", 'w') as file:
    file.write(createAbi(base, *E_conv, *K_conv))
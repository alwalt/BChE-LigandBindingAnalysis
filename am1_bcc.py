#!/usr/bin/env python
################################################################
# Copyright (C) 2015 OpenEye Scientific Software, Inc.
################################################################
# Generates canonical AM1-BCC charges
################################################################

import sys
from openeye.oechem import *
from openeye.oeomega import *
from openeye.oequacpac import *

def main(argv=[__name__]):
    if len(argv) != 3:
	OEThrow.Usage("%s <infile> <outfile>" % argv[0])

    ifs = oemolistream()
    if not ifs.open(argv[1]):
        OEThrow.Fatal("Unable to open %s for reading" % argv[1])

    if not OEIs3DFormat(ifs.GetFormat()):
        OEThrow.Fatal("Invalid input format: need 3D coordinates")

    ofs = oemolostream()
    if not ofs.open(argv[2]):
        OEThrow.Fatal("Unable to open %s for writing" % argv[2])

    if ofs.GetFormat() not in [OEFormat_MOL2, OEFormat_OEB]:
        OEThrow.Error("MOL2 or OEB output file is required!")

    omega = OEOmega()
    omega.SetIncludeInput(True)
    omega.SetCanonOrder(False)
    omega.SetSampleHydrogens(True)
    eWindow = 15.0
    omega.SetEnergyWindow(eWindow)
    omega.SetMaxConfs(1000)
    omega.SetRMSThreshold(1.0)

    for mol in ifs.GetOEMols():
        if omega(mol):
            OEAssignPartialCharges(mol, OECharges_AM1BCCSym)
            conf = mol.GetConf(OEHasConfIdx(0))
            absFCharge = 0
            sumFCharge = 0
            sumPCharge = 0.0
            for atm in mol.GetAtoms():
                sumFCharge += atm.GetFormalCharge()
                absFCharge += abs(atm.GetFormalCharge())
                sumPCharge += atm.GetPartialCharge()
            OEThrow.Info("%s: %d formal charges give total charge %d ; Sum of Partial Charges %5.4f"
                         % (mol.GetTitle(), absFCharge, sumFCharge, sumPCharge))
            OEWriteMolecule(ofs, conf)
        else:
            OEThrow.Warning("Failed to generate conformation(s) for molecule %s" % mol.GetTitle())

    return 0

if __name__ == "__main__":
    sys.exit(main(sys.argv))

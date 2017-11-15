#!/usr/bin/env python

#######################################################################################################
## Generates canonical AM1-BCC charges and assigns generic tags to partial charges for visualization ##
#######################################################################################################

import sys
from openeye.oechem import *
from openeye.oeomega import *
from openeye.oequacpac import *

def main(argv=[__name__]):
    
    itf = OEInterface(InterfaceData)

    if not OEParseCommandLine(itf, argv):
        return 1

    iname = itf.GetString("-in")
    oname = itf.GetString("-out")

    ifs = oemolistream()
    if not ifs.open(argv[1]):
        OEThrow.Fatal("Unable to open %s for reading" % argv[0])

    if not OEIs3DFormat(ifs.GetFormat()):
        OEThrow.Fatal("Invalid input format: need 3D coordinates")

    ofs = oemolostream()
    if not ofs.open(argv[2]):
        OEThrow.Fatal("Unable to open %s for writing" % argv[2])

    if ofs.GetFormat() not in [OEFormat_MOL2, OEFormat_OEB]:
        OEThrow.Error("MOL2 or OEB output file is required!")

    tagname = itf.GetString("-tagname")

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

 	    tag = OEGetTag(tagname)

            for atm in mol.GetAtoms():
                sumFCharge += atm.GetFormalCharge()
                absFCharge += abs(atm.GetFormalCharge())
                sumPCharge += atm.GetPartialCharge()
		atm.SetData(tag, atm.GetPartialCharge())
            OEThrow.Info("%s: %d formal charges give total charge %d ; Sum of Partial Charges %5.4f"
                         % (mol.GetTitle(), absFCharge, sumFCharge, sumPCharge))
            OEWriteMolecule(ofs, conf)
        else:
            OEThrow.Warning("Failed to generate conformation(s) for molecule %s" % mol.GetTitle())

    return 0

InterfaceData = """
!CATEGORY "input/output options"

  !PARAMETER -in
    !ALIAS -i
    !TYPE string
    !REQUIRED true
    !KEYLESS 1
    !VISIBILITY simple
    !BRIEF Input filename
  !END

  !PARAMETER -out
    !ALIAS -o
    !TYPE string
    !REQUIRED true
    !KEYLESS 2
    !VISIBILITY simple
    !BRIEF Output OEB filename
  !END

!CATEGORY "general options"

  !PARAMETER -tagname
    !ALIAS -tag
    !TYPE string
    !REQUIRED true
    !KEYLESS 3
    !VISIBILITY simple
    !BRIEF Generic data tag name for atom data
  !END

!END
"""

if __name__ == "__main__":
    sys.exit(main(sys.argv))


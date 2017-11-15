#!/usr/bin/env python

####################################################
##Visualizes generic atom properties read from OEB##
####################################################

import sys
from openeye.oechem import *
from openeye.oedepict import *
from openeye.oegrapheme import *


def main(argv=[__name__]):

    itf = OEInterface(InterfaceData)
    OEConfigureImageWidth(itf, 400.0)
    OEConfigureImageHeight(itf, 400.0)
    OEConfigure2DMolDisplayOptions(itf, OE2DMolDisplaySetup_AromaticStyle)
    OEConfigureColor(itf, "-negcolor", "-nc", "Color for negative values", "red")
    OEConfigureColor(itf, "-poscolor", "-pc", "Color for positive values", "blue")

    if not OEParseCommandLine(itf, argv):
        return 1

    iname = itf.GetString("-in")
    oname = itf.GetString("-out")
    style = itf.GetString("-style")

    # check input/output files

    ifs = oemolistream()
    if not ifs.open(iname):
        OEThrow.Fatal("Cannot open input file!")

    if ifs.GetFormat() != OEFormat_OEB:
        OEThrow.Fatal("Only works for oeb input file!")

    ext = OEGetFileExtension(oname)
    if not OEIsRegisteredImageFile(ext):
        OEThrow.Fatal("Unknown image type!")

    ofs = oeofstream()
    if not ofs.open(oname):
        OEThrow.Fatal("Cannot open output file!")

    # read a molecule

    mol = OEGraphMol()
    if not OEReadMolecule(ifs, mol):
        OEThrow.Fatal("Cannot read input file!")

    # check atom properties

    tagname = itf.GetString("-tagname")
    if not CheckAtomProperties(mol, tagname):
        OEThrow.Error("Cannot find tag %s on input molecule" % tagname)

    # prepare depiction

    clearcoords = itf.GetBool("-clearcoords")
    suppressH = itf.GetBool("-suppressH")
    OEPrepareDepiction(mol, clearcoords, suppressH)

    # create image

    width, height = OEGetImageWidth(itf), OEGetImageHeight(itf)
    image = OEImage(width, height)

    # setup depiction options

    opts = OE2DMolDisplayOptions(width, height, OEScale_AutoScale)
    OESetup2DMolDisplayOptions(opts, itf)
    opts.SetAtomColorStyle(OEAtomColorStyle_WhiteMonochrome)

    negcolor = OEGetColor(itf, "-negcolor")
    poscolor = OEGetColor(itf, "-poscolor")

    DepictAtomProperty(image, mol, opts, tagname, negcolor, poscolor, style)

    OEWriteImage(oname, image)

    return 0


def DepictAtomProperty(image, mol, opts, tagname, negcolor, poscolor, style):

    mwidth, mheight = image.GetWidth(), image.GetHeight() * 0.9
    cwidth, cheight = image.GetWidth(), image.GetHeight() * 0.1

    mframe = OEImageFrame(image, mwidth, mheight, OE2DPoint(0.0, 0.0))
    cframe = OEImageFrame(image, cwidth, cheight, OE2DPoint(0.0, mheight))

    opts.SetDimensions(mwidth, mheight, OEScale_AutoScale)
    opts.SetScale(OEGetMoleculeSurfaceScale(mol, opts))

    colorg = GetColorGradient(mol, tagname, negcolor, poscolor)

    disp = OE2DMolDisplay(mol, opts)

    if style == "propmap":
        DepictAtomPropertyPropMap(disp, tagname, negcolor, poscolor)

    if style == "atomglyph":
        DepictAtomPropertyAtomGlyph(disp, tagname, colorg)

    if style == "molsurface":
        DepictAtomPropertyMolSurface(disp, tagname, colorg)

    OERenderMolecule(mframe, disp)

    OEDrawColorGradient(cframe, colorg)

    font = OEFont(OEFontFamily_Default, OEFontStyle_Default, 14, OEAlignment_Left, OEBlack)
    cframe.DrawText(OE2DPoint(10.0, -10.0), tagname, font)


def CheckAtomProperties(mol, tagname):

    tag = OEGetTag(tagname)

    for atom in mol.GetAtoms():
        if atom.HasData(tag):
            return True
    return False


def GetMinMaxAtomProperty(mol, tagname):

    minvalue = float("inf")
    maxvalue = float("-inf")

    tag = OEGetTag(tagname)

    for atom in mol.GetAtoms():
        if atom.HasData(tag):
            val = atom.GetData(tag)
            minvalue = min(minvalue, val)
            maxvalue = max(maxvalue, val)

    return minvalue, maxvalue


def GetColorGradient(mol, tagname, ncolor, pcolor):

    minvalue, maxvalue = GetMinMaxAtomProperty(mol, tagname)

    colorg = OELinearColorGradient(OEColorStop(0.0, OEWhite))
    if minvalue < 0.0:
        colorg.AddStop(OEColorStop(minvalue, ncolor))
    if maxvalue > 0.0:
        colorg.AddStop(OEColorStop(maxvalue, pcolor))

    return colorg


def DepictAtomPropertyPropMap(disp, tagname, negcolor, poscolor):

    opts = disp.GetOptions()
    propmap = OE2DPropMap(opts.GetBackgroundColor())
    propmap.SetLegendLocation(OELegendLocation_Hidden)
    propmap.SetNegativeColor(negcolor)
    propmap.SetPositiveColor(poscolor)
    propmap.Render(disp, tagname)


def DepictAtomPropertyAtomGlyph(disp, tagname, colorg):

    tag = OEGetTag(tagname)
    mol = disp.GetMolecule()

    for atom in mol.GetAtoms():
        if atom.HasData(tag):
            value = atom.GetDoubleData(tag)
            color = colorg.GetColorAt(value)
            pen = OEPen(color, color, OEFill_Off, 3.0)
            glyph = OEAtomGlyphCircle(pen, OECircleStyle_Default, 1.2)
            OEAddGlyph(disp, glyph, OEHasAtomIdx(atom.GetIdx()))


def DepictAtomPropertyMolSurface(disp, tagname, colorg):

    tag = OEGetTag(tagname)
    mol = disp.GetMolecule()

    for atom in mol.GetAtoms():
        if atom.HasData(tag):
            value = atom.GetDoubleData(tag)
            color = colorg.GetColorAt(value)
            pen = OEPen(color, color, OEFill_Off, 4.0)
            OESetSurfaceArcFxn(mol, atom, OEDefaultArcFxn(pen))

    OEDraw2DSurface(disp)


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
    !BRIEF Output filename
  !END

!END

!CATEGORY "general options"

  !PARAMETER -tagname
    !ALIAS -tag
    !TYPE string
    !REQUIRED true
    !KEYLESS 3
    !VISIBILITY simple
    !BRIEF Generic data tag name for atomic data.
  !END

  !PARAMETER -dispstyle
    !ALIAS -style
    !TYPE string
    !REQUIRED false
    !DEFAULT propmap
    !VISIBILITY simple
    !LEGAL_VALUE atomglyph
    !LEGAL_VALUE propmap
    !LEGAL_VALUE molsurface
    !BRIEF Display style
    !DETAIL
        atomglyph - atom properties visualized by using atom glyphs
        propmap - atom properties visualized by property map
        molsurface - atom properties visualized on molecule surface
    !END

  !PARAMETER -clearcoords
    !ALIAS -clear
    !TYPE bool
    !REQUIRED false
    !DEFAULT false
    !VISIBILITY simple
    !BRIEF Clear 2D coordinates of input structure
  !END

  !PARAMETER -suppressH
    !ALIAS -sh
    !TYPE bool
    !REQUIRED false
    !DEFAULT false
    !VISIBILITY simple
    !BRIEF Suppress explicit hydrogens of input structure
    !DETAIL
      Explicit hydrogens that are necessary to represent tetrahedral stereochemistry are kept
  !END

!END

"""

if __name__ == "__main__":
    sys.exit(main(sys.argv))

#!/usr/bin/python

# foldxyastools.py

import yasara,string,disk,os,re
from python2to3 import *
from collections import namedtuple

# Zoom to mutation site and show sidechains of neighbourhood
# ==========================================================
def Zoom(residue,objnumbers,repairobjnumber,hbonds,vdwclashes):
  # make a list of all active objects to show hbonds
  if repairobjnumber != 0:
    objnumbers = str(residue.object.number.inyas) + " " + objnumbers + " " + str(repairobjnumber)
  else:
    objnumbers = str(residue.object.number.inyas) + " " + objnumbers
  # make a real list
  objnumberlist = objnumbers.split()
  # list residues that are in vicinity of wt or mt residue site
  contactreslist = []
  for objnum in objnumberlist:
    contactreslist += yasara.ListRes("Element !H Obj "+objnum+ " with distance < 5 from Element !H Res " + residue.number.inpdb + " Mol "+residue.molecule.name + " Obj "+objnum,"RESNUM MOLNAME")
  # plot all sidechains and CA that have anything to do with the mutation site (wt and mt) and plot them in both wt and mt to see all differences
  for contactres in contactreslist:
    contactresnumber = contactres.split()[0]
    # detect nameless molecules (these have a space and so [1] is not possible)
    try:
      contactresmolname = contactres.split()[1]
      yasara.ShowAtom("Sidechain CA Res "+contactresnumber+" Mol "+contactresmolname+ " Obj "+objnumbers)
    except:
      yasara.ShowAtom("Sidechain CA Res "+contactresnumber+" Obj "+objnumbers)
  # if we added hydrogen/bonds, hide the hydrogens
  if hbonds == 1:
    yasara.HideAtom("Element H")
  # zoom to mutation spot
  yasara.ZoomAtom("CA Res "+residue.number.inpdb+" Mol "+residue.molecule.name + " Obj "+objnum ,50,"Yes")
  # show message
  if hbonds == 1 and vdwclashes == 0:
    yasara.ShowMessage("Modeling complete. Visualized hydrogen bonds and zoomed in on mutation site.")
  elif hbonds == 1 and vdwclashes == 1:
    yasara.ShowMessage("Modeling complete. Visualized hydrogen bonds, steric clashes and zoomed in on mutation site.")
  elif vdwclashes == 1 and hbonds == 0:
    yasara.ShowMessage("Modeling complete. Visualized steric clashes and zoomed in on mutation site.")
  else:
    yasara.ShowMessage("Modeling complete. Zoomed in on mutation site.")
  yasara.Wait(3,"Seconds")
  yasara.HideMessage()

# Zoom to mutation site and show sidechains of neighbourhood after Mutate Multiple Residues
# =========================================================================================
def Zoom_MutateResidueStretch(residues,objnumbers,repairobjnumber,hbonds,vdwclashes):
  # make a list of all active objects to show hbonds
  if repairobjnumber != 0:
    objnumbers = str(residues.residue[0].object.number.inyas) + " " + objnumbers + " " + str(repairobjnumber)
  else:
    objnumbers = str(residues.residue[0].object.number.inyas) + " " + objnumbers
  # make a real list
  objnumberlist = objnumbers.split()
  # make a string with residuenumbers
  resnumbers = ""
  for i in range(residues.residues):
    resnumbers += str(residues.residue[i].number.inpdb) + " "
  # list residues that are in vicinity of wt or mt residue site
  contactreslist = []
  for objnum in objnumberlist:
    contactreslist += yasara.ListRes("Element !H Obj "+objnum+ " with distance < 5 from Element !H Res " + resnumbers + " Mol "+residues.residue[0].molecule.name + " Obj "+objnum,"RESNUM MOLNAME")
  # plot all sidechains and CA that have anything to do with the mutation site (wt and mt) and plot them in both wt and mt to see all differences
  for contactres in contactreslist:
    contactresnumber = contactres.split()[0]
    # detect nameless molecules (these have a space and so [1] is not possible)
    try:
      contactresmolname = contactres.split()[1]
      yasara.ShowAtom("Sidechain CA Res "+contactresnumber+" Mol "+contactresmolname+ " Obj "+objnumbers)
    except:
      yasara.ShowAtom("Sidechain CA Res "+contactresnumber+" Obj "+objnumbers)
  # if we added hydrogen/bonds, hide the hydrogens
  if hbonds == 1:
    yasara.HideAtom("Element H")
  # zoom to mutation spot
  yasara.ZoomRes(resnumbers+" Mol "+residues.residue[0].molecule.name + " Obj "+objnum ,50,"Yes")
  # show message
  if hbonds == 1 and vdwclashes == 0:
    yasara.ShowMessage("Modeling complete. Visualized hydrogen bonds and zoomed in on mutation site.")
  elif hbonds == 1 and vdwclashes == 1:
    yasara.ShowMessage("Modeling complete. Visualized hydrogen bonds, steric clashes and zoomed in on mutation site.")
  elif vdwclashes == 1 and hbonds == 0:
    yasara.ShowMessage("Modeling complete. Visualized steric clashes and zoomed in on mutation site.")
  else:
    yasara.ShowMessage("Modeling complete. Zoomed in on mutation site.")
  yasara.Wait(3,"Seconds")
  yasara.HideMessage()

# show Van der Waals clashes made by new residue
# ==============================================
def ShowVdwclashes(residue,objnumbers,repairobjnumber):
  # make a list of all active objects to show VdW clashes
  if repairobjnumber != 0:
    objnumbers = str(residue.object.number.inyas) + " " + objnumbers + " " + str(repairobjnumber)
  else:
    objnumbers = str(residue.object.number.inyas) + " " + objnumbers
  # make a real list
  objnumberlist = objnumbers.split()
  # loop over every object to show mutation site neighbourhood and VdW clashes between Carbon atoms (only C's for now)
  for objnum in objnumberlist:
    # show neighbourhood but do not take H's into account
    yasara.ShowRes("Element !H Obj "+objnum+ " with distance < 5 from Element !H Res " + residue.number.inpdb + " Mol "+residue.molecule.name + " Obj "+objnum)
    # get a list of residues that make clashes
    vdwlist = yasara.ListRes("Element C Obj "+objnum+ " with distance < 3 from Element C SideChain Res " + residue.number.inpdb + " Mol "+residue.molecule.name + " Obj "+objnum,"RESNUM")
    # only color clashes if we have any. This avoids coloring the mutated residue even if it doesn't make clashes
    if len(vdwlist) > 1:
      yasara.ColorRes("Element C Obj "+objnum+ " with distance < 3 from Element C SideChain Res " + residue.number.inpdb + " Mol "+residue.molecule.name + " Obj "+objnum,"Red")

# show Van der Waals clashes made by new residue stretch
# ======================================================
def ShowVdwclashes_MutateResidueStretch(residues,objnumbers,repairobjnumber):
  # make a list of all active objects to show hbonds
  if repairobjnumber != 0:
    objnumbers = str(residues.residue[0].object.number.inyas) + " " + objnumbers + " " + str(repairobjnumber)
  else:
    objnumbers = str(residues.residue[0].object.number.inyas) + " " + objnumbers
  # make a real list
  objnumberlist = objnumbers.split()
  # make a string with residuenumbers
  resnumbers = ""
  for i in range(residues.residues):
    resnumbers += str(residues.residue[i].number.inpdb) + " "
  # turn it into a list
  resnumberlist = resnumbers.split()
  # loop over every object to show mutation site neighbourhood and VdW clashes
  for objnum in objnumberlist:
    # and loop over every mutated residue to spot clashes per mutation
    for resnum in resnumberlist:
      # show neighbourhood but do not take H's into account
      yasara.ShowRes("Element !H Obj "+objnum+ " with distance < 5 from Element !H Res " + resnum + " Mol "+residues.residue[0].molecule.name + " Obj "+objnum)
      # get a list of residues that make clashes
      clashlist = yasara.ListRes("Element C Obj "+objnum+ " with distance < 3 from Element C SideChain Res " + resnum + " Mol "+residues.residue[0].molecule.name + " Obj "+objnum,"RESNUM")
      # only color clashes if we have any. This avoids coloring the mutated residue even if it doesn't make clashes
      if len(clashlist) > 1:
        yasara.ColorRes("Element C Obj "+objnum+ " with distance < 3 from Element C Sidechain Res " + resnum + " Mol "+residues.residue[0].molecule.name + " Obj "+objnum,"Red")

# Show new or disrupted hydrogen bonds around mutation site
# =========================================================
def ShowHbonds(residue,objnumbers,repairobjnumber,vdwclashes):
  # make a list of all active objects to show hbonds
  if repairobjnumber != 0:
    objnumbers = str(residue.object.number.inyas) + " " + objnumbers + " " + str(repairobjnumber)
  else:
    objnumbers = str(residue.object.number.inyas) + " " + objnumbers
  # make a real list
  objnumberlist = objnumbers.split()
  # list residues that make a contact with wt or mt residue site
  contactreslist = []
  for objnum in objnumberlist:
    contactreslist += yasara.ListRes("Element !H Obj "+objnum+ " with distance < 5 from Element !H Res " + residue.number.inpdb + " Mol "+residue.molecule.name + " Obj "+objnum,"RESNUM MOLNAME")
  # add new hydrogens to all these objects (from here the residue.number.inyas ID is broken) (this has to come after contactreslist to measure distances without hydrogens!)
  yasara.AddHydObj(objnumbers)
  # plot all sidechains that have anything to do with the mutation site (wt and mt) and plot them in both wt and mt to see all differences
  for contactres in contactreslist:
    contactresnumber = contactres.split()[0]
    # detect nameless molecules (these have a space and so [1] is not possible)
    try:
      contactresmolname = contactres.split()[1]
      yasara.ShowRes(contactresnumber+" Mol "+contactresmolname+ " Obj "+objnumbers)
    except:
      yasara.ShowRes(contactresnumber+" Obj "+objnumbers)
    # show hydrogen bonds of contacting residue
    yasara.ShowHBoAtom("SideChain Res "+contactresnumber+ " Obj "+objnumbers,"Yes")
  # hide hydrogens
  yasara.HideAtom("Element H")

# Show new or disrupted hydrogen bonds around mutation site after Mutate Multiple Residues
# ========================================================================================
def ShowHbonds_MutateResidueStretch(residues,objnumbers,repairobjnumber,vdwclashes):
  # make a list of all active objects to show hbonds
  if repairobjnumber != 0:
    objnumbers = str(residues.residue[0].object.number.inyas) + " " + objnumbers + " " + str(repairobjnumber)
  else:
    objnumbers = str(residues.residue[0].object.number.inyas) + " " + objnumbers
  # make a real list
  objnumberlist = objnumbers.split()
  # make a string with residuenumbers
  resnumbers = ""
  for i in range(residues.residues):
    resnumbers += str(residues.residue[i].number.inpdb) + " "
  # list residues that make a contact with wt or mt residue site
  contactreslist = []
  for objnum in objnumberlist:
    contactreslist += yasara.ListRes("Element !H Obj "+objnum+ " with distance < 5 from Element !H Res " + resnumbers + " Mol "+residues.residue[0].molecule.name + " Obj "+objnum,"RESNUM MOLNAME")
  # add new hydrogens to all these objects (from here the residue.number.inyas ID is broken) (this has to come after contactreslist to measure distances without hydrogens!)
  yasara.AddHydObj(objnumbers)
  # plot all sidechains that have anything to do with the mutation site (wt and mt) and plot them in both wt and mt to see all differences
  for contactres in contactreslist:
    contactresnumber = contactres.split()[0]
    # detect nameless molecules (these have a space and so [1] is not possible)
    try:
      contactresmolname = contactres.split()[1]
      yasara.ShowRes(contactresnumber+" Mol "+contactresmolname+ " Obj "+objnumbers)
    except:
      yasara.ShowRes(contactresnumber+" Obj "+objnumbers)
    # show hydrogen bonds of contacting residue
    yasara.ShowHBoAtom("SideChain Res "+contactresnumber+ " Obj "+objnumbers,"Yes")
  # hide hydrogens
  yasara.HideAtom("Element H")

# yasara.plugin.end("Stop it.")

def ColorAccordingToEnergy(allData):
  
  resColors={}
  
  #yasara.ColorRes("all","Gray")
  for energy in allData:
    #Parsing object number from file
    oNumber = re.search('./Object(?P<id>\w+).pdb',energy[0],re.IGNORECASE).groups('id')[0]
    chain=getattr(energy,'chain')
    name=getattr(energy,'three_letter')
    num=getattr(energy,'pdb_seq_num')
    vdwclash=float(getattr(energy,'energy_vdwclash'))

    key=str(name + " " + num + " mol " + chain + " obj " + str(oNumber))
    if not key in list(resColors.keys()):
      resColors[key]="None"
    #else:
	#  yasara.ShowMessage(resColors[key])
    #=====================================#
    if(1>=vdwclash>=0.6):
      if(resColors[key]!="Red"):
        #yasara.ColorRes(key,"Yellow")
        resColors[key]="Yellow"
    if(vdwclash>1.5):
      yasara.ColorRes(key,"Red")
      resColors[key]="Red"
    if(vdwclash<0.1):
      if(resColors[key]!="Red" and resColors[key]!="Yellow"):
        yasara.ColorRes(key,"Green")
        resColors[key]="Green"
  #yasara.ShowMessage(resColors)
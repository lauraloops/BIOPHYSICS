#!/usr/bin/python

# foldxrepair.py

import yasara,string,disk,os,foldxutilities
from python2to3 import *

# residue conversion dictionnary
aa_dict = {'ALA':'A','CYS':'C','ASP':'D','GLU':'E','PHE':'F','GLY':'G','HIS':'H','ILE':'I','LYS':'K','LEU':'L','MET':'M','ASN':'N','PRO':'P','GLN':'Q','ARG':'R','SER':'S','THR':'T','VAL':'V','TRP':'W','TYR':'Y','H1S':'H','H2S':'H'}

# Add stability option to cfile
# =============================
def AddOption(cachedir,option):
  optionsfile = open(os.path.join(cachedir,"config_RP.cfg"),"a")
  optionsfile.write(option+"\n")
  optionsfile.close()

# Add RepairPDB options to cfile
# =============================
def AddRepairOptions(cachedir,fixedreslist,objectnumber):
  AddOption(cachedir,"vdwDesign=2")
  AddOption(cachedir,"out-pdb=true")
  if fixedreslist != []:
    AddOption(cachedir,"fix-residues-file=fixed_residues")
    fixstring=foldxutilities.SetFixResidues(cachedir,fixedreslist,objectnumber)

# Write RepairPDB config file
# ============================
def WriteConfigFile(cachedir,objectnumber,fixedreslist):
  # do normal Repair
  configfile = open(os.path.join(cachedir,"config_RP.cfg"),"w")
  configfile.write("command=RepairPDB\n")
  configfile.close()

# RepairPDB
# =========
def Repair_Object(cachedir,foldxbin,objectnumber,fixedreslist):
  # save the object as PDB for FoldX input
  pdbfilename = "Object"+str(objectnumber)+".pdb"
  yasara.SavePDB(objectnumber,os.path.join(cachedir,pdbfilename))
  # remove the REMARK lines as FoldX crashed on those
  lines = open(os.path.join(cachedir,pdbfilename),"r").readlines()
  output = open(os.path.join(cachedir,pdbfilename),"w")
  for line in lines:
    if line[0:6] != "REMARK":
      output.write(line)
  output.close()
  # make the OPTIONS file
  WriteConfigFile(cachedir,objectnumber,fixedreslist)
  AddRepairOptions(cachedir,fixedreslist,objectnumber)
  AddOption(cachedir,"pdb="+pdbfilename)
  # build the RepairPDB command
  repaircommand = foldxbin + " -f config_RP.cfg"
  # change to temp folder
  os.chdir(cachedir)
  # citation to FoldX
  citation="Please cite Schymkowitz et al., Nucl. Acids Res. 2005; 33:W382-8."
  # show RepairPDB message
  yasara.ShowMessage("Repairing Object "+str(objectnumber)+". This can take quite a few minutes. "+citation)
  # run RepairPDB
  os.system(repaircommand)
  # load repaired structure in YASARA and put objectnumber in a list
  try:
    newobjlist = yasara.LoadPDB(os.path.join(cachedir,"Object"+str(objectnumber)+"_Repair.pdb"))
  except:
    yasara.plugin.end("Cannot read FoldX output. Make sure you have set the correct FoldX file locations in Analyze > FoldX > Configure plugin")
  # show superpose message
  yasara.ShowMessage("Superposing structures...")
  # superpose repaired PDB with original object
  yasara.AlignObj(newobjlist[0],objectnumber,"MUSTANGPP")
  # rename repaired object to RepairPDB
  yasara.NameObj(newobjlist[0],"RepairObj"+str(objectnumber))
  # show LoadPDB message of repaired structure
  yasara.ShowMessage("Repaired PDB loaded in YASARA soup as Object "+str(newobjlist[0])+". You can save this PDB when the plugin procedure ends. Please wait ...")
  yasara.Wait(5,"Seconds")
  # return the objectnumber of the repaired PDB
  return newobjlist[0]

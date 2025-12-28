#!/usr/bin/python

# foldxdihedrals.py

import yasara,string,disk,os
from python2to3 import *

# residue conversion dictionnary
aa_dict = {'ALA':'A','CYS':'C','ASP':'D','GLU':'E','PHE':'F','GLY':'G','HIS':'H','ILE':'I','LYS':'K','LEU':'L','MET':'M','ASN':'N','PRO':'P','GLN':'Q','ARG':'R','SER':'S','THR':'T','VAL':'V','TRP':'W','TYR':'Y','H1S':'H','H2S':'H'}

# Write dihedrals command file
# ============================
def WriteConfigFile(cachedir):
  optionsfile = open(os.path.join(cachedir,"config_DH.cfg"),"w")
  optionsfile.close()
  AddOption(cachedir,"command=Dihedrals")
  AddOption(cachedir,"noHeader=true")
  AddOption(cachedir,"overwriteBatch=false")

# Write dihedrals options file
# ============================
def AddOption(cachedir,option):
  optionsfile = open(os.path.join(cachedir,"config_DH.cfg"),"a")
  optionsfile.write(option+"\n")
  optionsfile.close()

# Dihedrals
# ==============
def Dihedrals_Object(cachedir,foldxbin,objectyas):
  # set pdb filename
  pdbfilename = "Object"+str(objectyas.number.inyas)+".pdb"
  # save the pdb
  yasara.SavePDB(objectyas.number.inyas,os.path.join(cachedir,pdbfilename))
  # change to temp folder
  os.chdir(cachedir)
  # citation to FoldX
  citation="Please cite Schymkowitz et al., Nucl. Acids Res. 2005; 33:W382-8."
  # show BuildModel message
  yasara.ShowMessage("Running FoldX Dihedrals on Object "+str(objectyas.number.inyas)+". "+citation)
  # make command file
  WriteConfigFile(cachedir)
  # run dihedrals for selected object
  AddOption(cachedir,"pdb="+pdbfilename)
  # build the FoldX command
  foldxcommand = foldxbin+ " -f config_DH.cfg"
  # run Dihedrals
  os.system(foldxcommand)
  # stop showing messages on screen        
  yasara.HideMessage()


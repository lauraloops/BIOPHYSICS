#!/usr/bin/python

# foldxstability.py

import yasara,string,disk,os
from python2to3 import *

# residue conversion dictionnary
aa_dict = {'ALA':'A','CYS':'C','ASP':'D','GLU':'E','PHE':'F','GLY':'G','HIS':'H','ILE':'I','LYS':'K','LEU':'L','MET':'M','ASN':'N','PRO':'P','GLN':'Q','ARG':'R','SER':'S','THR':'T','VAL':'V','TRP':'W','TYR':'Y','H1S':'H','H2S':'H'}

# Write stability command file
# ============================
def WriteConfigFile(cachedir):
  optionsfile = open(os.path.join(cachedir,"config_ST.cfg"),"w")
  optionsfile.close()
  AddOption(cachedir,"command=Stability")
  AddOption(cachedir,"out-pdb=false")
  AddOption(cachedir,"noHeader=true")
  AddOption(cachedir,"output-file=tmpST")
  AddOption(cachedir,"overwriteBatch=false")

# Write stability options file
# ============================
def AddOption(cachedir,option):
  optionsfile = open(os.path.join(cachedir,"config_ST.cfg"),"a")
  optionsfile.write(option+"\n")
  optionsfile.close()

# Write stability options file
# ============================
def Stability_LoopReconstruction(cachedir,foldxbin,objectnumber,repair,PDBDictionary):
  # set pdb filename
  for x in sorted(PDBDictionary.keys()):
    
    pdbfilename = PDBDictionary[x]
    # change to temp folder
    os.chdir(cachedir)
    # citation to FoldX
    citation="Please cite Schymkowitz et al., Nucl. Acids Res. 2005; 33:W382-8."
    # show BuildModel message
    yasara.ShowMessage("Running FoldX Stability on Object "+str(x)+". "+citation)
    # make command file
    WriteConfigFile(cachedir)
    # run stability for selected object
    AddOption(cachedir,"pdb="+os.path.basename(pdbfilename))
    # build the FoldX command
    foldxcommand = foldxbin+ " -f config_ST.cfg"
    # run Stability
    os.system(foldxcommand)
    # stop showing messages on screen        
    yasara.HideMessage()

# Stability
# ==============
def Stability_Object(cachedir,foldxbin,objectyas):
  # set pdb filename
  pdbfilename = "Object"+str(objectyas.number.inyas)+".pdb"
  # save the pdb
  yasara.SavePDB(objectyas.number.inyas,os.path.join(cachedir,pdbfilename))
  # change to temp folder
  os.chdir(cachedir)
  # citation to FoldX
  citation="Please cite Schymkowitz et al., Nucl. Acids Res. 2005; 33:W382-8."
  # show BuildModel message
  yasara.ShowMessage("Running FoldX Stability on Object "+str(objectyas.number.inyas)+". "+citation)
  # make command file
  WriteConfigFile(cachedir)
  # run stability for selected object
  AddOption(cachedir,"pdb="+pdbfilename)
  # build the FoldX command
  foldxcommand = foldxbin+ " -f config_ST.cfg"
  # run Stability
  os.system(foldxcommand)
  # stop showing messages on screen        
  yasara.HideMessage()


#!/usr/bin/python

# foldxoptimize.py

import yasara,string,disk,shutil,time,os
from python2to3 import *

# residue conversion dictionnary
aa_dict = {'ALA':'A','CYS':'C','ASP':'D','GLU':'E','PHE':'F','GLY':'G','HIS':'H','ILE':'I','LYS':'K','LEU':'L','MET':'M','ASN':'N','PRO':'P','GLN':'Q','ARG':'R','SER':'S','THR':'T','VAL':'V','TRP':'W','TYR':'Y','H1S':'H','H2S':'H'}

# Write optimize command file
# ============================
def WriteConfigFile(cachedir):
  optionsfile = open(os.path.join(cachedir,"config_OP.cfg"),"w")
  optionsfile.close()
  AddOption(cachedir,"command=Optimize")
  
# Write optimize options file
# ============================
def AddOption(cachedir,option):
  optionsfile = open(os.path.join(cachedir,"config_OP.cfg"),"a")
  optionsfile.write(option+"\n")
  optionsfile.close()

# Optimize
# ==============
def Optimize_Object(cachedir,foldxbin,objectyas):
  # set pdb filename
  pdbfilename = "Object"+str(objectyas.number.inyas)+".pdb"
  # save the pdb
  yasara.SavePDB(objectyas.number.inyas,os.path.join(cachedir,pdbfilename))
  # change to temp folder
  os.chdir(cachedir)
  # citation to FoldX
  citation="Please cite Schymkowitz et al., Nucl. Acids Res. 2005; 33:W382-8."
  # show BuildModel message
  yasara.ShowMessage("Running FoldX Optimize on Object "+str(objectyas.number.inyas)+". "+citation)
  # make command file
  WriteConfigFile(cachedir)
  # run optimize for selected object
  AddOption(cachedir,"pdb="+pdbfilename)
  # build the FoldX command
  foldxcommand = foldxbin+ " -f config_OP.cfg"
  # run Optimize
  os.system(foldxcommand)
  #yasara.DelObj(objectyas.number.inyas)
  shutil.move(os.path.join(cachedir,"Optimized_"+os.path.basename(pdbfilename)),os.path.join(cachedir,pdbfilename))
  objectTmp=yasara.LoadPDB(os.path.join(cachedir,pdbfilename))
  # stop showing messages on screen
  yasara.HideMessage()
  return objectTmp

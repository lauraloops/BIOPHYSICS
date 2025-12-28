#!/usr/bin/python

# foldxpssm.py

import yasara,string,disk,os
from python2to3 import *

# residue conversion dictionnary
aa_dict = {'ALA':'A','CYS':'C','ASP':'D','GLU':'E','PHE':'F','GLY':'G','HIS':'H','ILE':'I','LYS':'K','LEU':'L','MET':'M','ASN':'N','PRO':'P','GLN':'Q','ARG':'R','SER':'S','THR':'T','VAL':'V','TRP':'W','TYR':'Y','H1S':'H','H2S':'H'}

# make pssm command file
# ================================
def WriteConfigFilePM(cachedir):
  configfile = open(os.path.join(cachedir,"config_PM.cfg"),"w")
  configfile.write("command=Pssm\n")
  configfile.close()

def AddOption(cachedir,option):
  configfile = open(os.path.join(cachedir,"config_PM.cfg"),"a")
  configfile.write(option+"\n")
  configfile.close()

# InteractionPSSM
# ===============
def Pssm_Molecules(cachedir,foldxbin,mollist,aminoacids,positions,objectnumber):
  # set pdb filename
  pdbfilename = "Object"+str(objectnumber)+".pdb"
  # save pdb
  yasara.SavePDB(objectnumber,os.path.join(cachedir,pdbfilename))
  # change to temp folder
  os.chdir(cachedir)
  # citation to FoldX
  citation="Please cite Schymkowitz et al., Nucl. Acids Res. 2005; 33:W382-8."
  # show BuildModel message
  yasara.ShowMessage("Running FoldX PSSM. "+citation)
  # run pssm for selected molecule range
  WriteConfigFilePM(cachedir)
  AddOption(cachedir,"pdb="+pdbfilename)
  AddOption(cachedir,"output-file=yasaraPSSM.fxout")
  AddOption(cachedir,"analyseComplexChains="+mollist[0]+","+mollist[1])
  AddOption(cachedir,"positions="+positions)
  AddOption(cachedir,"noHeader=true")
  if(aminoacids!=""):
    AddOption(cachedir,"aminoacids="+aminoacids)
  # build the FoldX command
  foldxcommand = foldxbin+ " -f config_PM.cfg"
  # run AnalyseComplex
  os.system(foldxcommand)
  # stop showing messages on screen        
  yasara.HideMessage()

#!/usr/bin/python

# foldxpositionscan.py

import yasara,string,disk,os,foldxutilities
from python2to3 import *

# residue conversion dictionnary
aa_dict = {'ALA':'A','CYS':'C','ASP':'D','GLU':'E','PHE':'F','GLY':'G','HIS':'H','ILE':'I','LYS':'K','LEU':'L','MET':'M','ASN':'N','PRO':'P','GLN':'Q','ARG':'R','SER':'S','THR':'T','VAL':'V','TRP':'W','TYR':'Y','PTR':'y','TPO':'p','TYS':'z','SEP':'s','HYP':'h','H1S':'H','H2S':'H','H3S':'H'}
##'HYP':'h','MLZ':'k','MLY':'m','M3L':'l','H1S':'o','H2S':'e','H3S':'f'## these are not yet implemented in FoldX

# Write PositionScan OPTIONS file
# =============================
def AddOption(cachedir,option):
  configfile = open(os.path.join(cachedir,"config_PS.cfg"),"a")
  configfile.write(option+"\n")
  configfile.close()

# Write PositionScan COMMAND file
# =============================
def WriteConfigFile(cachedir):
  configfile = open(os.path.join(cachedir,"config_PS.cfg"),"w")
  configfile.write("command=PositionScan\n")
  configfile.close()
  AddOption(cachedir,"output-dir="+str(cachedir))
  AddOption(cachedir,"noHeader=true")

def PositionScanGetDictionary(cachedir,residues,objectyas):
  ResiduesDictionary = {}

  try: scanOut = open(os.path.join(cachedir,"PS_Object"+str(objectyas.number.inyas)+"_scanning_output.txt"),"r").readlines()
  except: yasara.plugin.end("Cannot read FoldX output. Make sure you have set the correct FoldX file locations in Analyze > FoldX > Configure plugin")

  for line in scanOut:
    res=line.split('\t')[0]
    resSeq=aa_dict[res[0:3]]
    resMol=res[3]
    resNum=res[4:-1]
    resMut=res[-1:]
    
    ddg=line.split('\t')[1][:-1]

    #key->objectyas+resSeq+resMol+resNum+resMut
    #key->objectyas+res
    ResiduesDictionary[str(resNum)+resMut] = [resSeq,resMol,resNum,resMut,ddg]
    
  return ResiduesDictionary

# PositionScan
# ==========
def PositionScan_Object(cachedir,foldxbin,residues,objectyas,residueType):
  positions = ""
  for residue in residues:
    residueSeq = residue.name1
    residueNumber = residue.number.inpdb
    residueMolecule = residue.molecule.name
    positions+=residueSeq+residueMolecule+residueNumber+residueType+","
  positions=positions[:-1]

  pdbfilename = "Object"+str(objectyas.number.inyas)+".pdb"
  # save the pdb
  yasara.SavePDB(objectyas.number.inyas,os.path.join(cachedir,pdbfilename))
  # change to temp folder
  os.chdir(cachedir)
  # citation to FoldX
  citation="Please cite Schymkowitz et al., Nucl. Acids Res. 2005; 33:W382-8."
  # show alascan message
  yasara.ShowMessage("Running FoldX PositionScan on Object "+str(objectyas.number.inyas)+". "+citation)
  # make CONFIG file
  WriteConfigFile(cachedir)
  AddOption(cachedir,"pdb="+pdbfilename)
  AddOption(cachedir,"positions="+positions)
  # build the FoldX command
  foldxcommand = foldxbin+ " -f config_PS.cfg"
  # run Stability
  os.system(foldxcommand)
  
  return PositionScanGetDictionary(cachedir,residues,objectyas)
  

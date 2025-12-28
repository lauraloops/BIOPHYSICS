#!/usr/bin/python

# foldxanalysecomplex.py


import yasara,string,disk,os,foldxutilities
from python2to3 import *

# residue conversion dictionnary
aa_dict = {'ALA':'A','CYS':'C','ASP':'D','GLU':'E','PHE':'F','GLY':'G','HIS':'H','ILE':'I','LYS':'K','LEU':'L','MET':'M','ASN':'N','PRO':'P','GLN':'Q','ARG':'R','SER':'S','THR':'T','VAL':'V','TRP':'W','TYR':'Y','H1S':'H','H2S':'H'}

# make analysecomplex command file
# ================================
def WriteConfigFileAC(cachedir):
  configfile = open(os.path.join(cachedir,"config_AC.cfg"),"w")
  configfile.write("command=AnalyseComplex\n")
  configfile.close()

def AddOption(cachedir,option):
  configfile = open(os.path.join(cachedir,"config_AC.cfg"),"a")
  configfile.write(option+"\n")
  configfile.close()

def LaunchAnalyseComplex(cachedir,pdb,targetmol,foldxbin):
  WriteConfigFileAC(cachedir)
  AddOption(cachedir,"pdb="+os.path.basename(pdb))
  AddOption(cachedir,"output-file=yasara")
  #AddOption(cachedir,"analyseComplexChains="+targetmol)
  AddOption(cachedir,"noHeader=true")
  
  foldxcommand = foldxbin+ " -f config_AC.cfg"
  os.system(foldxcommand) # run AnalyseComplex
  #Parse the output to tmpAC
  try: Alllines = open("Interaction_yasara_AC.fxout","r").readlines()
  except: yasara.plugin.end("Cannot read FoldX output. Make sure you have set the correct FoldX file locations in Analyze > FoldX > Configure plugin")
  appendfile = open("tmpAC","a")
  for line in Alllines:
    if( line.find("Group1	Group2") == -1): appendfile.write(line)
  appendfile.close()

# AnalyseComplex
# ==============
def AnalyseComplex_Mutate(cachedir,foldxbin,moleculename,mollist,objectnumber,mutations,runs,repair):
  # remove selected molecule from mollist
  mollist.remove(moleculename)
  # did we do a repair?
  if repair == 1: pdbfilename = "Object"+str(objectnumber)+"_Repair"
  else: pdbfilename = "Object"+str(objectnumber)
  
  os.chdir(cachedir) # change to temp folder
  
  # citation to FoldX
  citation="Please cite Schymkowitz et al., Nucl. Acids Res. 2005; 33:W382-8."
  yasara.ShowMessage("Running FoldX AnalyseComplex. "+citation)
  
  # loop over all WT and MT pdb's and appends to tmpAC
  if runs > 1:
    for i in range(runs):
      newpdbfilename = pdbfilename + "_1_" + str(i) + ".pdb"
      LaunchAnalyseComplex(cachedir,newpdbfilename,moleculename,foldxbin)
      WTnewpdbfilename = "WT_" + pdbfilename + "_1_" + str(i) + ".pdb"
      LaunchAnalyseComplex(cachedir,WTnewpdbfilename,moleculename,foldxbin)
  else:
    newpdbfilename = pdbfilename + "_1.pdb"
    LaunchAnalyseComplex(cachedir,newpdbfilename,moleculename,foldxbin)
    WTnewpdbfilename = "WT_" + pdbfilename + "_1.pdb"
    LaunchAnalyseComplex(cachedir,WTnewpdbfilename,moleculename,foldxbin)
  
  # stop showing messages on screen
  yasara.HideMessage()

# AnalyseComplex
# ==============
def AnalyseComplex_LoopReconstruction(cachedir,foldxbin,objectnumber,repair,PDBDictionary):
  
  os.chdir(cachedir) # change to temp folder
  
  # citation to FoldX
  citation="Please cite Schymkowitz et al., Nucl. Acids Res. 2005; 33:W382-8."
  yasara.ShowMessage("Running FoldX AnalyseComplex. "+citation)

  # loop over all WT and MT pdb's and appends to tmpAC
  moleculename=""
  for x in sorted(PDBDictionary.keys()):
    newpdbfilename = PDBDictionary[x]
    LaunchAnalyseComplex(cachedir,os.path.basename(newpdbfilename),moleculename,foldxbin)
  
  # stop showing messages on screen
  yasara.HideMessage()

# AnalyseComplex
# ==============
def AnalyseComplex_Molecules(cachedir,foldxbin,mollist,objectnumber):
  # set pdb filename
  pdbfilename = "Object"+str(objectnumber)+".pdb"
  # save pdb
  yasara.SavePDB(objectnumber,os.path.join(cachedir,pdbfilename))
  # change to temp folder
  os.chdir(cachedir)
  # citation to FoldX
  citation="Please cite Schymkowitz et al., Nucl. Acids Res. 2005; 33:W382-8."
  # show BuildModel message
  yasara.ShowMessage("Running FoldX AnalyseComplex. "+citation)
  # run analysecomplex for selected molecule range
  WriteConfigFileAC(cachedir)
  AddOption(cachedir,"pdb="+pdbfilename)
  AddOption(cachedir,"output-file=yasara.fxout")
  AddOption(cachedir,"analyseComplexChains="+mollist[0]+","+mollist[1])
  # build the FoldX command
  foldxcommand = foldxbin+ " -f config_AC.cfg"
  # run AnalyseComplex
  os.system(foldxcommand)
  # stop showing messages on screen        
  yasara.HideMessage()

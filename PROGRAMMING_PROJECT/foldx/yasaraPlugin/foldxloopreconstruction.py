#!/usr/bin/python

# foldxloopreconstruction.py

import yasara,string,disk,os,re,glob,shutil
from python2to3 import *

# residue conversion dictionnary
aa_dict = {'ALA':'A','CYS':'C','ASP':'D','GLU':'E','PHE':'F','GLY':'G','HIS':'H','ILE':'I','LYS':'K','LEU':'L','MET':'M','ASN':'N','PRO':'P','GLN':'Q','ARG':'R','SER':'S','THR':'T','VAL':'V','TRP':'W','TYR':'Y','PTR':'y','TPO':'p','TYS':'z','SEP':'s','HYP':'h'}
##'HYP':'h','MLZ':'k','MLY':'m','M3L':'l','H1S':'o','H2S':'e','H3S':'f'## these are not yet implemented in FoldX

# Write LoopReconstruction OPTIONS file
# =============================
def AddOption(cachedir,option):
  configfile = open(os.path.join(cachedir,"config_LR.cfg"),"a")
  configfile.write(option+"\n")
  configfile.close()

def AddOptions(cachedir,moleculename,anchor1,anchor2,loopSequence,dssp1,dssp2,scop,numClasses,pdbfilename,DBConfig,homologySequence,sequenceIdentityPercent,similaritymatrix):
  configfile = open(os.path.join(cachedir,"config_LR.cfg"),"a")
  AddOption(cachedir,"output-dir="+cachedir)
  AddOption(cachedir,"pdb="+pdbfilename)
  AddOption(cachedir,"LR_chain="+moleculename)
  AddOption(cachedir,"LR_pepStartRes="+str(anchor1))
  AddOption(cachedir,"LR_pepEndRes="+str(anchor2))
  AddOption(cachedir,"LR_sequence="+loopSequence)
  AddOption(cachedir,"LR_dsspBegin="+dssp1)
  AddOption(cachedir,"LR_dsspEnd="+dssp2)
  AddOption(cachedir,"LR_scop="+scop)
  #AddOption(cachedir,"LR_nbrClassesToEvaluate="+str(numClasses))
  AddOption(cachedir,"LR_nbrClassesToEvaluate=200")
  ##### Database options
  AddOption(cachedir,"BriX_database="+DBConfig[0].replace('"',""))
  AddOption(cachedir,"BriX_hostname="+DBConfig[1].replace('"',""))
  AddOption(cachedir,"BriX_user="+DBConfig[2].replace('"',""))
  AddOption(cachedir,"BriX_pass="+DBConfig[3].replace('"',""))
  AddOption(cachedir,"BriX_port="+DBConfig[4].replace('"',""))
  AddOption(cachedir,"BriX_PDBDir="+DBConfig[5].replace('"',"")[:-1])	
  AddOption(cachedir,"LoopX_database="+DBConfig[0].replace('"',""))
  AddOption(cachedir,"LoopX_hostname="+DBConfig[1].replace('"',""))
  AddOption(cachedir,"LoopX_user="+DBConfig[2].replace('"',""))
  AddOption(cachedir,"LoopX_pass="+DBConfig[3].replace('"',""))
  AddOption(cachedir,"LoopX_port="+DBConfig[4].replace('"',""))
  AddOption(cachedir,"LoopX_PDBDir="+DBConfig[5].replace('"',"")[:-1])
  ##### HARDCODED,could be a pannel
  AddOption(cachedir,"LR_screenOuput=1")
  AddOption(cachedir,"LR_backboneEntropyFilter = 0")
  AddOption(cachedir,"LR_backboneOmegaFilter = 0")
  AddOption(cachedir,"LR_clashFilter = 0")
  AddOption(cachedir,"LR_softEPFilter = 0")
  AddOption(cachedir,"LR_foldxStabilityFilter = 0")
  AddOption(cachedir,"LR_entropySortFilter = 0")
  homologyFilter="0"
  if(homologySequence != ""):
    homologyFilter="0"
    AddOption(cachedir,"LR_loopSequenceIdentity = " + sequenceIdentityPercent)
  AddOption(cachedir,"LR_sequenceSimilarityFilter = "+homologyFilter)
  AddOption(cachedir,"LR_similarityMatrix = "+ similaritymatrix)
  AddOption(cachedir,"LR_scopFilter = 0")

# Write BuildModel COMMAND file
# =============================
def WriteConfigFile(cachedir):
  configfile = open(os.path.join(cachedir,"config_LR.cfg"),"w")
  configfile.write("command=LoopReconstruction\n")
  configfile.close()

# Align models
# =============================
def AlignLoopsAndGetDictionary(cachedir,objectnumber,maxnumloops):

  PDBDictionary = {}
  objnumbers = ""
  loopNum = 0
  for infile in glob.glob( os.path.join(cachedir, 'Object*_loop_*.pdb') ):
    try: newobjlist = yasara.LoadPDB(infile)
    except: yasara.write("No loops builded")
    # add this object number to the growing object number list for later superposition with original structure
    PDBDictionary[newobjlist[0]] = infile
    objnumbers += " " + str(newobjlist[0])
    loopNum = loopNum + 1
    if(loopNum>=maxnumloops): break
  
  # add object numbers of loaded Mutated FoldX structures to the string (for superposition)
  yasara.ShowMessage("Superposing structures...")
  ### Superpose FoldX models with original object
  yasara.AlignObj(objnumbers,objectnumber,"MUSTANGPP")
  return PDBDictionary

def LoopOut(cachedir,pdbfilename,anchor1,anchor2,moleculename):
  lines = open(os.path.join(cachedir,pdbfilename),"r").readlines()
  output = open(os.path.join(cachedir,pdbfilename),"w")
  for line in lines:
    if line[0:6] != "REMARK":
      if line[0:4] == "ATOM":
        if (int(line[23:26]) > int(anchor1)) & (int(line[23:26]) < int(anchor2)):
          if str(line[21]) == moleculename: continue
      output.write(line)
  output.close()

def OptimizeLoops(PDBDictionary,cachedir,loopxbinpath,objectnumber,maxnumloops):
  for i in list(PDBDictionary.keys()):
    configfile = open(os.path.join(cachedir,"config_OP.cfg"),"w")
    configfile.write("command=Optimize\n")
    configfile.write("moveNeighbours=false\n")
    configfile.write("pdb="+os.path.basename(PDBDictionary[i])+"\n")
    configfile.write("output-dir="+cachedir+"\n")
    configfile.write("pdb-dir="+cachedir+"\n")
    configfile.close()
    
    os.chdir(cachedir)
    foldxcommand = loopxbinpath + " -f config_OP.cfg"
    # citation to FoldX
    citation="Please cite Schymkowitz et al., Nucl. Acids Res. 2005; 33:W382-8."
    yasara.ShowMessage("Running optimization on "+os.path.basename(PDBDictionary[i])+" This can take a few minutes. "+citation)
    # run Optimize
    os.system(foldxcommand)
    shutil.move(os.path.join(cachedir,"Optimized_"+os.path.basename(PDBDictionary[i])),PDBDictionary[i])

  for i in list(PDBDictionary.keys()): yasara.DelObj(i)

  return AlignLoopsAndGetDictionary(cachedir,objectnumber,maxnumloops)

#####################
# LoopReconstruction
# ===================
def LoopReconstruction(cachedir,loopxbinpath,repair,moleculename,anchor1,anchor2,loopSequence,dssp1,dssp2,scop,maxnumloops,objectnumber,DBConfig,homologySequence,sequenceIdentityPercent,similaritymatrix):
  
  if repair == 1: pdbfilename = "Object"+str(objectnumber)+"_Repair.pdb"
  else: pdbfilename = "Object"+str(objectnumber)+".pdb"

  # make the CONFIG file
  WriteConfigFile(cachedir)

  # add OPTIONS to the config file
  if(re.match("[A-Za-z]", moleculename)):
    AddOptions(cachedir,moleculename,anchor1,anchor2,loopSequence,dssp1,dssp2,scop,200,pdbfilename,DBConfig,homologySequence,sequenceIdentityPercent,similaritymatrix)
    yasara.SavePDB(objectnumber,os.path.join(cachedir,pdbfilename))

    #######################
    #HERE CUT THE LOOP OFF#
    #######################
    LoopOut(cachedir,pdbfilename,anchor1,anchor2,moleculename)
  
    os.chdir(cachedir)
    foldxcommand = loopxbinpath + " -f config_LR.cfg"
  
    # citation to FoldX
    citation="Please cite Schymkowitz et al., Nucl. Acids Res. 2005; 33:W382-8."
    yasara.ShowMessage("Running FoldX LoopReconstruction. This can take a few minutes. "+citation)
    # run BuildModel
    os.system(foldxcommand)
    # stop showing messages on screen        
    yasara.HideMessage()
    return AlignLoopsAndGetDictionary(cachedir,objectnumber,maxnumloops)
  else:
    yasara.plugin.end("Rename the molecule names in your object to alphabetic characters")

#!/usr/bin/python

# foldxbuildmodel.py

import yasara,string,disk,os,foldxutilities
from python2to3 import *

# residue conversion dictionnary
aa_dict = {'ALA':'A','CYS':'C','ASP':'D','GLU':'E','PHE':'F','GLY':'G','HIS':'H','ILE':'I','LYS':'K','LEU':'L','MET':'M','ASN':'N','PRO':'P','GLN':'Q','ARG':'R','SER':'S','THR':'T','VAL':'V','TRP':'W','TYR':'Y','PTR':'y','TPO':'p','TYS':'z','SEP':'s','HYP':'h'}
##'HYP':'h','MLZ':'k','MLY':'m','M3L':'l','H1S':'o','H2S':'e','H3S':'f'## these are not yet implemented in FoldX

# Write BuildModel OPTIONS file
# =============================
def AddOption(cachedir,option):
  configfile = open(os.path.join(cachedir,"config_BM.cfg"),"a")
  configfile.write(option+"\n")
  configfile.close()

def AddOptions(cachedir,moveneighbours,runs,temperature,ph,ionicstrength,vdwdesign,outpdb,fixedreslist,objectnumber):
  configfile = open(os.path.join(cachedir,"config_BM.cfg"),"a")
  AddOption(cachedir,"output-dir="+str(cachedir))
  AddOption(cachedir,"temperature="+str(temperature))
  AddOption(cachedir,"ionStrength="+str(ionicstrength))
  AddOption(cachedir,"pH="+str(ph))
  AddOption(cachedir,"moveNeighbours="+str(bool(moveneighbours)).lower())
  AddOption(cachedir,"vdwDesign="+str(vdwdesign))
  AddOption(cachedir,"numberOfRuns="+str(runs))
  AddOption(cachedir,"out-pdb="+outpdb)
  AddOption(cachedir,"noHeader=true")
  if fixedreslist != []:
    AddOption(cachedir,"fix-residues-file=fixed_residues")
    fixresiduesFile=open(os.path.join(cachedir,"fixed_residues"),"w")
    fixstring=foldxutilities.SetFixResidues(cachedir,fixedreslist,objectnumber)
    fixresiduesFile.write(fixstring)
    fixresiduesFile.close()
  configfile.close()

# Write BuildModel COMMAND file
# =============================
def WriteConfigFile(cachedir,fixstring,mutantfilekind,command):
  configfile = open(os.path.join(cachedir,"config_BM.cfg"),"w")
  if command == "BuildModel": configfile.write("command=BuildModel\n")
  if command == "AlaScan": commandfile.write("command=AlaScan\n")
  if mutantfilekind == "individual": configfile.write("mutant-file=individual_list.txt\n")
  if mutantfilekind == "mutant": configfile.write("mutant-file=mutant_file.txt\n")
  # we have residues to fix
  configfile.close()

# Write BuildModel mutant file
# ============================
def WriteMutantFile(cachedir,templatesequence,targetsequence):
  mutantfile = open(os.path.join(cachedir,"mutant_file.txt"),"w")
  mutantfile.write(templatesequence+"\n")
  mutantfile.write(targetsequence+"\n")
  mutantfile.close()

def SetFixResidues(cachedir,fixedreslist,objectnumber):
  fixstring = ""
  fixedresfile = open(os.path.join(cachedir,"fixed_residues"),"w")
  for res in fixedreslist:
    objnumber = yasara.ListObj(res,"OBJNUM")
    if objnumber[0] == int(objectnumber):
      fixres = yasara.ListRes(res,"RESNAME1 MOLNAME RESNUM")
      fixstring += fixres[0].replace(" ","")+","
  if(fixstring!=""): fixedresfile.write(fixstring[:-1]+";")
  fixedresfile.close()
  return fixstring

def RemoveREMARKs(cachedir,pdbfilename):
  # remove the REMARK lines as FoldX crashed on those
  lines = open(os.path.join(cachedir,pdbfilename),"r").readlines()
  output = open(os.path.join(cachedir,pdbfilename),"w")
  for line in lines:
    if line[0:6] != "REMARK":
      output.write(line)
  output.close()

# BuildModel
# ==========
def BuildModel_MutateResidue(cachedir,foldxbin,residue,mutations,objectnumber,repair,moveneighbours,runs,temperature,ph,ionicstrength,vdwdesign,fixedreslist):
  # make the fixed_residues file
  fixstring = SetFixResidues(cachedir,fixedreslist,objectnumber)
  # make the CONFIG file
  WriteConfigFile(cachedir,fixstring,"individual","BuildModel")
  # add OPTIONS to the config file
  AddOptions(cachedir,moveneighbours,runs,temperature,ph,ionicstrength,vdwdesign,"true",fixedreslist,objectnumber)
  # make the individual list of the point mutation(s)
  indivlist = open(os.path.join(cachedir,"individual_list.txt"),"w")
  for i in range(mutations):
    # build FoldX mutation format
    foldxmutation = residue.name1+residue.molecule.name+str(residue.number.inpdb)+aa_dict[yasara.selection[2].listentry[i].upper()]
    indivlist.write(foldxmutation+";\n")
  indivlist.close()
  # grab the object number
  objectnumber = residue.object.number.inyas
  # Save the object (if no repair was done) and build the FoldX command
  if repair == 1:
    pdbfilename = "Object"+str(objectnumber)+"_Repair.pdb"
  else:
    pdbfilename = "Object"+str(objectnumber)+".pdb"
    yasara.SavePDB(objectnumber,os.path.join(cachedir,pdbfilename))
  # remove the REMARK lines as FoldX crashed on those
  lines = open(os.path.join(cachedir,pdbfilename),"r").readlines()
  output = open(os.path.join(cachedir,pdbfilename),"w")
  for line in lines:
    if line[0:6] != "REMARK":
      output.write(line)
  output.close()
  # build the FoldX command
  AddOption(cachedir,"pdb="+pdbfilename)

  foldxcommand = foldxbin+ " -f config_BM.cfg"
  # change to temp folder
  os.chdir(cachedir)
  # citation to FoldX
  citation="Please cite Schymkowitz et al., Nucl. Acids Res. 2005; 33:W382-8."
  # show BuildModel message
  yasara.ShowMessage("Running FoldX BuildModel. This can take a few minutes. "+citation)
  # run BuildModel
  os.system(foldxcommand)
  # stop showing messages on screen        
  yasara.HideMessage()
  # Load the FoldX pdb's into YASARA
  objnumberlist = []
  for j in range(mutations):
    # FoldX filename creation is different between 1 run or multiple runs
    if runs > 1:
      for i in range(runs):
        # grab object number
        try: newobjlist = yasara.LoadPDB(os.path.join(cachedir,pdbfilename[:-4]+"_"+str(j+1)+"_"+str(i)+".pdb"))
        except: yasara.plugin.end("Cannot read FoldX output. Make sure you have set the correct FoldX file locations in Analyze > FoldX > Configure plugin")
        # build FoldX mutation format to rename object
        foldxmutation = residue.name1+residue.molecule.name+str(residue.number.inpdb)+aa_dict[yasara.selection[2].listentry[j].upper()]
        # rename object
        yasara.NameObj(newobjlist[0],foldxmutation+"_"+str(i+1).replace(" ",""))
        # add this object number to the growing object number list for later superposition with original structure
        objnumberlist += newobjlist
    else:
      # grab object number
      try: newobjlist = yasara.LoadPDB(os.path.join(cachedir,pdbfilename[:-4]+"_"+str(j+1)))
      except: yasara.plugin.end("Cannot read FoldX output. Make sure you have set the correct FoldX file locations in Analyze > FoldX > Configure plugin")
      # build FoldX mutation format to rename object
      foldxmutation = residue.name1+residue.molecule.name+str(residue.number.inpdb)+aa_dict[yasara.selection[2].listentry[j].upper()]
      # rename object
      yasara.NameObj(newobjlist[0],foldxmutation.replace(" ",""))
      # add this object number to the growing object number list for later superposition with original structure
      objnumberlist += newobjlist
  # create empty object string
  objnumbers = ""
  # add object numbers of loaded Mutated FoldX structures to the string (for superposition)
  for objnumber in objnumberlist:
    objnumbers += " " + str(objnumber)
  # show superpose message
  yasara.ShowMessage("Superposing structures...")
  # and superpose FoldX models with original object
  yasara.AlignObj(objnumbers,objectnumber,"MUSTANGPP")
  # return objectnumbers
  return objnumbers

# BuildModel
# ==========
def BuildModel_MutateResidueStretch(cachedir,foldxbin,residues,newsequence,objectnumber,repair,moveneighbours,runs,temperature,ph,ionicstrength,vdwdesign,fixedreslist):
  # make the fixed_residues file
  fixstring = SetFixResidues(cachedir,fixedreslist,objectnumber)
  
  # make the individual list of the mutations to be made
  mutationstring = ""
  indivlist = open(os.path.join(cachedir,"individual_list.txt"),"w")
  for i in range(len(newsequence)):
    # build FoldX mutation format
    foldxmutation = residues.residue[i].name1+residues.residue[i].molecule.name+str(residues.residue[i].number.inpdb)+newsequence[i]
    mutationstring += foldxmutation+","
  indivlist.write(mutationstring[:-1]+";")
  indivlist.close()

  # grab the object number
  objectnumber = residues.residue[0].object.number.inyas
  # Save the object (if no repair was done) and build the FoldX command
  if repair == 1:
    pdbfilename = "Object"+str(objectnumber)+"_Repair.pdb"
  else:
    pdbfilename = "Object"+str(objectnumber)+".pdb"
  yasara.SavePDB(objectnumber,os.path.join(cachedir,pdbfilename))
  
  RemoveREMARKs(cachedir,pdbfilename)
  # build the FoldX command
  # make the CONFIG file
  WriteConfigFile(cachedir,fixstring,"individual","BuildModel")
  AddOptions(cachedir,moveneighbours,runs,temperature,ph,ionicstrength,vdwdesign,"true",fixedreslist,objectnumber)
  AddOption(cachedir,"pdb="+pdbfilename)
  foldxcommand = foldxbin+ " -f config_BM.cfg"
  os.chdir(cachedir)
  citation="Please cite Schymkowitz et al., Nucl. Acids Res. 2005; 33:W382-8."
  yasara.ShowMessage("Running FoldX BuildModel. This can take a few minutes. "+citation)
  os.system(foldxcommand)
  yasara.HideMessage()
  
  # Load the FoldX pdb's into YASARA
  objnumberlist = []
  # FoldX filename creation is different between 1 run or multiple runs
  if runs > 1:
    for i in range(runs):
      # grab object number
      pdbOutFile=pdbfilename.split('/')[-1].split('.')[0] +"_1_"+str(i)+".pdb"
      try: newobjlist = yasara.LoadPDB(os.path.join(cachedir,pdbOutFile))
      except: yasara.plugin.end("Cannot read FoldX output. Make sure you have set the correct FoldX file locations in Analyze > FoldX > Configure plugin")
      # label mutation
      foldxmutation = "Model"
      # rename object
      yasara.NameObj(newobjlist[0],foldxmutation+"_"+str(i+1).replace(" ",""))
      # add this object number to the growing object number list for later superposition with original structure
      objnumberlist += newobjlist
  else:
    # grab object number
    pdbOutFile=pdbfilename.split('/')[-1].split('.')[0] +"_1.pdb"
    try: newobjlist = yasara.LoadPDB(os.path.join(cachedir,pdbOutFile))
    except: yasara.plugin.end("Cannot read FoldX output. Make sure you have set the correct FoldX file locations in Analyze > FoldX > Configure plugin")
    # label mutation
    foldxmutation = "Model"
    # rename object
    yasara.NameObj(newobjlist[0],foldxmutation.replace(" ",""))
    # add this object number to the growing object number list for later superposition with original structure
    objnumberlist += newobjlist
  # create empty object string
  objnumbers = ""
  # add object numbers of loaded Mutated FoldX structures to the string (for superposition)
  for objnumber in objnumberlist:
    objnumbers += " " + str(objnumber)
  
  # show superpose message
  yasara.ShowMessage("Superposing structures...")
  # and superpose FoldX models with original object
  if repair == 1:
    objectnumber=int(objectnumber)+1
  if(yasara.request!="MutateDNA"):
    yasara.AlignObj(objnumbers,objectnumber,"MUSTANGPP")
  # return objnumbers
  return objnumbers

# Write BuildModel COMMAND file
# =============================
def WriteASConfigFile(cachedir,pdbfilename,mutantfilekind,command):
  configfile = open(os.path.join(cachedir,"config_AS.cfg"),"w")
  configfile.write("command=AlaScan\n")
  configfile.write("pdb="+pdbfilename+"\n")
  configfile.write("overwriteBatch=false\n")
  configfile.write("output-file=tmp\n")
  configfile.close()

# Alanine scan
# ============
def Alascan_Object(cachedir,foldxbin,objectyas,i):
  # set pdb filename
  pdbfilename = "Object"+str(objectyas.number.inyas)+".pdb"
  # save the pdb
  yasara.SavePDB(objectyas.number.inyas,os.path.join(cachedir,pdbfilename))
  # change to temp folder
  os.chdir(cachedir)
  # citation to FoldX
  citation="Please cite Schymkowitz et al., Nucl. Acids Res. 2005; 33:W382-8."
  # show alascan message
  yasara.ShowMessage("Running FoldX alascan on Object "+str(objectyas.number.inyas)+". "+citation)
  # make CONFIG file
  WriteASConfigFile(cachedir,pdbfilename,"mutant","AlaScan")
  # run alascan for selected object
  AddOptions(cachedir,0,1,298,7,0.05,2,"false","","")
  # build the FoldX command
  foldxcommand = foldxbin+ " -f config_AS.cfg"
  # run Stability
  os.system(foldxcommand)



#!/usr/bin/python

# foldxutilities.py

import yasara,string,disk,os
from python2to3 import *

# residue conversion dictionnary
aa_dict = {'ALA':'A',
           'CYS':'C',
           'ASP':'D',
           'GLU':'E',
           'PHE':'F',
           'GLY':'G',
           'HIS':'H',
           'ILE':'I',
           'LYS':'K',
           'LEU':'L',
           'MET':'M',
           'ASN':'N',
           'PRO':'P',
           'GLN':'Q',
           'ARG':'R',
           'SER':'S',
           'THR':'T',
           'VAL':'V',
           'TRP':'W',
           'TYR':'Y',
           'PTR':'y',
           'TPO':'p',
           'TYS':'z',
           'SEP':'s',
           'HYP':'h'}

def RemoveREMARKs(cachedir,pdbfilename):
  # remove the REMARK lines as FoldX crashed on those
  lines = open(os.path.join(cachedir,pdbfilename),"r").readlines()
  output = open(os.path.join(cachedir,pdbfilename),"w")
  for line in lines:
    if line[0:6] != "REMARK":
      output.write(line)
  output.close()

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

# Print nice table
# ================
def PrintTable(rows):
  #if len(rows) > 1:
  headers = rows[0]._fields
  lens = []
  for i in range(len(rows[0])):
    lens.append(len(max([x[i] for x in rows] + [headers[i]],key=lambda x:len(str(x)))))
  formats = []
  hformats = []
  for i in range(len(rows[0])):
    if isinstance(rows[0][i], int):
      formats.append("%%%dd" % lens[i])
    else:
      formats.append("%%-%ds" % lens[i])
    hformats.append("%%-%ds" % lens[i])
  pattern = " | ".join(formats)
  hpattern = " | ".join(hformats)
  separator = "-+-".join(['-' * n for n in lens])
  print(hpattern % tuple(headers))
  print(separator)
  for line in rows:
    print(pattern % tuple(line))

def PrintTableInFile(cachedir,rows):
  summary = open(os.path.join(cachedir,"FOLDX_SUMMARY.fxout"),"a")
  headers = rows[0]._fields
  lens = []
  for i in range(len(rows[0])):
    lens.append(len(max([x[i] for x in rows] + [headers[i]],key=lambda x:len(str(x)))))
  
  formats = []
  hformats = []
  for i in range(len(rows[0])):
    if isinstance(rows[0][i], int):
      formats.append("%%%dd" % lens[i])
    else:
      formats.append("%%-%ds" % lens[i])
    hformats.append("%%-%ds" % lens[i])
  pattern = " | ".join(formats)
  hpattern = " | ".join(hformats)
  separator = "-+-".join(['-' * n for n in lens])
  summary.write(separator.strip(" ")+'\n')
  summary.write(hpattern % tuple(headers)+'\n')
  summary.write(separator+'\n')
  for line in rows:
    summary.write((pattern % tuple(line)).strip(" ")+"\n")
    summary.write(str(separator).strip(" ")+'\n')

  summary.close()

def PrintTabSeparatedSummary(cachedir,rows):
  summary = open(os.path.join(cachedir,"FOLDX_SUMMARY.tsv"),"a")
  headers = rows[0]._fields
  summary.write(str("\t".join(headers))+"\n")
  for residueRow in rows:
    #for residueField in residueRow:
    summary.write("\t".join(residueRow))
  summary.close()

# Print interface residues between molecules
# ==========================================
def ListInterfaceResidues(cachedir,objectnumber):
  #read the FoldX output files to extract interface residues
  try: lines = open(os.path.join(cachedir,"Interface_Residues_yasara_AC.fxout"),"r").readlines()
  except: yasara.plugin.end("Cannot read FoldX output. Make sure you have set the correct FoldX file locations in Analyze > FoldX > Configure plugin")
  intResidues = lines[9]+lines[10].replace('\t',' ')
  return intResidues

# Create the interaction dictionary
# =================================
def CreateInteractionDictionary(moleculename):
  lines = open("tmpAC","r").readlines()
  EnergyDict = {}
  for line in lines:
    if (line.find("\t"+moleculename+"\t")==-1): continue
    energy=line.split('\t')
    nameObj=energy[0][2:].split('.')[0]

    mol = energy[1]
    targetmol = energy[2]
    intraclashes_mol = energy[3]
    intraclashes_targetmol = energy[4]
    interaction = energy[5]

    key=nameObj+"_"+energy[1]+"_"+energy[2]
    delta_interaction=0
    corrected_interaction=0
    EnergyDict[key] = []
    EnergyDict[key].append([mol,targetmol,intraclashes_mol,intraclashes_targetmol,interaction,delta_interaction,corrected_interaction])

  return EnergyDict

# validate the selected molecules
# ===============================
def ValidateMolecules(mollist):
  if len(mollist) < 2:
    yasara.plugin.end("To calculate interaction energy, the object should contain more than one molecule (= PDB chain). If you are sure the original structure contains more than one molecule, FoldX might be incapable of recognizing certain molecules for energy calculation. See the manual for supported molecules.")
  # if analysecomplex was selected, each molecule in the object must have a unique name
  for mol in mollist:
    molcount = mollist.count(mol)
    if molcount > 1:
      yasara.plugin.end("To calculate FoldX interaction energy, all molecules in the object should have unique names. Molecule "+mol+" occurs "+str(molcount)+" times.")
    if mol == " ":
      yasara.plugin.end("To calculate FoldX interaction energy, all molecules in the object should have unique names. Please rename all nameless molecules with Edit > Rename > Molecule.")

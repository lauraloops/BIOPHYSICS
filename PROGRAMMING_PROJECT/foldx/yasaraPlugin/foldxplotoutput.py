#!/usr/bin/python

# foldxplotoutput.py

import yasara,string,disk,os,math,foldxutilities,foldxyastools
from python2to3 import *
from collections import namedtuple

# residue conversion dictionnary
aa_dict = {'ALA':'A','CYS':'C','ASP':'D','GLU':'E','PHE':'F','GLY':'G','HIS':'H','ILE':'I','LYS':'K','LEU':'L','MET':'M','ASN':'N','PRO':'P','GLN':'Q','ARG':'R','SER':'S','THR':'T','VAL':'V','TRP':'W','TYR':'Y','H1S':'H','H2S':'H'}

# Plot stability energy change output
# =====================================
def PlotStability_Mutate(cachedir,mollist,mutations,runs,repair,objectnumber,mutations3lett):
  ouptufile1="Average_Object"+objectnumber
  ouptufile2="Dif_Object"+objectnumber
  if(repair==1): 
    ouptufile1+="_Repair"
    ouptufile2+="_Repair"
  ouptufile1+=".fxout"
  ouptufile2+=".fxout"
  Av_lines = open(os.path.join(cachedir,ouptufile1),"r").readlines()
  Dif_lines = open(os.path.join(cachedir,ouptufile2),"r").readlines()
  
  Row=namedtuple('Row',['nameObj','total_energy'])
  AvData = []
  for line in Av_lines:
    if (line.find("Pdb\tSD\ttotal energy")!=-1): continue
    energy=line.split('\t')
    nameObj=mutations3lett+":"+energy[0].split('.')[0]
    tenergy = energy[2]
    data=Row(nameObj,tenergy)
    AvData.append(data)

  DifData = []
  for line in Dif_lines:
    if (line.find("Pdb\ttotal energy")!=-1): continue
    energy=line.split('\t')
    nameObj=mutations3lett+":"+energy[0].split('.')[0]
    tenergy = energy[1]
    data=Row(nameObj,tenergy)
    DifData.append(data)
  
  foldxutilities.PrintTable(DifData)
  foldxutilities.PrintTableInFile(cachedir,DifData)
  #foldxutilities.PrintTabSeparatedSummary(cachedir,DifData)
  if runs > 1:
    foldxutilities.PrintTable(AvData)
    foldxutilities.PrintTableInFile(cachedir,AvData)
    #foldxutilities.PrintTabSeparatedSummary(cachedir,AvData)

# Plot interaction energy change output
# =====================================
def PlotAnalyseComplex_LoopReconstruction(cachedir,foldxbin,objectnumber,repair,PDBDictionary,moleculename,numberOfMolecules):
  Row=namedtuple('Row',['Loop','mol','targetmol','intraclashes_mol','intraclashes_targetmol','interaction','delta_interaction','corrected_interaction'])
  Data = []

  Combinations = math.factorial(numberOfMolecules)/(2*math.factorial(numberOfMolecules-2))

  WTypes = {}
  WTenergies=open("tmpAC","r").readlines()[-Combinations:]
  for line in WTenergies:
    WTenergy=line.split('\t')
    Data.append(Row("WT",WTenergy[1],WTenergy[2],WTenergy[3],WTenergy[4],WTenergy[5],"0","0"))
    WTypes[WTenergy[1]+WTenergy[2]] = [WTenergy[3],WTenergy[4],WTenergy[5]]
  
  lines = open("tmpAC","r").readlines()[:-Combinations]
  for line in lines:
    energy=line.split('\t')
    if(line[0].find("_loop_")!=-1): continue

    wildtype=WTypes[energy[1]+energy[2]]
    loop=energy
    
    # calculate interaction energy difference
    delta_interaction = float(loop[5]) - float(wildtype[2])
    # correct interaction energy difference for intraclashes
    delta_intraclashes_mol = float(loop[3]) - float(wildtype[0])
      
    if delta_intraclashes_mol > 0.5: corrected_interaction = delta_interaction + delta_intraclashes_mol
    else: corrected_interaction = delta_interaction

    delta_intraclashes_targetmol = float(loop[4]) - float(wildtype[1])
    if delta_intraclashes_targetmol > 0.5: corrected_interaction = corrected_interaction + delta_intraclashes_targetmol
    Data.append(Row(os.path.basename(energy[0]),energy[1],energy[2],energy[3],energy[4],energy[5],str(delta_interaction),str(corrected_interaction)))
    
  foldxutilities.PrintTable(Data)
  foldxutilities.PrintTabSeparatedSummary(cachedir,Data)
  foldxutilities.PrintTableInFile(cachedir,Data)

# Plot stability energy of loop reconst
# =====================================
def PlotStability_LoopReconstruction(cachedir,foldxbin,objectnumber,repair,PDBDictionary):
  Row=namedtuple('Row',['nameObj','Stability','DDG_Stability'])
  
  WTenergy=open("tmpST_ST.fxout","r").readlines()
  for line in WTenergy:
    energy=line.split('\t')
    if(energy[0].find("_loop_")!=-1): continue
    WTenergy = energy[1]

  Data = []
  Data.append(Row("WT",WTenergy,"0"))
  lines = open("tmpST_ST.fxout","r").readlines()[:-1]
  for line in lines:
    energy=line.split('\t')
    if(energy[0].find("_loop_")==-1): continue
    DG = str(float(energy[1])-float(WTenergy))
    data=Row(os.path.basename(energy[0]),energy[1],DG)
    Data.append(data)

  foldxutilities.PrintTable(Data)
  foldxutilities.PrintTableInFile(cachedir,Data)
  foldxutilities.PrintTabSeparatedSummary(cachedir,Data)

# Plot interaction energy upon mutation
# =====================================
def PlotPSSM_Interaction(cachedir,aminoacids,objectyas):
  aaVector = ['Res']
  for j in range(len(aminoacids)):
    aaVector.append(aminoacids[j])
  
  #yasara.plugin.end(aaVector)
  Row=namedtuple('Row',aaVector)

  try: lines = open(os.path.join(cachedir,"PSSM_Object"+objectyas+".txt"),"r").readlines()[1:]
  except: yasara.plugin.end("Cannot reared FoldX output. Make sure you have set the correct FoldX file locations in Analyze > FoldX > Configure plugin")

  PSSM = []
  for line in lines:
    residueEnergies = line.split("\t")
    rows = []
    residueName = 1
    for s in residueEnergies:
      if(residueName==1):
        residueName=0
        rows.append(s.strip('\n'))
        continue
      rows.append(str(100*round(float(s.strip('\n'))/100,4)))
    data=Row(*rows)
    PSSM.append(data)
  
  return PSSM

# Plot interaction energy upon mutation
# =====================================
def PlotAnalyseComplex_Mutate(cachedir,moleculename,mollist,mutations,runs,repair,mutations3lett,objectnumber):
  
  EnergyDict=foldxutilities.CreateInteractionDictionary(moleculename)
  allData = []
  
  Row=namedtuple('Row',['Mutation','mol','targetmol','intraclashes_mol','intraclashes_targetmol','interaction','delta_interaction','corrected_interaction'])
  for i in sorted(EnergyDict.keys()):
    if (i.find("WT_")!=-1): continue
    else:
      mutant=EnergyDict[i][0]
      wildtype=EnergyDict["WT_"+i][0]
      # calculate interaction energy difference
      delta_interaction = float(mutant[4]) - float(wildtype[4])
      # correct interaction energy difference for intraclashes
      delta_intraclashes_mol = float(mutant[2]) - float(wildtype[2])
      
      if delta_intraclashes_mol > 0.5: corrected_interaction = delta_interaction + delta_intraclashes_mol
      else: corrected_interaction = delta_interaction
      delta_intraclashes_targetmol = float(mutant[3]) - float(wildtype[3])
      if delta_intraclashes_targetmol > 0.5: corrected_interaction = corrected_interaction + delta_intraclashes_targetmol
      EnergyDict[i][0][5]=delta_interaction
      EnergyDict[i][0][6]=corrected_interaction
      
      data=Row(mutations3lett+": "+str(i),mutant[0],mutant[1],mutant[2],mutant[3],mutant[4],mutant[5],mutant[6])
      dataWT=Row(mutations3lett+": WT_"+str(i),wildtype[0],wildtype[1],wildtype[2],wildtype[3],wildtype[4],wildtype[5],wildtype[6])
      allData.append(data)
      allData.append(dataWT)
  foldxutilities.PrintTable(allData)
  foldxutilities.PrintTableInFile(cachedir,allData)
  if runs==1: return
  
  allData = []
  averageEnergy = {}
  if(repair==1): keyRoot="Object"+objectnumber+"_Repair_1_"
  else: keyRoot="Object"+objectnumber+"_1_"
  for mol in sorted(mollist):
    sortedmolecules=sorted([str(moleculename),str(mol)])
    avKeyName=keyRoot+sortedmolecules[0]+"_"+sortedmolecules[1]
    averageEnergy[avKeyName] = []
    averageEnergy[avKeyName].append(sortedmolecules[0])
    averageEnergy[avKeyName].append(sortedmolecules[1])
    for k in range(2,7):
      for o in range(runs):
        keyName=keyRoot+str(o)+"_"+sortedmolecules[0]+"_"+sortedmolecules[1]
        if o==0: averageEnergy[avKeyName].append(EnergyDict[keyName][0][k])
        else:
          alreadyIn=float(averageEnergy[avKeyName][k])
          averageEnergy[avKeyName][k] = float(EnergyDict[keyName][0][k]) + alreadyIn
      
      averageEnergy[avKeyName][k]=float(averageEnergy[avKeyName][k])/runs
    
    avE=averageEnergy[avKeyName]
    data=Row(mutations3lett+": "+keyRoot+sortedmolecules[0]+"_"+sortedmolecules[1],str(avE[0]),str(avE[1]),str(avE[2]),str(avE[3]),str(avE[4]),str(avE[5]),str(avE[6]))
    allData.append(data)
    
  foldxutilities.PrintTable(allData)
  foldxutilities.PrintTableInFile(cachedir,allData)

# Print interaction energy between molecules
# =========================================
def PlotDihedrals_Object(cachedir,foldxbin,obj):
  try: lines = open(os.path.join(cachedir,"DH_Object"+obj.number.inyas+".fxout"),"r").readlines()
  except: yasara.plugin.end("Cannot read FoldX output. Make sure you have set the correct FoldX file locations in Analyze > FoldX > Configure plugin")

  Row=namedtuple('Row',['Object','threeLetter','seq','omega','phi','psi'])
  dihedrals = []
  metEnergy=""
  for line in lines:
    if line[0:12].find("./Object"+str(obj.number.inyas)) != -1:
      Energy = line.split("\t")
      data=Row(*Energy[0:6])
      dihedrals.append(data)

  return dihedrals

def PlotDihedrals_Residue(cachedir,foldxbin,residue):
  try: lines = open(os.path.join(cachedir,"DH_Object"+residue.object.number.inyas+".fxout"),"r").readlines()
  except: yasara.plugin.end("Cannot read FoldX output. Make sure you have set the correct FoldX file locations in Analyze > FoldX > Configure plugin")

  Row=namedtuple('Row',['Object','threeLetter','seq','omega','phi','psi'])
  dihedrals = []
  metEnergy=""
  for line in lines:
    if line[0:12].find("./Object"+str(residue.object.number.inyas)) != -1:
      residueInFile = line.split("\t")
      if(residue.name3 == residueInFile[1] and residue.molecule.name == residueInFile[2] and str(residue.number.inpdb) == residueInFile[3]):
        dihedrals.append(Row(*residueInFile[0:6]))

  return dihedrals

# Print interaction energy between molecules
# =========================================
def PlotMetalBinding_Object(cachedir,foldxbin,obj,metal):
  try: lines = open(os.path.join(cachedir,"MB_Object"+obj.number.inyas+".fxout"),"r").readlines()
  except: yasara.plugin.end("Cannot read FoldX output. Make sure you have set the correct FoldX file locations in Analyze > FoldX > Configure plugin")
  
  Row=namedtuple('Row',['Object','metal','energy'])
  metalsPredicted = []
  metEnergy=""
  for line in lines:
    if line[0:15].find("./"+metal+"_MB_Object"+str(obj.number.inyas)) != -1:
      Energy = line.split("\t")
      data=Row(*Energy[0:3])
      metalsPredicted.append(data)
      break
  
  return metalsPredicted

# Print interaction energy between molecules
# =========================================
def PlotAnalyseComplex_Molecules(cachedir,mollist,objectnumber):
  #read the FoldX output files to extract average energies per mutation
  intEnergy = ""
  try: lines = open(os.path.join(cachedir,"Interaction_yasara_AC.fxout"),"r").readlines()
  except: yasara.plugin.end("Cannot read FoldX output. Make sure you have set the correct FoldX file locations in Analyze > FoldX > Configure plugin")
  for line in lines:
    if line[0:8].find("./Object") != -1:
      # put the average energies in a list
      Energy = line.split("\t")
      intEnergy += "Interaction energy between molecule(s) "+mollist[0]+" and "+mollist[1]+" in object "+str(objectnumber)+" = " + str(round(float(Energy[5]),2)) + " (kcal/mol)\n"
  return intEnergy

# Print stability energy between molecules
# =========================================
def PlotStability_Object(cachedir,foldxbin,objectyas,i,mutation3letts):
  #read the FoldX output files to extract average energies per mutation
  stability = ""
  try: line = open(os.path.join(cachedir,"tmpST_ST.fxout"),"r").readlines()[-1]
  except: yasara.plugin.end("Cannot read FoldX output. Make sure you have set the correct FoldX file locations in Analyze > FoldX > Configure plugin")
  if line[0:8].find("./Object") != -1:
    # put the average energies in a list
    stability = line.split('\t')
  # and put the output in a formatted string table
  # start with table for all individual run energies
  stability = "Stability of"+mutation3letts+": "+str(objectyas.name)+" = " + str(round(float(stability[1]),2)) + " (kcal/mol)"
  return stability

# Print stability energy between molecules
# =========================================
def PlotSequenceDetail_Object(cachedir,foldxbin,objectyas,i,mutation3letts):
  #read the FoldX output files
  try: lines = open(os.path.join(cachedir,"SD_Object"+str(objectyas.number.inyas)+".fxout"),"r").readlines()
  except: yasara.plugin.end("Cannot read FoldX output. Make sure you have set the correct FoldX file locations in Analyze > FoldX > Configure plugin")
  
  fields = ['obj',
            'three_letter','chain','pdb_seq_num','omega',
            'phi','psi','sec_struct','electro','entrop_mc',
            'entrop_sc','sideHbond','backHbond','energy_VdW',
            'energy_SolvP','energy_SolvH','energy_vdwclash',
            'sloop_entropy','mloop_entropy','energy_torsion',
            'backbone_vdwclash','total','cis_bond','water',
            'disulfide','energy_dipole','energy_kon','partcov',
            'sc_partcov','energy_sc_SolvH','energy_sc_VdW',
            'energy_sc_SolvP','entr_complex','electroSc',
            'energyKonSc','energyIonisation','ab_index']

  Row=namedtuple('Row',fields)
  
  allData = []
  for line in lines:
    if line[0:8].find("./Object") != -1:
      residue=line.split('\t')
      if(residue[1]=="GAP"): continue
      allData.append(Row(*residue))

  foldxutilities.PrintTable(allData)
  #foldxutilities.PrintTableInFile(cachedir,allData)
  foldxutilities.PrintTabSeparatedSummary(cachedir,allData)
  return allData

def PlotSequenceDetail_Residue(cachedir,foldxbin,residue):
  #read the FoldX output files to extract average energies per mutation
  
  fields = ['obj',
            'three_letter','chain','pdb_seq_num','omega',
            'phi','psi','sec_struct','electro','entrop_mc',
            'entrop_sc','sideHbond','backHbond','energy_VdW',
            'energy_SolvP','energy_SolvH','energy_vdwclash',
            'sloop_entropy','mloop_entropy','energy_torsion',
            'backbone_vdwclash','total','cis_bond','water',
            'disulfide','energy_dipole','energy_kon','partcov',
            'sc_partcov','energy_sc_SolvH','energy_sc_VdW',
            'energy_sc_SolvP','entr_complex','electroSc',
            'energyKonSc','energyIonisation','ab_index']

  Row=namedtuple('Row',fields)

  allData = []
  try: lines = open(os.path.join(cachedir,"SD_Object"+str(residue.object.number.inyas)+".fxout"),"r").readlines()
  except: yasara.plugin.end("Cannot read FoldX output. Make sure you have set the correct FoldX file locations in Analyze > FoldX > Configure plugin")
  
  for line in lines:
    if line[0:8].find("./Object") != -1:
      residueInFile=line.split('\t')
      if(residue.name3 == residueInFile[1] and residue.molecule.name == residueInFile[2] and str(residue.number.inpdb) == residueInFile[3]):
        allData.append(Row(*residueInFile))

  foldxutilities.PrintTable(allData)
  #foldxutilities.PrintTableInFile(cachedir,allData)
  foldxutilities.PrintTabSeparatedSummary(cachedir,allData)
  
  return allData

# Print positionScan results per object
# =====================================

def PlotPositionScan(cachedir,foldxbin,ResDictionary,objectyas):
  allData = []
  Row=namedtuple('Row',['WTSeq','Mol','resNumber','MTseq','ddg'])
  
  for key in sorted(ResDictionary.keys()):
    data=Row(ResDictionary[key][0],ResDictionary[key][1],ResDictionary[key][2],ResDictionary[key][3],ResDictionary[key][4])
    allData.append(data)
  foldxutilities.PrintTable(allData)
  #foldxutilities.PrintTableInFile(cachedir,allData)
  foldxutilities.PrintTabSeparatedSummary(cachedir,allData)
  return 

# Print alascan results per object
# ================================
def PlotAlascan_Object(cachedir,foldxbin,objectyas,i):
  #read the FoldX output files to extract alascan results
  allData = []
  try: lines = open(os.path.join(cachedir,"tmp_AS.fxout"),"r").readlines()
  except: yasara.plugin.end("Cannot read FoldX output. Make sure you have set the correct FoldX file locations in Analyze > FoldX > Configure plugin")
  
  Row=namedtuple('Row',['Residue','Alascan_ddG'])
  for line in lines:
    if line.find("ALA energy change") != -1:
      # grab details
      resname = line.split()[0]
      resnum = line.split()[1]
      energy = line.split()[-1]
      energy = energy.replace("\n","")
      energy = energy.replace("\r","")
      data=Row(str(resname)+str(resnum),energy)
      allData.append(data)
  # print nice table with all results from all objects
  yasara.write("ALANINE SCAN RESULTS FOR OBJECT "+str(objectyas.number.inyas))
  foldxutilities.PrintTable(allData)
  foldxutilities.PrintTableInFile(cachedir,allData)
  foldxutilities.PrintTabSeparatedSummary(cachedir,allData)
  return 


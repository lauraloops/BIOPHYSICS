# YASARA PLUGIN
# TOPIC:       Molecular modeling
# TITLE:       FoldX forcefield plugin
# AUTHOR:      Javier Delgado, Joost Van Durme. CRG-VIB
# LICENSE:     GPL (www.gnu.org)
# DESCRIPTION: This plugin provides access to various tools from the FoldX molecular modeling and design program in YASARA.
#              FoldX can be obtained from http://foldx.crg.es.
###########################################################################################################################
# This is a YASARA plugin to be placed in the /plg subdirectory #
# Go to www.yasara.org/plugins for documentation and downloads  #
#################################################################
# Web address of the FoldX site to download the program: http://foldx.switchlab.org or http://foldx.crg.es
# Please cite the FoldX article when publishing results from this plugin:
# Schymkowitz J, Borg J, Stricher F, Nys R, Rousseau F, Serrano L. 
# The FoldX web server: an online force field.
# Nucleic Acids Res. 2005 Jul 1;33(Web Server issue):W382-8.
# Version 2.1.12
#
"""
MainMenu: Analyze
  PullDownMenu after Energy: FoldX
    PopUpMenu: Fix residues
      ResidueSelectionMenu: Select residue(s) not to move in FoldX analysis
      Request: FixResidues

    PopUpMenu: Free residues
      ResidueSelectionMenu: Select residue(s) to unfix for FoldX analysis
      Request: FreeResidues

    PopUpMenu: Optimize object sidechains
      ObjectSelectionMenu: Select object(s) to optimize
      Request: Optimize

    PopUpMenu: Repair object
      ObjectSelectionMenu: Select the object(s) for FoldX RepairPDB
      Request: RepairObject

    PopUpMenu: Stability of object
      ObjectSelectionMenu: Select the object(s) for stability calculation
      Request: Stability

    PopUpMenu: Alanine scan of object
      ObjectSelectionMenu: Select the object(s) for FoldX alanine scan
      Request: Alascan

    PopUpMenu: Interaction energy of molecules
      MoleculeSelectionMenu: Select first molecule range
      MoleculeSelectionMenu: Select second molecule range
      Request: AnalyseComplex

    PopUpMenu: Interaction PSSM of molecules
      MoleculeSelectionMenu: Select first molecule range
      MoleculeSelectionMenu: Select second molecule range
      ResidueSelectionMenu: Select a residue for FoldX analysis	
      TextInputMenu1: New amino acids
        Text: Type the new amino acids in the box in correct order
        Text: New residue list (case sensitive!):
      Request: InteractionPSSM

    PopUpMenu: Mutate residue
      ResidueSelectionMenu: Select a residue for FoldX analysis
      CheckBoxMenu3: Select FoldX routines
        Text: Select FoldX actions
        Text: FoldX _R_epairPDB
        Text: Calculate _s_tability change (Checked)
        Text: Calculate _i_nteraction energy change
      ListMenu: Select new amino acid residue(s)
        MultipleSelections: Yes
        Text: Select new amino acid(s)
        TextFile: foldxaminoacids.txt
      CheckBoxMenu4: Set FoldX options (1)
        Text: Enable or disable FoldX options and visualization (default is shown):
        Text: Move neighbours (Checked)
        Text: Zoom to mutation site
        Text: Show disrupted and new hydrogen bonds
        Text: Show VdW clashes in WT and mutant
      NumberInputMenu5: Set FoldX options (2)
        Text: Set values for FoldX options (default is shown):
        Number: Number of runs,1,1,5
        Number: Temperature (K),298,0,500
        Number: pH,7,1,14
        Number: Ionic strength (x100),5,0,10
        Number: VdW design,2,0,2
      Request: MutateResidue

    PopUpMenu: Mutate multiple residues
      ResidueSelectionMenu: Select residues for FoldX analysis
      CheckBoxMenu3: Select FoldX routines
        Text: Select FoldX actions
        Text: FoldX _R_epairPDB
        Text: Calculate _s_tability change (Checked)
        Text: Calculate _i_nteraction energy change
      TextInputMenu1: New amino acids
        Text: Type the new amino acids in the box in correct order
        Text: New residues (case sensitive!):
      CheckBoxMenu4: Set FoldX options (1)
        Text: Enable or disable FoldX options and visualization (default is shown):
        Text: Move neighbours (Checked)
        Text: Zoom to mutation site
        Text: Show disrupted and new hydrogen bonds
        Text: Show VdW clashes in WT and mutant
      NumberInputMenu5: Set FoldX options (2)
        Text: Set values for FoldX options (default is shown):
        Number: Number of runs,1,1,5
        Number: Temperature (K),298,0,500
        Number: pH,7,1,14
        Number: Ionic strength (x100),5,0,10
        Number: VdW design,2,0,2
      Request: MutateResidueStretch

    PopUpMenu: Mutate DNA
      ResidueSelectionMenu: Select bases for FoldX analysis
      CheckBoxMenu3: Select FoldX routines
        Text: Select FoldX actions
        Text: FoldX _R_epairPDB
        Text: Calculate _s_tability change (Checked)
        Text: Calculate _i_nteraction energy change
      TextInputMenu1: Base Code
        Text: Type the new base code in the box in correct order
        Text: New residues (case sensitive!):
      CheckBoxMenu4: Set FoldX options (1)
        Text: Enable or disable FoldX options and visualization (default is shown):
        Text: Move neighbours (Checked)
        Text: Zoom to mutation site
        Text: Show disrupted and new hydrogen bonds
        Text: Show VdW clashes in WT and mutant
      NumberInputMenu5: Set FoldX options (2)
        Text: Set values for FoldX options (default is shown):
        Number: Number of runs,1,1,5
        Number: Temperature (K),298,0,500
        Number: pH,7,1,14
        Number: Ionic strength (x100),5,0,10
        Number: VdW design,2,0,2
      Request: MutateDNA
    
    PopUpMenu: Sequence detail of object
      ObjectSelectionMenu: Select the object(s) for sequence detail calculation
      Request: SequenceDetail
    
    PopUpMenu: Metal Binding prediction on object
      ObjectSelectionMenu: Select metal for FoldX Metal Binding prediction
      ListMenu: Select metal type for Metal Binding prediction
        MultipleSelections: No
        Text: Select metal for Metal Binding prediction
        TextFile: foldxmetals.txt
      Request: MetalBinding
    
    PopUpMenu: Position scan
      ResidueSelectionMenu: Select residues for FoldX position scan analysis
      ListMenu: Select aa types for position scan, look manual for descriptions
        MultipleSelections: No
        Text: Select aa types for position scan
        TextFile: foldxaatypes.txt
      Request: PositionScan
    
    PopUpMenu: Loop reconstruction
      ResidueSelectionMenu: Select anchor residue 1
      ResidueSelectionMenu: Select anchor residue 2
      CheckBoxMenu4: Select FoldX routines
        Text: Select FoldX actions
        Text: FoldX _R_epairPDB
        Text: Calculate _s_tability change (Checked)
        Text: Calculate _i_nteraction energy change
        Text: Reconstruct _s_ide chains
      TextInputMenu2: Loop sequence for sequence identity filtering
        Text: Loop sequence identity
        Text: Insert sequence:
        Text: Percentage of identity between loops,0
      ListMenu: Select dssp for anchor residue 1 (optional)
        MultipleSelections: No
        Text: Select dssp for anchor residue 1
        TextFile: foldxdssptypes.txt
      ListMenu: Select dssp for anchor2 (optional)
        MultipleSelections: No
        Text: Select dssp for anchor residue 2
        TextFile: foldxdssptypes.txt
      ListMenu: Select scop family (optional)
        MultipleSelections: Yes
        Text: Select scop families
        TextFile: foldxscopfamilies.txt
      NumberInputMenu1: Max number of loops to return
        Text: Max number of loops to load on scene:
        Number: Number of loops to load on scene,1,1,1000
      Request: LoopReconstruction
    
    PopUpMenu: Color residues by FoldX energy
      ObjectSelectionMenu: Select the object(s) for foldx energy color
      Request: ColorResidues

    PopUpMenu: Backbone Dihedrals on object
      ObjectSelectionMenu: Select the object(s) for dihedrals calculation
      Request: Dihedrals

    PopUpMenu: Save last calculation
      FileSelectionMenu: Select a folder to save files, it must exists
        MultipleSelections: No
        Filename: *
      Request: SaveCalc

    PopUpMenu: Configure plugin
      FileSelectionMenu: Select your FoldX executable file
        MultipleSelections: No
        Filename: *
      FileSelectionMenu: Select your FoldX rotabase.txt file
        MultipleSelections: No
        Filename: *
      Request: FileLocations

    PopUpMenu: Configure LoopX plugin
      FileSelectionMenu: Select your FoldX with LoopX executable file
        MultipleSelections: No
        Filename: *
      TextInputMenu1: Set up LoopX database
        Text: LoopX Database name
        Text: LoopX Database name:, LoopX
      TextInputMenu4: Set up LoopX
        Text: LoopX parameters:
        Text: LoopX host name, localhost
        Text: LoopX database user, loopx
        Text: Select your LoopX database password, loopx
        Text: Select your LoopX port, 3306
      FileSelectionMenu: Select your loop database pdb directory
        MultipleSelections: No
        Filename: *
      Request: LoopXFileLocations

AtomContextMenu after Swap: FoldX
  PopUpMenu: Fix residue
    Request: FixResidues_FromAtom

  PopUpMenu: Free residue
    Request: FreeResidues_FromAtom

  PopUpMenu: Repair object
    Request: RepairObject_FromAtom

  PopUpMenu: Stability of object
    Request: Stability_FromAtom

  PopUpMenu: Optimize object
    Request: Optimize_FromAtom

  PopUpMenu: Mutate residue
    CheckBoxMenu3: Select FoldX routines
      Text: Select FoldX actions
      Text: FoldX _R_epairPDB
      Text: Calculate _s_tability change (Checked)
      Text: Calculate _i_nteraction energy change
    ListMenu: Select new amino acid residue(s)
      MultipleSelections: Yes
      Text: Select new amino acid(s)
      TextFile: foldxaminoacids.txt
    CheckBoxMenu4: Set FoldX options (1)
      Text: Enable or disable FoldX options and visualization (default is shown):
      Text: Move neighbours (Checked)
      Text: Zoom to mutation site
      Text: Show disrupted and new hydrogen bonds
      Text: Show VdW clashes in WT and mutant
    NumberInputMenu5: Set FoldX options (2)
      Text: Set values for FoldX options (default is shown):
      Number: Number of runs,1,1,5
      Number: Temperature (K),298,0,500
      Number: pH,7,1,14
      Number: Ionic strength (x100),5,0,10
      Number: VdW design,2,0,2
    Request: MutateResidue_FromAtom

  PopUpMenu: Sequence Detail
    Request: SequenceDetail_FromAtom

  PopUpMenu: Metal Binding prediction
    ListMenu: Select metal type for Metal Binding prediction
      MultipleSelections: No
      Text: Select metal for Metal Binding prediction
      TextFile: foldxmetals.txt
    Request: MetalBinding_FromAtom

  PopUpMenu: Color residue by FoldX energy
    Request: ColorResidues_FromAtom

  PopUpMenu: Optimize object sidechains
    Request: Optimize_FromAtom

  PopUpMenu: Backbone Dihedrals on object
    Request: Dihedrals_FromAtom

HUDMoleculeContextMenu after Split: FoldX
  PopUpMenu: Repair object
    Request: RepairObject_FromMolecule

  PopUpMenu: Stability object
    Request: Stability_FromMolecule

  PopUpMenu: Optimize object
    Request: Optimize_FromMolecule

  PopUpMenu: Metal Binding
    ListMenu: Select metal type for Metal Binding prediction
      MultipleSelections: No
      Text: Select metal for Metal Binding prediction
      TextFile: foldxmetals.txt
    Request: MetalBinding_FromMolecule

  PopUpMenu: Color residue by FoldX energy
    Request: ColorResidues_FromMolecule

HUDObjectContextMenu after Split: FoldX
  PopUpMenu: Repair object
    Request: RepairObject_FromObject

  PopUpMenu: Stability object
    Request: Stability_FromObject

  PopUpMenu: Optimize object
    Request: Optimize_FromObject

  PopUpMenu: Metal Binding
    ListMenu: Select metal type for Metal Binding prediction
      MultipleSelections: No
      Text: Select metal for Metal Binding prediction
      TextFile: foldxmetals.txt
    Request: MetalBinding_FromObject

  PopUpMenu: Backbone Dihedrals on object
    Request: Dihedrals_FromObject

  PopUpMenu: Color residue by FoldX energy
    Request: ColorResidues_FromObject

"""
#!/usr/bin/python
import yasara,string,disk,os,platform
import foldxrepair,foldxoptimize,foldxbuildmodel,foldxanalysecomplex,foldxplotoutput
import foldxsequencedetail,foldxloopreconstruction,foldxyastools,foldxpositionscan
import foldxdihedrals,foldxmetalbinding,foldxstability,foldxpssm,foldxutilities

from python2to3 import *

# amino acid residue list
aminoacidDict = {}
aas = ['A','C','D','E','F','G','H','I','K','L','M','N','P','Q','R','S','T','V','W','Y','y','p','s','z']

aas2 = ['ALA','CYS','ASP','GLU','PHE','GLY','HIS','ILE','LYS','LEU','MET','ASN','PRO','GLN','ARG','SER','THR','VAL','TRP','TYR','PTR','TPO','SEP','TYS']

nucs=['a','c','g','t','u']

for index in range( len(aas) ):
  aminoacidDict[ aas2[index ] ] = aas[index]

# set non-calculation requests so we don't clean the cache folder after these. 
noncalc = ['FixResidues','FixResidues_FromResidue','FixResidues_FromAtom','FreeResidues','FreeResidues_FromResidue','FreeResidues_FromAtom','SaveCalc','FileLocations']
# redundant files that can be removed when saving FoldX output to new folder
redundantfiles = ['rotabase.txt','PdbList_']

# Main program
# ============

# check if operating system is OK, otherwise disable the plugin
if (yasara.request=="CheckIfDisabled"): yasara.plugin.exitcode=0
  # ASSIGN A 1 TO yasara.plugin.exitcode IF THIS PLUGIN CANNOT WORK AND SHOULD
  # BE DISABLED (DATA MISSING, WRONG OPERATING SYSTEM ETC.)
  # check operating system

# *************************************************************************************************************************************************
# * set file locations of binary and rotabase, this section needs to be the first one, otherwise even this section will raise an error at the end *
# *************************************************************************************************************************************************
elif (yasara.request=="FileLocations"):
  foldxbinpath = yasara.selection[0].filename[0]
  rotabasepath = yasara.selection[1].filename[0]
  # stop if wrong rotabase file was selected (just by filename)
  if platform.system() == "Windows":
    if not os.path.split(rotabasepath)[-1] == 'rotabase.txt':
      yasara.plugin.end("Invalid rotabase.txt file. Download FoldX from foldx.crg.es")
  lines = open("foldx.cnf","r").readlines()
  output = open("foldx.cnf","w")
  for line in lines:
    if line[0:9] == "FOLDX_BIN": output.write("FOLDX_BIN = "+foldxbinpath+"\n")
    elif line[0:14] == "FOLDX_ROTABASE": output.write("FOLDX_ROTABASE = "+rotabasepath+"\n")
    else: output.write(line)
  output.close()
elif (yasara.request=="LoopXFileLocations"):
  loopxbinpath = yasara.selection[0].filename[0]
  databasename = yasara.selection[1].text[0]
  hostname = yasara.selection[2].text[0]
  username = yasara.selection[2].text[1]
  password = yasara.selection[2].text[2]
  yasara.write(password)
  port = yasara.selection[2].text[3]
  loopxpdbdir = yasara.selection[3].filename[0]
  # stop if wrong rotabase file was selected (just by filename)
  
  lines = open("foldx.cnf","r").readlines()
  output = open("foldx.cnf","w")
  for line in lines:
    if line[0:9] == "LOOPX_BIN": output.write("LOOPX_BIN = "+loopxbinpath+"\n")
    elif line[0:8] == "LOOPX_DB": output.write("LOOPX_DB = "+databasename+"\n")
    elif line[0:10] == "LOOPX_HOST": output.write("LOOPX_HOST = "+hostname+"\n")
    elif line[0:10] == "LOOPX_USER": output.write("LOOPX_USER = "+username+"\n")
    elif line[0:10] == "LOOPX_PASS": output.write("LOOPX_PASS = "+password+"\n")
    elif line[0:10] == "LOOPX_PORT": output.write("LOOPX_PORT = "+port+"\n")
    elif line[0:12] == "LOOPX_PDBDIR": output.write("LOOPX_PDBDIR = "+loopxpdbdir+"\n")
    else: output.write(line)
  output.close()
# continue if OK, do some basic initializing
else:
  # check for path of FoldX binary
  if not os.path.exists(yasara.plugin.config["FOLDX_BIN"]):
    yasara.plugin.end("Cannot locate the FoldX executable file. Please click Analyze > FoldX > Configure plugin to select the FoldX executable file. Download FoldX from foldx.crg.es")
  # check for path of rotabase
  if platform.system() == "Windows":
    if not os.path.exists(yasara.plugin.config["FOLDX_ROTABASE"]) or not os.path.split(yasara.plugin.config["FOLDX_ROTABASE"])[-1] == 'rotabase.txt':
      yasara.plugin.end("Cannot locate the FoldX rotabase.txt file. Please click Analyze > FoldX > Configure plugin to select the FoldX rotabase.txt file. Download FoldX from foldx.crg.es")
  # hide Console for speed gain
  yasara.Console("Hidden")
  # set the cache dir
  cachedir = os.path.join(os.getcwd(),yasara.plugin.config["FOLDX_TMP"])
  # Delete temp folder only if a calculation is chosen. This saves the 'Save last calculation' plugin command from deletion after eg. FixResidues ...
  if yasara.request not in noncalc: disk.rmdir(cachedir)
  # make the cachedir if not already present
  if (not os.path.exists(cachedir)): disk.makedirs(cachedir,yasara.permissions)
  # copy the rotabase file to the cachedir
  if platform.system() == "Windows":
    disk.copy(yasara.plugin.config["FOLDX_ROTABASE"],os.path.join(cachedir,"rotabase.txt"))
  # get the foldx executable absolute path and put double quotes around it to handle paths with spaces
  foldxbin='"'+yasara.plugin.config["FOLDX_BIN"]+'"'
  database='"'+yasara.plugin.config["LOOPX_DB"]+'"'
  host='"'+yasara.plugin.config["LOOPX_HOST"]+'"'
  user='"'+yasara.plugin.config["LOOPX_USER"]+'"'
  #password='"'+yasara.plugin.config["LOOPX_PASS"]+'"'
  password='"'+yasara.plugin.config["LOOPX_PASS"]+'"'
  port='"'+yasara.plugin.config["LOOPX_PORT"]+'"'
  pdbdir='"'+yasara.plugin.config["LOOPX_PDBDIR"]+'"'
  # put fixed residues in a list, this yasara command returns the unique atom number string (same as residue.number.inyas)
  fixedreslist = yasara.ListRes("Property=10000")

# ************************************************************************************
# * save the last calculation in a specified folder (copy entire cache folder there) *
# ************************************************************************************
if (yasara.request=="SaveCalc"):
  # grab the selected copy folder
  savetarget = yasara.selection[0].filename[0]
  # set the cache dir from the config file
  cachedir = os.path.join(os.getcwd(),yasara.plugin.config["FOLDX_TMP"])
  # remove redundant cache files. don't copy these to save folder
  for filename in redundantfiles:
    disk.remove(os.path.join(cachedir,filename))
  # check if the path is real
  try:
    # here we assume a folder and not a prefix was selected
    os.listdir(savetarget)
    targetfolder = savetarget
    targetprefix = ""
  except OSError:
    # here we assume a file prefix was given. check if underlying folder is real
    os.listdir(os.path.split(savetarget)[0])
    targetfolder = os.path.split(savetarget)[0]
    targetprefix = os.path.split(savetarget)[1]
  except:
    yasara.plugin.end("No proper folder was selected. You probably specified an unexisting folder.")
  # which files to copy?
  # all files
  cachelist = os.listdir(cachedir)
  # see if we don't overwrite any files. if so, do not continue
  copylist = os.listdir(targetfolder)
  for filename in cachelist:
    if targetprefix != "" and targetprefix != "*" and targetprefix+"_"+filename in copylist:
      yasara.plugin.end("Identical filenames found in source folder and selected save folder. Please empty the save folder or select a new save folder or choose another filename prefix to avoid overwriting.")
    elif (targetprefix == "" or targetprefix == "*") and filename in copylist:
      yasara.plugin.end("Identical filenames found in source folder and selected save folder. Please empty the save folder or select a new save folder or choose another filename prefix to avoid overwriting.")    
  # raise error when cache folder is empty
  if cachelist == []:
    yasara.plugin.end("No recent calculation found.")
  # copy the entire cache dir to the copy folder
  for filename in cachelist:
    disk.copy(os.path.join(cachedir,filename),os.path.join(targetfolder,filename))
    yasara.write("Downloaded "+filename+" from "+cachedir+" to "+targetfolder)

# ***************************************************************************
# * Fix residues, give them a property value of 10000 and color them yellow *
# ***************************************************************************
if (yasara.request=="FixResidues"):
  
  # put all residues to be fixed in a string to pass to YASARA. This increases speed when fixing loads of residues instead of looping over all residues.
  fixstring = ""
  for i in range(yasara.selection[0].residues):
    fixstring += str(yasara.selection[0].residue[i].number.inyas) + " "
  yasara.PropRes(fixstring,10000)
  yasara.ColorRes(fixstring,"Yellow")
elif (yasara.request=="FixResidues_FromAtom"):
  yasara.PropRes(yasara.selection[0].atom[0].residue.number.inyas,10000)
  yasara.ColorRes(yasara.selection[0].atom[0].residue.number.inyas,"Yellow")
elif (yasara.request=="FixResidues_FromResidue"):
  yasara.PropRes(yasara.selection[0].residue[0].number.inyas,10000)
  yasara.ColorRes(yasara.selection[0].residue[0].number.inyas,"Yellow")

# *********************************************************************************
# * Free residues, give them a property value of 0 and color them back to element *
# *********************************************************************************
if (yasara.request=="FreeResidues"):
	# put all residues to be freed in a string to pass to YASARA. This increases speed when freeing loads of residues instead of looping over all residues.
  freestring = ""
  for i in range(yasara.selection[0].residues):
    freestring += str(yasara.selection[0].residue[i].number.inyas) + " "
  yasara.PropRes(freestring,0)
  yasara.ColorRes(freestring,"Element")
elif (yasara.request=="FreeResidues_FromAtom"):
  yasara.PropRes(yasara.selection[0].atom[0].residue.number.inyas,0)
  yasara.ColorRes(yasara.selection[0].atom[0].residue.number.inyas,"Element")
elif (yasara.request=="FreeResidues_FromResidue"):
  yasara.PropRes(yasara.selection[0].residue[0].number.inyas,0)
  yasara.ColorRes(yasara.selection[0].residue[0].number.inyas,"Element")

# **********
# * Repair *
# **********
if (yasara.request=="RepairObject") \
or (yasara.request=="RepairObject_FromAtom") \
or (yasara.request=="RepairObject_FromResidue") \
or (yasara.request=="RepairObject_FromMolecule") \
or (yasara.request=="RepairObject_FromObject"):
  # count selected objects
  if (yasara.request=="RepairObject"):
    objcount = yasara.selection[0].objects
    # loop over the selected objects
    for i in range(objcount):
      objectnumber = yasara.selection[0].object[i].number.inyas
      foldxrepair.Repair_Object(cachedir,foldxbin,objectnumber,fixedreslist)
  # grab the object number, depends on selected menu
  elif (yasara.request=="RepairObject_FromAtom"): objectnumber = yasara.selection[0].atom[0].object.number.inyas
  elif (yasara.request=="RepairObject_FromResidue"): objectnumber = yasara.selection[0].residue[0].object.number.inyas
  elif (yasara.request=="RepairObject_FromMolecule"): objectnumber = yasara.selection[0].molecule[0].object.number.inyas
  elif (yasara.request=="RepairObject_FromObject"): objectnumber = yasara.selection[0].object[0].number.inyas
  # Repair the PDB if there's certainly only 1 object to repair
  if yasara.request!="RepairObject": foldxrepair.Repair_Object(cachedir,foldxbin,objectnumber,fixedreslist)
  # stop showing messages on screen        
  yasara.HideMessage()

# ******************
# * Mutate residue *
# ******************
elif (yasara.request=="MutateResidue") \
or  (yasara.request=="MutateResidue_FromAtom") \
or  (yasara.request=="MutateResidue_FromResidue"):
  
  # grab residue details, depends on selected menu
  if (yasara.request=="MutateResidue"): residue = yasara.selection[0].residue[0]
  elif (yasara.request=="MutateResidue_FromAtom"): residue = yasara.selection[0].atom[0].residue
  elif (yasara.request=="MutateResidue_FromResidue"): 
    residue = yasara.selection[0].residue[0]

  # count how many residues to mutate
  #rescount = 
  resWT1Lett=residue.name1
  # grab the new residues and remove spaces if present
  newsequence = yasara.selection[2].listentries
  # grab routine FoldX options
  if yasara.selection[1].checkbox[0] == 1: repair = 1
  else: repair = 0
  if yasara.selection[1].checkbox[1] == 1: stability = 1
  else: stability = 0
  if yasara.selection[1].checkbox[2] == 1: analysecomplex = 1
  else: analysecomplex = 0
  
  # see how many residues were selected
  if yasara.selection[0].residues > 1: yasara.plugin.end("The FoldX plugin only analyses 1 residue per run. Please try again and select one residue from the leftmost list.")
  
  # grab the object number
  objectnumber = residue.object.number.inyas
  
  # if analysecomplex was selected validate molecules for FoldX
  mollist = yasara.ListMol("Obj "+objectnumber,"MOLNAME")
  
  if analysecomplex == 1: foldxutilities.ValidateMolecules(mollist)

  # run RepairPDB if needed
  if repair == 1: repairobjnumber = foldxrepair.Repair_Object(cachedir,foldxbin,objectnumber,fixedreslist)
  else: repairobjnumber = 0

  # run BuildModel to make the mutation(s)
  # grab mutation(s)
  mutations = yasara.selection[2].listentries
  # grab true or false OPTIONS
  if yasara.selection[3].checkbox[0] == 1: moveneighbours = 1
  else: moveneighbours = 0
  if yasara.selection[3].checkbox[1] == 1: zoom = 1
  else: zoom = 0
  if yasara.selection[3].checkbox[2] == 1: hbonds = 1
  else: hbonds = 0
  if yasara.selection[3].checkbox[3] == 1: vdwclashes = 1
  else: vdwclashes = 0
  
  # grab numerical OPTIONS
  runs = yasara.selection[4].number[0]
  temperature = yasara.selection[4].number[1]
  ph = yasara.selection[4].number[2]
  ionicstrength = str(round(yasara.selection[4].number[3]/100.0,2))
  vdwdesign = yasara.selection[4].number[4]
  
  mutation3letts = []
  mutations3lett = ""
  for entry in yasara.selection[2].listentry:
    mutations3lett+=resWT1Lett+residue.molecule.name+residue.number.inobj+aminoacidDict[entry.upper()]
  # run BuildModel
  objnumbers = foldxbuildmodel.BuildModel_MutateResidue(cachedir,foldxbin,residue,mutations,objectnumber,repair,moveneighbours,runs,temperature,ph,ionicstrength,vdwdesign,fixedreslist)

  # make a new mollist of a newly loaded foldx structure, because could be that FoldX deleted some molecules it didn't recognize
  foldxmollist = yasara.ListMol("Obj "+objnumbers.split()[0],"MOLNAME")
  # do AnalyseComplex if it was selected
  if analysecomplex == 1:
    # see if we still have more than 1 molecule after FoldX might discard some
    foldxutilities.ValidateMolecules(foldxmollist)
    moleculename = residue.molecule.name
    foldxanalysecomplex.AnalyseComplex_Mutate(cachedir,foldxbin,moleculename,foldxmollist,objectnumber,mutations,runs,repair)

  # show vdwclashes if needed (this has to come before hbond visual because otherwise we calculate H-H distance and not C-C e.g. and we delete H's anyway)
  if vdwclashes == 1: foldxyastools.ShowVdwclashes(residue,objnumbers,repairobjnumber)
  # show hbonds if needed
  if hbonds == 1: foldxyastools.ShowHbonds(residue,objnumbers,repairobjnumber,vdwclashes)
  # zoom if needed
  if zoom == 1: foldxyastools.Zoom(residue,objnumbers,repairobjnumber,hbonds,vdwclashes)
  # plot the Stability energy output if selected
  if stability == 1: foldxplotoutput.PlotStability_Mutate(cachedir,mollist,mutations,runs,repair,objectnumber,mutations3lett)

  # plot the AnalyseComplex energy output if selected
  if analysecomplex == 1: foldxplotoutput.PlotAnalyseComplex_Mutate(cachedir,moleculename,foldxmollist,mutations,runs,repair,mutations3lett,objectnumber)
  # hide all messages
  yasara.HideMessage()

  # open the console
  yasara.Console("Open")

# ****************************
# * Mutate multiple residues *
# ****************************
elif (yasara.request=="MutateResidueStretch") or\
(yasara.request=="MutateDNA"):
  # count how many residues to mutate
  rescount = yasara.selection[0].residues
  # grab the new residues and remove spaces if present
  newsequence = yasara.selection[2].text[0].replace(" ","")
  # check the validity of the new sequence
  for i in range(len(newsequence)):
    if(yasara.request=="MutateDNA"):
      if newsequence[i] not in nucs:
        yasara.plugin.end("The character "+newsequence[i]+" is not a standard FoldX base code. Code is case sensitive, please check the documentation for more info.")
    else:
      if newsequence[i] not in aas:
        yasara.plugin.end("The character "+newsequence[i]+" is not a standard FoldX amino acid residue code. Code is case sensitive, please check the documentation for more info.")
  # check if stretch length and new sequence length is the same
  if len(newsequence) != rescount:
    yasara.plugin.end("The number of selected residues and the length of the new sequence do not match. These should be of equal length.")
  # check if all selected residues are in same object
  for i in range(rescount):
    objnum = yasara.selection[0].residue[i].object.number.inall
    # don't go until the very last, we can't count past the end
    if i < rescount-1:
     # check if residues are in same molecule
     if yasara.selection[0].residue[i+1].object.number.inall != objnum:
       yasara.plugin.end("Not all selected residues belong to the same object. Please select residues from the same object")
  
  # grab routine FoldX options
  if yasara.selection[1].checkbox[0] == 1: repair = 1
  else: repair = 0
  if yasara.selection[1].checkbox[1] == 1: stability = 1
  else: stability = 0
  if yasara.selection[1].checkbox[2] == 1: analysecomplex = 1
  else: analysecomplex = 0
  # grab objectnumber
  objectnumber = yasara.selection[0].residue[0].object.number.inyas
  # if analysecomplex was selected we need to have more than 1 molecule
  mollist = yasara.ListMol("Obj "+objectnumber,"MOLNAME")
  if analysecomplex == 1: foldxutilities.ValidateMolecules(mollist)

  # run RepairPDB if needed
  if repair == 1: repairobjnumber = foldxrepair.Repair_Object(cachedir,foldxbin,objectnumber,fixedreslist)
  else: repairobjnumber = 0

  # run BuildModel to make the mutation(s)
  # grab mutation(s)
  mutations = yasara.selection[2].listentries

  # grab true or false OPTIONS
  if yasara.selection[3].checkbox[0] == 1: moveneighbours = 1
  else: moveneighbours = 0
  if yasara.selection[3].checkbox[1] == 1: zoom = 1
  else: zoom = 0
  if yasara.selection[3].checkbox[2] == 1: hbonds = 1
  else: hbonds = 0
  if yasara.selection[3].checkbox[3] == 1: vdwclashes = 1
  else: vdwclashes = 0

  # grab numerical OPTIONS
  runs = yasara.selection[4].number[0]
  temperature = yasara.selection[4].number[1]
  ph = yasara.selection[4].number[2]
  ionicstrength = str(round(yasara.selection[4].number[3]/100.0,2))
  vdwdesign = yasara.selection[4].number[4]

  # make a readable variable with the residue descriptors of the selected stretch
  residues = yasara.selection[0]

  # make a readable variable with the residue descriptors of the selected stretch
  mutations3lett=""
  newsequence = yasara.selection[2].text[0]

  for j in range(yasara.selection[0].residues):
    res=yasara.selection[0].residue[j]
    yasara.write(newsequence[j])
    yasara.write(res.number.inobj)
    yasara.write(newsequence[j])
    mutations3lett+=res.name1+res.molecule.name+res.number.inobj+newsequence[j]+","
  mutations3lett=mutations3lett[:-1]
  # run buildmodel and pass new sequence
  objnumbers = foldxbuildmodel.BuildModel_MutateResidueStretch(cachedir,foldxbin,residues,newsequence,objectnumber,repair,moveneighbours,runs,temperature,ph,ionicstrength,vdwdesign,fixedreslist)

  # make a new mollist of a newly loaded foldx structure, because could be that FoldX deleted some molecules it didn't recognize
  foldxmollist = yasara.ListMol("Obj "+objnumbers.split()[0],"MOLNAME")

  # do AnalyseComplex if it was selected, check how many different molecules we have and calculate interaction energy of all combinations
  if analysecomplex == 1:
    moleculename = residues.residue[0].molecule.name
    foldxanalysecomplex.AnalyseComplex_Mutate(cachedir,foldxbin,moleculename,foldxmollist,objectnumber,mutations,runs,repair)

  # show vdwclashes if needed (this has to come before hbond visual because otherwise we calculate H-H distance and not C-C e.g.
  if vdwclashes == 1: foldxyastools.ShowVdwclashes_MutateResidueStretch(residues,objnumbers,repairobjnumber)
  # show hbonds if needed
  if hbonds == 1: foldxyastools.ShowHbonds_MutateResidueStretch(residues,objnumbers,repairobjnumber,vdwclashes)
  # zoom if needed
  if zoom == 1: foldxyastools.Zoom_MutateResidueStretch(residues,objnumbers,repairobjnumber,hbonds,vdwclashes)
 
  # plot the Stability energy output if selected
  if stability == 1: foldxplotoutput.PlotStability_Mutate(cachedir,mollist,mutations,runs,repair,objectnumber,mutations3lett)

  # plot the AnalyseComplex energy output if selected
  if analysecomplex == 1: foldxplotoutput.PlotAnalyseComplex_Mutate(cachedir,moleculename,foldxmollist,mutations,runs,repair,mutations3lett,objectnumber)

  # hide all messages
  yasara.HideMessage()   
  
  # open the console
  yasara.Console("Open")

# ***********************
# * Loop Reconstruction *
# ***********************
elif (yasara.request=="LoopReconstruction"):

  #If the molecules of both residues is the same we pick the molecule name
  if(yasara.selection[0].residue[0].molecule.name!=yasara.selection[1].residue[0].molecule.name):
    yasara.plugin.end("Anchor residues must be at the same molecule")
  moleculename=yasara.selection[0].residue[0].molecule.name
  
  #Database parameters came from the foldx.cnf
  DBConfig=([database,host,user,password,str(port),pdbdir])
  
  objectnumber=yasara.selection[1].residue[0].object.number.inyas
  foldxmollist = yasara.ListMol("Obj "+objectnumber,"MOLNAME")
  numberOfMolecules = len(foldxmollist)
  
  #Anchor residues
  anchor1 = yasara.selection[0].residue[0].number.inpdb
  anchor2 = yasara.selection[1].residue[0].number.inpdb
  if(anchor1>anchor2):
    anchor2=yasara.selection[0].residue[0].number.inpdb
    anchor1=yasara.selection[1].residue[0].number.inpdb
  
  homologySequence = ""
  homologySequence=yasara.selection[3].text[0].upper()
  sequenceIdentityPercent=yasara.selection[3].text[1]

  #Loop sequence must include anchor residues
  loopSequence=yasara.selection[0].residue[0].name1
  loopSequence+=yasara.selection[3].text[0].upper()
  loopSequence+=yasara.selection[1].residue[0].name1
  loopSequence=loopSequence.upper()
  
  #DSSP info of the anchor points
  dssp1=dssp2=""
  for entry in range(yasara.selection[4].listentries):
    dssp1=yasara.selection[4].listentry[entry].upper()
  for entry in range(yasara.selection[5].listentries):
    dssp2=yasara.selection[5].listentry[entry].upper()
  #Scop family is optional
  scop=""
  for entry in range(yasara.selection[6].listentries):
    scop+=yasara.selection[6].listentry[entry]+","
  scop=scop[:-1]
  #Classes to evaluate
  numClasses=yasara.selection[7].number[0]

  # grab routine FoldX options
  if yasara.selection[2].checkbox[0] == 1: repair = 1
  else: repair = 0
  if yasara.selection[2].checkbox[1] == 1: stability = 1
  else: stability = 0
  if yasara.selection[2].checkbox[2] == 1: analysecomplex = 1
  else: analysecomplex = 0
  optimizeSideChains = yasara.selection[2].checkbox[3]
  # if analysecomplex was selected we need to have more than 1 molecule
  mollist = yasara.ListMol("Obj "+objectnumber,"MOLNAME")

  # run RepairPDB if needed
  if repair == 1: repairobjnumber = foldxrepair.Repair_Object(cachedir,foldxbin,objectnumber,fixedreslist)
  
  PDBDictionary=foldxloopreconstruction.LoopReconstruction(cachedir,
                                                            foldxbin,
                                                            repair,
                                                            moleculename,
                                                            anchor1,
                                                            anchor2,
                                                            loopSequence,
                                                            dssp1,
                                                            dssp2,
                                                            scop,
                                                            numClasses,
                                                            objectnumber,
                                                            DBConfig,
                                                            homologySequence,
                                                            sequenceIdentityPercent)
  if(optimizeSideChains == 1):
    PDBDictionary=foldxloopreconstruction.OptimizeLoops(PDBDictionary,cachedir,foldxbin,objectnumber,numClasses)
  # did we do a repair?
  if repair == 1: pdbfilename = "Object"+str(objectnumber)+"_Repair.pdb"
  else: pdbfilename = "Object"+str(objectnumber)+".pdb"
  os.chdir(cachedir) # change to temp folder
  PDBDictionary[objectnumber] = pdbfilename

  if analysecomplex == 1:
    foldxanalysecomplex.AnalyseComplex_LoopReconstruction(cachedir,foldxbin,objectnumber,repair,PDBDictionary)
  if stability == 1:
    foldxstability.Stability_LoopReconstruction(cachedir,foldxbin,objectnumber,repair,PDBDictionary)
  if analysecomplex == 1:
    foldxplotoutput.PlotAnalyseComplex_LoopReconstruction(cachedir,foldxbin,objectnumber,repair,PDBDictionary,moleculename,numberOfMolecules)
  if stability == 1:
    foldxplotoutput.PlotStability_LoopReconstruction(cachedir,foldxbin,objectnumber,repair,PDBDictionary)

  # open the console
  yasara.Console("Open")

# ***********************
# *   Position Scan     *
# ***********************
elif (yasara.request=="PositionScan"):

  residueList=yasara.selection[0].residue
  objectyas=yasara.selection[0].residue[0].object
  residueType=yasara.selection[1].listentry[0]
  objectnumber=objectyas.number.inyas

  ResDictionary=foldxpositionscan.PositionScan_Object(cachedir,foldxbin,residueList,objectyas,residueType)
  foldxplotoutput.PlotPositionScan(cachedir,foldxbin,ResDictionary,objectyas)

  # open the console
  yasara.Console("Open")

# ******************
# * AnalyseComplex *
# ******************
elif (yasara.request=="AnalyseComplex"):
  # count molecules
  molcount = yasara.selection[0].molecules + yasara.selection[1].molecules
  # open a list to store all molecule descriptors
  allmols = []
  # make the list of molname selections, put all molecules from each selection in 1 string for FoldX
  mollist = []
  firstrange = ""
  for i in range(yasara.selection[0].molecules):
    firstrange += yasara.selection[0].molecule[i].name
    allmols.append(yasara.selection[0].molecule[i])
  secondrange = ""
  for i in range(yasara.selection[1].molecules):
    secondrange += yasara.selection[1].molecule[i].name
    allmols.append(yasara.selection[1].molecule[i])
  mollist.append(firstrange)
  mollist.append(secondrange)
  # molecules should belong to same object
  objectnumber = allmols[0].object.number.inyas
  for molecule in allmols:
    nextobjectnum = molecule.object.number.inyas
    if nextobjectnum != objectnumber:
      yasara.plugin.end("To calculate interaction energy, the molecules must be in the same object. Please try again and select molecules from one object.")
  # same molecule cannot occur in both selections
  for mol in firstrange:
    if mol in secondrange:
      yasara.plugin.end("First and second selection cannot contain the same molecule "+mol)
  # no nameless molecules
  if firstrange.find(" ") != -1 or secondrange.find(" ") != -1:
    yasara.plugin.end("To calculate FoldX interaction energy, all molecules in the object should have unique names. Please rename all nameless molecules with Edit > Rename > Molecule.")
  # run analysecomplex
  foldxanalysecomplex.AnalyseComplex_Molecules(cachedir,foldxbin,mollist,objectnumber)
  # make plot string of energy
  intEnergy = foldxplotoutput.PlotAnalyseComplex_Molecules(cachedir,mollist,objectnumber)
  # make plot string of interface residues
  intResidues = foldxutilities.ListInterfaceResidues(cachedir,objectnumber)
  # plot interface residues
  yasara.write(intResidues)
  # plot interaction energy
  yasara.write(intEnergy)
  
  # open the console
  yasara.Console("Open")

# ******************
# *      Pssm      *
# ******************
elif (yasara.request=="InteractionPSSM"):
  # count molecules
  molcount = yasara.selection[0].molecules + yasara.selection[1].molecules
  # open a list to store all molecule descriptors
  allmols = []
  # make the list of molname selections, put all molecules from each selection in 1 string for FoldX
  mollist = []
  firstrange = ""
  for i in range(yasara.selection[0].molecules):
    firstrange += yasara.selection[0].molecule[i].name
    allmols.append(yasara.selection[0].molecule[i])
  secondrange = ""
  for i in range(yasara.selection[1].molecules):
    secondrange += yasara.selection[1].molecule[i].name
    allmols.append(yasara.selection[1].molecule[i])
  mollist.append(firstrange)
  mollist.append(secondrange)
  
  # make a readable variable with the residue descriptors of the selected stretch
  positions=""
  aminoacids = yasara.selection[3].text[0]

  for j in range(yasara.selection[2].residues):
    res=yasara.selection[2].residue[j]
    yasara.write(res.number.inobj)
    yasara.write(aminoacids[j])
    positions+=res.name1+res.molecule.name+res.number.inpdb+aminoacids[j]+","
  positions=positions[:-1]
  
  # molecules should belong to same object
  objectnumber = allmols[0].object.number.inyas
  for molecule in allmols:
    nextobjectnum = molecule.object.number.inyas
    if nextobjectnum != objectnumber:
      yasara.plugin.end("To calculate interaction energy, the molecules must be in the same object. Please try again and select molecules from one object.")
  # same molecule cannot occur in both selections
  for mol in firstrange:
    if mol in secondrange:
      yasara.plugin.end("First and second selection cannot contain the same molecule "+mol)
  # no nameless molecules
  if firstrange.find(" ") != -1 or secondrange.find(" ") != -1:
    yasara.plugin.end("To calculate FoldX interaction energy, all molecules in the object should have unique names. Please rename all nameless molecules with Edit > Rename > Molecule.")
  # run analysecomplex
  foldxpssm.Pssm_Molecules(cachedir,foldxbin,mollist,aminoacids,positions,objectnumber)

  # make plot string of energy
  intEnergy = foldxplotoutput.PlotPSSM_Interaction(cachedir,aminoacids,objectnumber)
  foldxutilities.PrintTable(intEnergy)

  # open the console
  yasara.Console("Open")

# *************
# * Stability *
# *************
elif (yasara.request=="Stability") \
or   (yasara.request=="Stability_FromAtom") \
or   (yasara.request=="Stability_FromResidue") \
or   (yasara.request=="Stability_FromMolecule") \
or   (yasara.request=="Stability_FromObject"):
  # make a readable variable with the residue descriptors of the selected stretch

  # create an empty list
  allstabilities = []
  # grab the number of objects
  if yasara.request=="Stability":
    objcount = yasara.selection[0].objects
    # loop over the selected objects
    for i in range(objcount):
      foldxstability.Stability_Object(cachedir,foldxbin,yasara.selection[0].object[i])
      stability = foldxplotoutput.PlotStability_Object(cachedir,foldxbin,yasara.selection[0].object[i],i,"")
      #yasara.plugin.end(stability)
      allstabilities.append(stability)
  elif (yasara.request=="Stability_FromAtom"):
    foldxstability.Stability_Object(cachedir,foldxbin,yasara.selection[0].atom[0].object)
    stability = foldxplotoutput.PlotStability_Object(cachedir,foldxbin,yasara.selection[0].atom[0].object,1,"")
    allstabilities.append(stability)
  elif (yasara.request=="Stability_FromResidue"):
    foldxstability.Stability_Object(cachedir,foldxbin,yasara.selection[0].residue[0].object)
    stability = foldxplotoutput.PlotStability_Object(cachedir,foldxbin,yasara.selection[0].residue[0].object,1,"")
    allstabilities.append(stability)
  elif (yasara.request=="Stability_FromMolecule"):
    foldxstability.Stability_Object(cachedir,foldxbin,yasara.selection[0].molecule[0].object)
    stability = foldxplotoutput.PlotStability_Object(cachedir,foldxbin,yasara.selection[0].molecule[0].object,1,"")
    allstabilities.append(stability)
  elif (yasara.request=="Stability_FromObject"):
    foldxstability.Stability_Object(cachedir,foldxbin,yasara.selection[0].object[0])
    stability = foldxplotoutput.PlotStability_Object(cachedir,foldxbin,yasara.selection[0].object[0],1,"")
    allstabilities.append(stability)

  yasara.HideMessage()
  # print all stabilities of all objects
  for stability in allstabilities:
    print(stability)
    
  # open the console
  yasara.Console("Open")    

# ****************
# * MetalBinding *
# ****************
elif (yasara.request=="MetalBinding") \
or   (yasara.request=="MetalBinding_FromAtom") \
or   (yasara.request=="MetalBinding_FromResidue") \
or   (yasara.request=="MetalBinding_FromMolecule") \
or   (yasara.request=="MetalBinding_FromObject"):
  
  # make a readable variable with the residue descriptors of the selected stretch
  metal=yasara.selection[1].listentry[0]
  # create an empty list
  allmetalbinding = []
  objnumberlist = []
  # grab the number of objects
  if yasara.request=="MetalBinding":
    objcount = yasara.selection[0].objects
    # loop over the selected objects
    for i in range(objcount):
      foldxmetalbinding.MetalBinding_Object(cachedir,foldxbin,yasara.selection[0].object[i],metal)
      metalbinding = foldxplotoutput.PlotMetalBinding_Object(cachedir,foldxbin,yasara.selection[0].object[i],metal)
      allmetalbinding += metalbinding
      pdbfilename = metal + "_MB_Object" + yasara.selection[0].object[i].number.inyas + ".pdb"
      objnumberlist += yasara.LoadPDB(os.path.join(cachedir,pdbfilename))
      
      yasara.ShowMessage("Superposing structures...")
      yasara.AlignObj(objnumberlist,yasara.selection[0].object[i].number.inyas,"MUSTANGPP")
      
  elif (yasara.request=="MetalBinding_FromAtom"):
    foldxmetalbinding.MetalBinding_Object(cachedir,foldxbin,yasara.selection[0].atom[0].object,metal)
    metalbinding = foldxplotoutput.PlotMetalBinding_Object(cachedir,foldxbin,yasara.selection[0].atom[0].object,metal)
    allmetalbinding = metalbinding
    pdbfilename = metal + "_MB_Object" + yasara.selection[0].atom[0].object.number.inyas + ".pdb"
    objnumberlist += yasara.LoadPDB(os.path.join(cachedir,pdbfilename))
    yasara.ShowMessage("Superposing structures...")
    yasara.AlignObj(objnumberlist,yasara.selection[0].atom[0].object.number.inyas,"MUSTANGPP")
  elif (yasara.request=="MetalBinding_FromResidue"):
    foldxmetalbinding.MetalBinding_Object(cachedir,foldxbin,yasara.selection[0].residue[0].object,metal)
    metalbinding = foldxplotoutput.PlotMetalBinding_Object(cachedir,foldxbin,yasara.selection[0].residue[0].object,metal)
    allmetalbinding = metalbinding
    pdbfilename = metal + "_MB_Object" + yasara.selection[0].atom[0].object.number.inyas + ".pdb"
    objnumberlist += yasara.LoadPDB(os.path.join(cachedir,pdbfilename))
    yasara.ShowMessage("Superposing structures...")
    yasara.AlignObj(objnumberlist,yasara.selection[0].residue[0].object.number.inyas,"MUSTANGPP")
  elif (yasara.request=="MetalBinding_FromMolecule"):
    foldxmetalbinding.MetalBinding_Object(cachedir,foldxbin,yasara.selection[0].molecule[0].object,metal)
    metalbinding = foldxplotoutput.PlotMetalBinding_Object(cachedir,foldxbin,yasara.selection[0].molecule[0].object,metal)
    allmetalbinding = metalbinding
    pdbfilename = metal + "_MB_Object" + yasara.selection[0].molecule[0].object.number.inyas + ".pdb"
    objnumberlist += yasara.LoadPDB(os.path.join(cachedir,pdbfilename))
    yasara.ShowMessage("Superposing structures...")
    yasara.AlignObj(objnumberlist,yasara.selection[0].molecule[0].object.number.inyas,"MUSTANGPP")
  elif (yasara.request=="MetalBinding_FromObject"):
    foldxmetalbinding.MetalBinding_Object(cachedir,foldxbin,yasara.selection[0].object[0],metal)
    metalbinding = foldxplotoutput.PlotMetalBinding_Object(cachedir,foldxbin,yasara.selection[0].object[0],metal)
    allmetalbinding += metalbinding
    pdbfilename = metal + "_MB_Object" + yasara.selection[0].object[0].number.inyas + ".pdb"
    objnumberlist += yasara.LoadPDB(os.path.join(cachedir,pdbfilename))
    yasara.ShowMessage("Superposing structures...")
    yasara.AlignObj(objnumberlist,yasara.selection[0].object[0].number.inyas,"MUSTANGPP")
  yasara.HideMessage()
  
  foldxutilities.PrintTable(allmetalbinding)
    
  # open the console
  yasara.Console("Open")    


# ******************
# * SequenceDetail *
# ******************
elif (yasara.request=="SequenceDetail") \
or   (yasara.request=="SequenceDetail_FromAtom") \
or   (yasara.request=="SequenceDetail_FromResidue") \
or   (yasara.request=="SequenceDetail_FromMolecule") \
or   (yasara.request=="SequenceDetail_FromObject"):

  # grab the number of objects
  if yasara.request=="SequenceDetail":
    objcount = yasara.selection[0].objects
    # loop over the selected objects
    for i in range(objcount):
      foldxsequencedetail.SequenceDetail_Object(cachedir,foldxbin,yasara.selection[0].object[i])
      sequencedetail = foldxplotoutput.PlotSequenceDetail_Object(cachedir,foldxbin,yasara.selection[0].object[i],i,"")
  elif (yasara.request=="SequenceDetail_FromAtom"):
    foldxsequencedetail.SequenceDetail_Object(cachedir,foldxbin,yasara.selection[0].atom[0].object)
    sequencedetail = foldxplotoutput.PlotSequenceDetail_Residue(cachedir,foldxbin,yasara.selection[0].atom[0].residue)
  elif (yasara.request=="SequenceDetail_FromResidue"):
    foldxsequencedetail.SequenceDetail_Object(cachedir,foldxbin,yasara.selection[0].residue[0].object)
    sequencedetail = foldxplotoutput.PlotSequenceDetail_Residue(cachedir,foldxbin,yasara.selection[0].residue)
  elif (yasara.request=="SequenceDetail_FromMolecule"):
    foldxsequencedetail.SequenceDetail_Object(cachedir,foldxbin,yasara.selection[0].molecule[0].object)
    sequencedetail = foldxplotoutput.PlotSequenceDetail_Object(cachedir,foldxbin,yasara.selection[0].molecule[0].object,1,"")
  elif (yasara.request=="SequenceDetail_FromObject"):
    foldxsequencedetail.SequenceDetail_Object(cachedir,foldxbin,yasara.selection[0].object[0])
    sequencedetail = foldxplotoutput.PlotSequenceDetail_Object(cachedir,foldxbin,yasara.selection[0].object[0],1,"")

  yasara.HideMessage()
  # open the console
  yasara.Console("Open")

# *************
# * Optimize *
# *************
elif (yasara.request=="Optimize") \
or   (yasara.request=="Optimize_FromAtom") \
or   (yasara.request=="Optimize_FromResidue") \
or   (yasara.request=="Optimize_FromMolecule") \
or   (yasara.request=="Optimize_FromObject"):
  
  # make a list for the objects
  objnumberlist = []
  
  # create an empty list
  if yasara.request=="Optimize":
    objectnumber=yasara.selection[0].object[0].number.inyas
    objcount = yasara.selection[0].objects
    # loop over the selected objects
    for i in range(objcount):
      objnumberlist += foldxoptimize.Optimize_Object(cachedir,foldxbin,yasara.selection[0].object[i])
  elif (yasara.request=="Optimize_FromAtom"):
    objectnumber=yasara.selection[0].atom[0].object.number.inyas
    objnumberlist += foldxoptimize.Optimize_Object(cachedir,foldxbin,yasara.selection[0].atom[0].object)
  elif (yasara.request=="Optimize_FromResidue"):
    objectnumber=yasara.selection[0].residue[0].object.number.inyas
    objnumberlist += foldxoptimize.Optimize_Object(cachedir,foldxbin,yasara.selection[0].residue[0].object)
  elif (yasara.request=="Optimize_FromMolecule"):
    objnumberlist += foldxoptimize.Optimize_Object(cachedir,foldxbin,yasara.selection[0].molecule[0].object)
  elif (yasara.request=="Optimize_FromObject"):
    objectnumber=yasara.selection[0].object[0].number.inyas
    objnumberlist += foldxoptimize.Optimize_Object(cachedir,foldxbin,yasara.selection[0].object[0])

  objnumbers = ""
  for objnumber in objnumberlist:
    objnumbers += " " + str(objnumber)
  
  ### Superpose FoldX models with original object
  #yasara.plugin.end(objnumbers,objectnumber)
  yasara.ShowMessage("Superposing structures...")
  yasara.AlignObj(objnumbers,objectnumber,"MUSTANGPP")
  yasara.HideMessage()

# ******************
# * Color Residues *
# ******************
elif (yasara.request=="ColorResidues") \
or   (yasara.request=="ColorResidues_FromAtom") \
or   (yasara.request=="ColorResidues_FromResidue") \
or   (yasara.request=="ColorResidues_FromMolecule") \
or   (yasara.request=="ColorResidues_FromObject"):

  # create an empty list
  allData = []
  if yasara.request=="ColorResidues":
    objcount = yasara.selection[0].objects
    # loop over the selected objects
    for i in range(objcount):
      foldxsequencedetail.SequenceDetail_Object(cachedir,foldxbin,yasara.selection[0].object[i])
      allData+=foldxplotoutput.PlotSequenceDetail_Object(cachedir,foldxbin,yasara.selection[0].object[i],i,"")
  elif (yasara.request=="ColorResidues_FromAtom"):
    foldxsequencedetail.SequenceDetail_Object(cachedir,foldxbin,yasara.selection[0].atom[0].object)
    allData+=foldxplotoutput.PlotSequenceDetail_Residue(cachedir,foldxbin,yasara.selection[0].atom[0].residue)
  elif (yasara.request=="ColorResidues_FromResidue"):
    foldxsequencedetail.SequenceDetail_Object(cachedir,foldxbin,yasara.selection[0].object)
    allData+=foldxplotoutput.PlotSequenceDetail_Residue(cachedir,foldxbin,yasara.selection[0].residue)
  elif (yasara.request=="ColorResidues_FromMolecule"):
    foldxsequencedetail.SequenceDetail_Object(cachedir,foldxbin,yasara.selection[0].molecule[0].object)
    #yasara.plugin.end(yasara.selection[0].residue)
    residues = []
    print(("res"+str(yasara.selection[0])))
    for res in yasara.selection[0].residue:
      residues.append(res)
    allData+=foldxplotoutput.PlotSequenceDetail_Residue(cachedir,foldxbin,residues)
  elif (yasara.request=="ColorResidues_FromObject"):
    for obj in yasara.selection[0].object:
      foldxsequencedetail.SequenceDetail_Object(cachedir,foldxbin,obj)
      allData+=foldxplotoutput.PlotSequenceDetail_Object(cachedir,foldxbin,obj,1,"")
  print(allData)

  foldxyastools.ColorAccordingToEnergy(allData)

  yasara.HideMessage()

# ****************
# * Alanine scan *
# ****************
elif (yasara.request=="Alascan") \
or   (yasara.request=="Alascan_FromAtom") \
or   (yasara.request=="Alascan_FromResidue") \
or   (yasara.request=="Alascan_FromMolecule") \
or   (yasara.request=="Alascan_FromObject"):
  
  # create an empty list
  allalascans = []
  # grab the number of objects
  if yasara.request=="Alascan":
    objcount = yasara.selection[0].objects
    # loop over the selected objects
    for i in range(objcount):
      foldxbuildmodel.Alascan_Object(cachedir,foldxbin,yasara.selection[0].object[i],i)
      AlascanEnergyOut = foldxplotoutput.PlotAlascan_Object(cachedir,foldxbin,yasara.selection[0].object[i],i)
      allalascans.append(AlascanEnergyOut)
  elif (yasara.request=="Alascan_FromAtom"):
    foldxbuildmodel.Alascan_Object(cachedir,foldxbin,yasara.selection[0].atom[0].object,1)
    AlascanEnergyOut = foldxplotoutput.PlotAlascan_Object(cachedir,foldxbin,yasara.selection[0].atom[0].object,1)
    allalascans.append(AlascanEnergyOut)
  elif (yasara.request=="Alascan_FromResidue"):
    foldxbuildmodel.Alascan_Object(cachedir,foldxbin,yasara.selection[0].residue[0].object,1)
    AlascanEnergyOut = foldxplotoutput.PlotAlascan_Object(cachedir,foldxbin,yasara.selection[0].residue[0].object,1)
    allalascans.append(AlascanEnergyOut)
  elif (yasara.request=="Alascan_FromMolecule"):
    foldxbuildmodel.Alascan_Object(cachedir,foldxbin,yasara.selection[0].molecule[0].object,1)
    AlascanEnergyOut = foldxplotoutput.PlotAlascan_Object(cachedir,foldxbin,yasara.selection[0].molecule[0].object,1)
    allalascans.append(AlascanEnergyOut)
  elif (yasara.request=="Alascan_FromObject"):
    foldxbuildmodel.Alascan_Object(cachedir,foldxbin,yasara.selection[0].object[0],1)
    AlascanEnergyOut = foldxplotoutput.PlotAlascan_Object(cachedir,foldxbin,yasara.selection[0].object[0],1)
    allalascans.append(AlascanEnergyOut)
  # stop showing messages on screen        
  yasara.HideMessage()
  # print all stabilities of all objects
  for AlascanEnergyOut in allalascans:
    yasara.write(AlascanEnergyOut)
    
  # open the console
  yasara.Console("Open")    

# ****************
# *   Dihedrals  *
# ****************
elif (yasara.request=="Dihedrals") \
or   (yasara.request=="Dihedrals_FromAtom") \
or   (yasara.request=="Dihedrals_FromResidue") \
or   (yasara.request=="Dihedrals_FromMolecule") \
or   (yasara.request=="Dihedrals_FromObject"):
  
  # create an empty list
  alldihedrals = []
  # grab the number of objects
  if yasara.request=="Dihedrals":
    objcount = yasara.selection[0].objects
    # loop over the selected objects
    for i in range(objcount):
      foldxdihedrals.Dihedrals_Object(cachedir,foldxbin,yasara.selection[0].object[i])
      DihedralsEnergyOut = foldxplotoutput.PlotDihedrals_Object(cachedir,foldxbin,yasara.selection[0].object[i])
      alldihedrals += DihedralsEnergyOut
  elif (yasara.request=="Dihedrals_FromAtom"):
    foldxdihedrals.Dihedrals_Object(cachedir,foldxbin,yasara.selection[0].atom[0].object)
    DihedralsEnergyOut = foldxplotoutput.PlotDihedrals_Residue(cachedir,foldxbin,yasara.selection[0].atom[0].residue)
    alldihedrals += DihedralsEnergyOut
  elif (yasara.request=="Dihedrals_FromResidue"):
    foldxdihedrals.Dihedrals_Object(cachedir,foldxbin,yasara.selection[0].residue[0].object)
    #DihedralsEnergyOut = foldxplotoutput.PlotDihedrals_Object(cachedir,foldxbin,yasara.selection[0].residue[0].object)
    DihedralsEnergyOut = foldxplotoutput.PlotDihedrals_Residue(cachedir,foldxbin,yasara.selection[0].residue[0])
    alldihedrals += DihedralsEnergyOut
  elif (yasara.request=="Dihedrals_FromMolecule"):
    foldxdihedrals.Dihedrals_Object(cachedir,foldxbin,yasara.selection[0].molecule[0].object)
    DihedralsEnergyOut = foldxplotoutput.PlotDihedrals_Object(cachedir,foldxbin,yasara.selection[0].molecule[0].object)
    alldihedrals += DihedralsEnergyOut
  elif (yasara.request=="Dihedrals_FromObject"):
    foldxdihedrals.Dihedrals_Object(cachedir,foldxbin,yasara.selection[0].object[0])
    DihedralsEnergyOut = foldxplotoutput.PlotDihedrals_Object(cachedir,foldxbin,yasara.selection[0].object[0])
    alldihedrals += DihedralsEnergyOut
  # stop showing messages on screen        
  yasara.HideMessage()

  # print all dihedrals of all objects
  foldxutilities.PrintTable(alldihedrals)
  # open the console
  yasara.Console("Open")    

#===================================================================================================================#

# stop the plugin, einde.
yasara.plugin.end()

#===================================================================================================================#

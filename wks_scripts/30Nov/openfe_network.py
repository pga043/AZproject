import openfe
from openfe import SmallMoleculeComponent
import networkx as nx
import glob, sys, os

import rdkit
from rdkit import Chem
from rdkit.Chem import AllChem
from rdkit.Chem import rdMolAlign
from rdkit.Chem import Recap,BRICS
from rdkit import DataStructs
from rdkit import RDLogger

import time
import math
import numpy as np
p = AllChem.ETKDGv2()
p.verbose = False

#------------------------------------------------------

## 1. open all sdf files ------------------------------
sdfs = []
for f in glob.glob('*.sdf') :
    tmp = Chem.SDMolSupplier(f, removeHs=False)
    sdfs.append(tmp)

molecules = []
for sdf in sdfs:
    smi = [Chem.MolToSmiles(x) for x in sdf]
    mol = [Chem.MolFromSmiles(x) for x in smi]
    molecules.append(mol[0])


## 2. ------------------------------------------
# generate dummy atom
du = Chem.MolFromSmiles('*')

all_smiles = []
all_smiles.clear()

##************************
## suppress rdkit warnings
RDLogger.DisableLog('rdApp.*')
##************************

## fragmenting using BRICSDecompose
for mol in molecules:
    frags = list(Chem.BRICS.BRICSDecompose(mol, minFragmentSize=1, keepNonLeafNodes=True, returnMols=False))
    ##if returnMols = False
    mols = [Chem.MolFromSmiles(x) for x in frags]
    ## replace dummy atoms with hydrogens
    molh = [AllChem.ReplaceSubstructs(x,du,Chem.MolFromSmiles('[H]'), True)[0] for x in mols]
    ## convert mols to smiles
    all_smiles.append([Chem.MolToSmiles(x) for x in molh])

##***************************
## enable default rdkit warnings again
RDLogger.EnableLog('rdApp.info')
##****************************

## 3. -----------------------------------------------------

## save fragemnets using similarity as a filter
## tanimoto similarity cut-off  = 0.35
## Dice threshold = 0.45
## the cut-offs are somewhat arbitrary for now.
x = 1
y = 1

cutoff = 0.35

os.makedirs('1-tmp', exist_ok=True)
##  all_smiles is a list of list so, need to use a nested loop
for i in all_smiles:
    for smi in i:
        mol   = Chem.MolFromSmiles(smi)
        query = mol
        mol   = Chem.AddHs(mol)
        ref   = molecules[all_smiles.index(i)] ## similarity wrt reference parent molecule
        #display(query, ref) # <<- jupyter-notebook
        #print(all_smiles.index(i))
        fp_q  = Chem.RDKFingerprint(query)
        fp_r  = Chem.RDKFingerprint(ref)
        sim   = DataStructs.TanimotoSimilarity(fp_r,fp_q)
        print(f'similarity = {sim}, reference sdf index = {all_smiles.index(i)}')
        print('========')
        if round(sim,1) > 0.9:
            #print(sim)
            #print('*************')
           ## save parent molecules with separate name
            ## save 3d sdf file
            with Chem.SDWriter('1-tmp/parent_' + str(y) + '.sdf') as writer:
                 ## set names for each molecule (to be used with openfe)
                 confIds = Chem.AllChem.EmbedMultipleConfs(mol, 1)
                 for confId in range(1):
                     Chem.AllChem.UFFOptimizeMolecule(mol, confId=confId)
                     writer.write(mol, confId=confId)
        else:
            if round(sim,2) > cutoff:
                ## save 3d sdf file
                with Chem.SDWriter('1-tmp/tmp_' +str(y) + '_' +str(x) + '.sdf') as writer:
                     confIds = Chem.AllChem.EmbedMultipleConfs(mol, 1)
                     for confId in range(1):
                         Chem.AllChem.UFFOptimizeMolecule(mol, confId=confId)
                         writer.write(mol, confId=confId)
                x += 1
            else:
                x += 0
    y +=1
    x = 1

## 4. -----------------------------------------------------------------
## align all the fragments to one parent as a reference molecule

## 4.1. open sdf files, one as a reference
##************************************************************
mol1 = Chem.SDMolSupplier('1-tmp/parent_1.sdf', removeHs=False)
for ref in mol1:
    ref1_smi = Chem.MolToSmiles(ref)
    ref1_mol = Chem.MolFromSmiles(ref1_smi)
    ref1_molh = Chem.AddHs(ref1_mol)

AllChem.EmbedMultipleConfs(ref1_molh)

## get MMFF parameters for the reference molecule
mmff_ref_param = AllChem.MMFFGetMoleculeProperties(ref1_molh)
ref_mol2 = ref1_molh

##*****************************************************

sdfs = []
names = []
for f in glob.glob("1-tmp/*.sdf"):
    tmp = Chem.SDMolSupplier(f, removeHs=False)
    sdfs.append(tmp[0])
    names.append(f[6:-4])

molecules = []
molecules.clear()
for sdf in sdfs:
    smi = Chem.MolToSmiles(sdf) #[Chem.MolToSmiles(x) for x in sdf]
    mol = Chem.MolFromSmiles(smi) #[Chem.MolFromSmiles(x) for x in smi]
    molecules.append(Chem.AddHs(mol))

st = time.time()

print('Generating 100 conformers for each molecule')
for mol in molecules[0:]:
    AllChem.EmbedMultipleConfs(mol, 100, p)

elapsed_time = time.time() - st
print('*************************')
print('Conformer generation time spent: ', time.strftime("%H:%M:%S", time.gmtime(elapsed_time)))
print('************************')

## 4.2. get MMFF parameters for each molecule
mmff_params = [AllChem.MMFFGetMoleculeProperties(mol) for mol in molecules]
#mmff_ref_param = mmff_params[0]
mmff_prob_params = mmff_params[0:]
#ref_mol2 = molecules[0]
prob_mols_2 = molecules[0:]


## 4.3. perform alignment and save the best scored conformer
print('Performing alignment and saving the conformer with best score only.')

with Chem.SDWriter('aligned.sdf') as writer:
    #writer.write(ref_mol2)
    pyO3A_score = []
    for idx, mol in enumerate(prob_mols_2):
        tempscore = []
        for cid in range(100):
            pyO3A = rdMolAlign.GetO3A(mol, ref_mol2, mmff_prob_params[idx], mmff_ref_param, cid, 0)
            pyO3A.Align()
            tempscore.append(pyO3A.Score())
        best = np.argmax(tempscore)
        ## save sdf file
        ## write names, to be used with network in openFE
        mol.SetProp('_Name', str(names[prob_mols_2.index(mol)]))
        mol.SetProp('ID', str(int(best)))
        writer.write(mol, confId=int(best))
        pyO3A_score.append(tempscore[best])

print('pyO3A_score: ' + str(pyO3A_score))


quit()

## 5. --------------------------------------------------------- 
##  openFE netwok generation 

sdfs = Chem.SDMolSupplier('aligned.sdf', removeHs=False)
#ligands = [SmallMoleculeComponent(sdf, name='junk') for sdf in sdfs]
ligands = [SmallMoleculeComponent(sdf) for sdf in sdfs]

mapper = openfe.LomapAtomMapper(max3d=1.0, element_change=True, threed=True)
scorer = openfe.lomap_scorers.default_lomap_score
network_planner = openfe.ligand_network_planning.generate_minimal_spanning_network

ligand_network = network_planner(
                 ligands = ligands[0:],
                 mappers = [mapper],
                 scorer = scorer)

nodes = list(ligand_network.nodes)

G = nx.Graph()
G



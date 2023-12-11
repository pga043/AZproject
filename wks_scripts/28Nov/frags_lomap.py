import rdkit
from rdkit import Chem
from rdkit.Chem import AllChem
from rdkit.Chem import Recap,BRICS
from rdkit import DataStructs
import lomap
import networkx as nx
import sys
from rdkit import RDLogger
import math
from rdkit.Chem import rdMolAlign
import numpy as np
p = AllChem.ETKDGv2()
p.verbose = False

#==========================================================

mol1 = Chem.SDMolSupplier(sys.argv[1], removeHs=False)
mol2 =  Chem.SDMolSupplier(sys.argv[2], removeHs=False)

sdfs = [mol1, mol2]
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
            with Chem.SDWriter('tmp/parent_' + str(y) + '.sdf') as writer:
                 writer.write(mol)
            ## save 3d sdf file
            with Chem.SDWriter('1-tmp/parent_' + str(y) + '.sdf') as writer:
                 ## set names for each molecule (to be used with openfe)
                 mol.SetProp('ID', 'parent_' + str(y))
                 confIds = Chem.AllChem.EmbedMultipleConfs(mol, 1)
                 for confId in range(1):
                     Chem.AllChem.UFFOptimizeMolecule(mol, confId=confId)
                     writer.write(mol, confId=confId)
        else:
            if round(sim,2) > cutoff:           
              ## save 2D sdf file fragments
                with Chem.SDWriter('tmp/tmp_' +str(y) + '_' +str(x) + '.sdf') as writer:
                     writer.write(mol)
                ## save 3d sdf file
                with Chem.SDWriter('1-tmp/tmp_' +str(y) + '_' +str(x) + '.sdf') as writer:
                     mol.SetProp('ID', 'tmp_' + str(y))
                     confIds = Chem.AllChem.EmbedMultipleConfs(mol, 1)
                     for confId in range(1):
                         Chem.AllChem.UFFOptimizeMolecule(mol, confId=confId)
                         writer.write(mol, confId=confId) 
                x += 1
            else:
                x += 0
    y +=1
    x = 1


## 4. -------------------------------------------------
## LOMAP network
db_mol = lomap.DBMolecules("tmp/", output=True)

strict, loose = db_mol.build_matrices()

nx_graph = db_mol.build_graph()


## 4.1 -----------------------------
## minimal spanning network (openFE)



## 5. -----------------------------------------------------
# get node values for the parent molecules
for key, value in nx.get_node_attributes(nx_graph, 'fname_comp').items():
    if 'parent_1.sdf' == value:
         node1 = key
    else:
        None
    if 'parent_2.sdf' == value:
         node2 = key
    else:
        None
        
try:
   node1, node2
except NameError:
     print('One or both of the main node(s) are undefined (sdf writing failed maybe).')
     print('Main nodes refers to the input molecules.')
     print('Exiting...')
     quit()


## 6. ------------------------------------------------------
## extract shortest path to connect parent molecules from lomap graph

## copy the original graph
G = nx_graph.copy()

## change the edge values 
for u, v, data in G.edges(data=True):
    #print(u, v, data)
    try:
       #data['weight'] = 1 / data['weight']         # one way
       data['weight'] = -math.log(data['similarity'])  # second way
    #except ZeroDivisionError:
    except ValueError:
        print("Warning ! You are dividing by zero ")
        print("or something wrong in the logarithm")

# get the shortest path from netowrk
try:
   sp = nx.shortest_path(G, source=node1, target=node2, weight='weight', method='dijkstra')
   #print(sp)
except nx.NetworkXNoPath:
    print(f"No path between {node1} and {node2}.")
    print('Exiting since no path found.')
    print('Exit....')
    quit()

print('Shortest path length: ')
print(nx.shortest_path_length(G, source=node1, target=node2, weight='weight', method='dijkstra'))

pathGraph = nx.path_graph(sp)


## 7. -----------------------------------------------------------
## get the path and align the 3D sdf files on a reference molecule
## save a sdf file for each path containing the node molecules


##************************************
## open a reference sturcture upon which others will be aligned
ref1 = Chem.SDMolSupplier('1-tmp/parent_1.sdf', removeHs=False)
for ref in ref1:
    ref1_smi = Chem.MolToSmiles(ref)
    ref1_mol = Chem.MolFromSmiles(ref1_smi)
    ref1_molh = Chem.AddHs(ref1_mol)

AllChem.EmbedMultipleConfs(ref1_molh)

## get MMFF parameters for the reference molecule
mmff_ref_param = AllChem.MMFFGetMoleculeProperties(ref1_molh)
ref_mol2 = ref1_molh
##*************************************

## :== do the alignment ==:
x = 0
for ea in pathGraph.edges():
    print(ea, nx_graph.edges[ea[0], ea[1]])
    #print(nx_graph.nodes[ea[0]]['fname_comp'])
    #print(nx_graph.nodes[ea[1]]['fname_comp'])
    print('paths: ' +str(nx_graph.nodes[ea[0]]['fname_comp']) + ' -> ' + str(nx_graph.nodes[ea[1]]['fname_comp']))

    mol1 = Chem.SDMolSupplier('1-tmp/' + str(nx_graph.nodes[ea[0]]['fname_comp']), removeHs=False)
    mol2 = Chem.SDMolSupplier('1-tmp/' + str(nx_graph.nodes[ea[1]]['fname_comp']), removeHs=False)
 
    sdfs = [mol1, mol2]
    
    molecules = []
    molecules.clear()
    for sdf in sdfs:
        smi = [Chem.MolToSmiles(x) for x in sdf]
        mol = [Chem.MolFromSmiles(x) for x in smi]
        molecules.append(Chem.AddHs(mol[0]))
 
    print('Generating 100 conformers for each molecule')
    for mol in molecules[0:]:
        AllChem.EmbedMultipleConfs(mol, 100, p)

    ## get MMFF parameters for each molecule
    mmff_params = [AllChem.MMFFGetMoleculeProperties(mol) for mol in molecules]
    #mmff_ref_param = mmff_params[0]
    mmff_prob_params = mmff_params[0:]
    #ref_mol2 = molecules[0]
    prob_mols_2 = molecules[0:]


    ## perform alignment and save the best scored conformer
    print('Performing alignment and saving the conformer with best score only.')

    with Chem.SDWriter(f'aligned_{x}.sdf') as writer:
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
             mol.SetProp('_Name', str(nx_graph.nodes[ea[prob_mols_2.index(mol)]]['fname_comp']))
             mol.SetProp('ID', str(int(best)))
             writer.write(mol, confId=int(best))
             pyO3A_score.append(tempscore[best])

    #print('pyO3A_score: ' + str(pyO3A_score))
    x +=1

print('Finished. Next step OpenFE')
quit()


import lomap
import networkx as nx

import os, sys, glob, subprocess, re
import math
import numpy as np
import matplotlib.pyplot as plt

import rdkit
from rdkit import Chem
from rdkit.Chem import AllChem
from rdkit.Chem import Recap,BRICS
from rdkit import DataStructs
from rdkit import RDLogger
from rdkit.Chem import rdMolAlign

p = AllChem.ETKDGv2()
p.verbose = False

import openfe
from openfe import SmallMoleculeComponent

##=============================================================
##*************************************************************
input_dir = sys.argv[1]   ## input directory containing sdf files
ref_mol   = sys.argv[2]   ## reference molecule for alignment
work_dir  = sys.argv[3]   ## working directory to save output files

'''
Required inputs:
1. input directory containing the sdf files for LOMAP netowrk generation
2. Reference molecule that will be used for alignment (required)

Output directories or files:
a. frags : fragments from RDKit
b. fromLomap : Lomap edges with good scores, aligned to a given reference structure.
c. forOpenfe : Lomap edges with bad scores and their fragments, all aligned to a given reference structure.
'''
##*************************************************************
##==============================================================

## load all the sdf files and generate lomap network

db_mol = lomap.DBMolecules(input_dir, output=True, name=str(work_dir)+'/outx')
strict, loose = db_mol.build_matrices()
nx_graph = db_mol.build_graph()

## get edges whose score is less than cutoff
lomap_cutoff = 0.4 ## I defined this cutoff to filter out bad edges from lomap graph

for_frag = [] ## these molecules will be fragmented later on with RDKit
for_frag.clear()
for edge in nx_graph.edges(data=True):
    score = nx_graph.edges[[edge[0], edge[1]]]['similarity']
    if score > lomap_cutoff:
        print('transformation: ' + str(nx_graph.nodes[edge[0]]['fname_comp']) + ' <--> ' + str(nx_graph.nodes[edge[1]]['fname_comp']))
    else:
        for_frag.append(nx_graph.nodes[edge[0]]['fname_comp'])
        for_frag.append(nx_graph.nodes[edge[1]]['fname_comp'])
        #print('Fragmentation required: ' + str(nx_graph.nodes[edge[0]]['fname_comp']) + ' <--> ' + str(nx_graph.nodes[edge[1]]['fname_comp']))

##----------------------------------

## generate a directory to store fragments:
os.makedirs(str(work_dir)+'/frags', exist_ok=True)

## get unique molecules from edges which have low lomap scores
molecules = []
molecules.clear()
for m in range(len(list(set(for_frag)))):
    tmp = Chem.SDMolSupplier(f'{input_dir}/{list(set(for_frag))[m]}', removeHs=False)
    smi = Chem.MolToSmiles(tmp[0])
    mol = Chem.MolFromSmiles(smi)
    molecules.append(mol)


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

## save fragments using similarity as a filter
## tanimoto similarity cut-off  = 0.35 (arbitrary)
## or Dice threshold = 0.45
## the cut-offs are somewhat arbitrary for now.
x = 1
y = 1

## Different parent might generate same fragments so,
fps = [] ## save all fragment's fingerprint to filter duplicate fragments from different parents
fps.clear()

cutoff = 0.30  ## arbitrary similarity cutoff, fragments with lower similarity than this won't be saved.

##  all_smiles is a list of list so, need to use a nested loop
for i in all_smiles:
    ref   = molecules[all_smiles.index(i)] ## similarity wrt reference parent molecule
    fp_r  = Chem.RDKFingerprint(ref)
    #fps.append(fp_r)
    for smi in i:
        mol   = Chem.MolFromSmiles(smi)
        query = mol
        mol   = Chem.AddHs(mol)        
        fp_q  = Chem.RDKFingerprint(query)
        sim   = DataStructs.TanimotoSimilarity(fp_r,fp_q)
        #print(f'similarity = {sim}, reference sdf index = {all_smiles.index(i)}')
        #print('========')
        if cutoff < round(sim,1) < 0.9: ## it means it's a fragment
            ## check whether we have a similar fragment from another parent
            frag_fps = []
            frag_fps.clear()
            for fp in range(len(fps)):
                frag_fps.append(DataStructs.TanimotoSimilarity(fps[fp], fp_q))
            if any(z==1 for z in frag_fps):
                x+=0
            else:
                fps.append(fp_q)
                ## save fragment molecules             
                with Chem.SDWriter(str(work_dir)+'/frags/tmp_' + str(y) + '_' +str(x) + '.sdf') as writer:
                    ## settings names in the output sdf file
                    ## important for OpenFE, nodes are referred to as their names in the sdf files not the sdf_file_name
                    mol.SetProp('_Name', 'tmp_' + str(y) + '_' + str(x)) 
                    mol.SetProp('ID', 'tmp_' + str(y) + '_' + str(x))
                    confIds = Chem.AllChem.EmbedMultipleConfs(mol, 1)
                    for confId in range(1):
                        Chem.AllChem.UFFOptimizeMolecule(mol, confId=confId)
                        writer.write(mol, confId=confId)
                        x += 1
        else:
            x += 0
    y +=1
    x = 1

## file name: tmp_y_x.sdf
## y = parent, x = fragment

##=============================================

## generate and align 3D structures of the edges with good lomap scores

##***************************************
## open a reference sturcture upon which others will be aligned
## preserving input coordinates (must be 3D structure data file)
ref1 = Chem.SDMolSupplier(ref_mol, removeHs=False)
for ref in ref1:
    #ref1_smi = Chem.MolToSmiles(ref)
    #ref1_mol = Chem.MolFromSmiles(ref1_smi)
    ref1_molh = ref #Chem.AddHs(ref1_mol)

AllChem.EmbedMultipleConfs(ref1_molh, 1)

## get MMFF parameters for the reference molecule
mmff_ref_param = AllChem.MMFFGetMoleculeProperties(ref1_molh)
ref_mol2 = ref1_molh

##=================================================


## LOMAP edges ====>>
## write out one aligned (to a given reference) sdf files for each lomap edge with good score.
## These can then be directly used with OpenFE.
## change the reference molecule if you're not satisfied with the alignment (RMSD).

os.makedirs(str(work_dir)+'/fromLomap', exist_ok=True)

for edge in nx_graph.edges(data=True):
    score = nx_graph.edges[[edge[0], edge[1]]]['similarity']
    if score > lomap_cutoff:
        print('transformation: ' + str(nx_graph.nodes[edge[0]]['fname_comp']) + ' <--> ' + str(nx_graph.nodes[edge[1]]['fname_comp']))
        
        mol1 = Chem.SDMolSupplier(str(input_dir)+'/'+str(nx_graph.nodes[edge[0]]['fname_comp']), removeHs=False)
        mol2 = Chem.SDMolSupplier(str(input_dir)+'/'+str(nx_graph.nodes[edge[1]]['fname_comp']), removeHs=False)
        
        sdfs = [mol1, mol2]
        molecules = []
        molecules.clear()
        for sdf in sdfs:
            smi = [Chem.MolToSmiles(x) for x in sdf]
            mol = [Chem.MolFromSmiles(x) for x in smi]
            molecules.append(Chem.AddHs(mol[0]))

        confs = 100 ## number of conformers to generate
        print(f'Generating {confs} conformers for each molecule')
        for mol in molecules[0:]:
            AllChem.EmbedMultipleConfs(mol, confs, p)

        ## get MMFF parameters for each molecule
        mmff_params = [AllChem.MMFFGetMoleculeProperties(mol) for mol in molecules]
        #mmff_ref_param = mmff_params[0]
        mmff_prob_params = mmff_params[0:]
        #ref_mol2 = molecules[0]
        prob_mols_2 = molecules[0:]


        ## perform alignment and save the best scored conformer
        print('Performing alignment and saving the conformer with best score only.')

        with Chem.SDWriter(str(work_dir)+'/fromLomap/' +str(nx_graph.nodes[edge[0]]['fname_comp']) +'_' + str(nx_graph.nodes[edge[1]]['fname_comp'])) as writer:
            #writer.write(ref_mol2) ## save the reference if you want to look at the alignment
            pyO3A_score = []
            rmsd = []
            for idx, mol in enumerate(prob_mols_2):
                tempscore = []
                temprmsd = []
                for cid in range(confs):
                    pyO3A = rdMolAlign.GetO3A(mol, ref_mol2, mmff_prob_params[idx], mmff_ref_param, cid, 0)
                    temprmsd.append(pyO3A.Align())
                    tempscore.append(pyO3A.Score())
                best = np.argmax(tempscore)
                rmsd.append(temprmsd[best])
                ## save sdf file
                ## write names, to be used with network in openFE
                mol.SetProp('_Name', str(nx_graph.nodes[edge[prob_mols_2.index(mol)]]['fname_comp']))
                mol.SetProp('ID', str(int(best)))
                mol.SetProp('Open3D alignment rmsd wrt reference: ', str(round(rmsd[idx], 4)))
                writer.write(mol, confId=int(best))
                pyO3A_score.append(tempscore[best])
                #print('rmsd of best(saved) conformer wrt reference: ' + str(round(rmsd[idx], 4)))
                
                
        
        #print('pyO3A_score: ' + str(pyO3A_score))
        print('===================================')

print('Finished LOMAP edges. Next step OpenFE')

##============================================================

## openFE requires 3D sdf files and also they need to be aligned
## align fragments and their parent molecules to the same reference and save as individual sdf file

os.makedirs(str(work_dir)+'/forOpenfe', exist_ok=True)

file_names = []
sdfs = []
sdfs.clear()
for edge in nx_graph.edges(data=True):
    score = nx_graph.edges[[edge[0], edge[1]]]['similarity']
    if score < lomap_cutoff:
        print('poor edge: ' + str(nx_graph.nodes[edge[0]]['fname_comp']) + ' <--> ' + str(nx_graph.nodes[edge[1]]['fname_comp']))
        
        for i in range(2):  ## the length of an edge is 2.
            if str(nx_graph.nodes[edge[i]]['fname_comp']) not in file_names[0:]:
                file_names.append(str(nx_graph.nodes[edge[i]]['fname_comp']))
                sdfs.append(Chem.SDMolSupplier(str(input_dir)+'/'+str(nx_graph.nodes[edge[i]]['fname_comp']), removeHs=False))

if not sdfs:
    print("No bad edges found")
    print("OpenFE network generation is not required.")
    print("Finished, go for RBFE with LOMAP network only.")
    quit()

for f in glob.glob(str(work_dir)+"/frags/*.sdf"):
    sdf = Chem.SDMolSupplier(f, removeHs=False)
    sdfs.append(sdf)
    file_names.append(re.search('tmp_(.+?).sdf', f).group(0)) ## I just want to get 'tmp_x_y.sdf' and store it

#print(file_names)
print(f'------------------------------')

molecules = []
molecules.clear()
for sdf in sdfs:
    smi = [Chem.MolToSmiles(x) for x in sdf]
    mol = [Chem.MolFromSmiles(x) for x in smi]
    molecules.append(Chem.AddHs(mol[0]))

confs = 100 ## number of conformers to generate
print(f'Generating {confs} conformers for each molecule')
for mol in molecules[0:]:
    AllChem.EmbedMultipleConfs(mol, confs, p)

## get MMFF parameters for each molecule
mmff_params = [AllChem.MMFFGetMoleculeProperties(mol) for mol in molecules]
#mmff_ref_param = mmff_params[0]
mmff_prob_params = mmff_params[0:]
#ref_mol2 = molecules[0]
prob_mols_2 = molecules[0:]


## perform alignment and save the best scored conformer
print('Performing alignment and saving the conformer with best score only.')

#with Chem.SDWriter(f'forOpenfe/.sdf') as writer:
    #writer.write(ref_mol2)

pyO3A_score = []
rmsd = []
for idx, mol in enumerate(prob_mols_2):   
    tempscore = []
    temprmsd = []
    for cid in range(confs):
        pyO3A = rdMolAlign.GetO3A(mol, ref_mol2, mmff_prob_params[idx], mmff_ref_param, cid, 0)
        #pyO3A.Align()
        temprmsd.append(pyO3A.Align())
        tempscore.append(pyO3A.Score())
    best = np.argmax(tempscore)
    rmsd.append(temprmsd[best])
    ## save sdf file
    print(f'writing sdf file: {file_names[prob_mols_2.index(mol)]}')
    with Chem.SDWriter(f'{work_dir}/forOpenfe/{file_names[prob_mols_2.index(mol)]}') as writer:
        #writer.write(ref_mol2) ## save the reference if you want to look at the alignment
        ## write names, to be used with network in openFE
        mol.SetProp('_Name', str(file_names[prob_mols_2.index(mol)]))
        mol.SetProp('ID', str(int(best)))
        mol.SetProp('Open3D alignment rmsd wrt reference: ', str(round(rmsd[idx], 4)))
        writer.write(mol, confId=int(best))
        
    #print('rmsd of best(saved) conformer wrt reference: ' + str(round(rmsd[idx], 4)))
    pyO3A_score.append(tempscore[best])

#print('pyO3A_score: ' + str(pyO3A_score))

print('Finished with aligned file, go to OpenFE for network generation.')

##=============================================

# openFE maximal_netwok

sdfs = []
for f in glob.glob(str(work_dir)+'/forOpenfe/*.sdf'):
    sdf = Chem.SDMolSupplier(f, removeHs=False)
    sdfs.append(sdf[0])

ligands = [SmallMoleculeComponent(sdf) for sdf in sdfs]

mapper = openfe.LomapAtomMapper(max3d=100.0, element_change=True)
scorer = openfe.lomap_scorers.default_lomap_score
network_planner = openfe.ligand_network_planning.generate_maximal_network


ligand_network = network_planner(
ligands=ligands[0:],
     mappers=[mapper],
     scorer=scorer
)

##==================================================

## openFE network is a frozen set and, I am not sure how to work with that. so,
## I create a new graph and we can now do wahtever we want to this graph.
## The new graph contains node names and lomap scores from the OpenFE network.

G = nx.Graph()

## low lomap score == considered as good edge weight by networkx or netowrkx edges are closer
## high lomap score == considered as bad edge weight by networkx or networkx edges are further 

for edge in ligand_network.edges:
    try:
        G.add_edge(edge.componentA.name, edge.componentB.name, weight=round(-math.log(edge.annotations['score']), 4), lomap_score=round(edge.annotations['score'], 4))
    except ValueError: ## if the edge score is zero, assign a high score (here randomly I chose 5)
        G.add_edge(edge.componentA.name, edge.componentB.name, weight=5, lomap_score=round(edge.annotations['score'], 4))


##=================================
## from LOMAP network get poor edges since those are the edges that we are interseted in.
## for each edge, define one node as source and the other as target.
## Then we can find all paths (within a cutoff) connecting the source and target from the OpenFE maximal network.
## From all the paths, we can then define a shortest path which connects the source and target and also has good lomap scores along the way.
## Similarity scores were converted to edge weights in graph G as: W = -log(similarity)
## Function to find the shortest path = minimum{[(-logW1) + (-logW2) + .... + (-logWn) ]/n} where, n = number of edges

for edge in nx_graph.edges(data=True): 
    score = nx_graph.edges[[edge[0], edge[1]]]['similarity']
    if score < lomap_cutoff:
        source = nx_graph.nodes[edge[0]]['fname_comp']
        target = nx_graph.nodes[edge[1]]['fname_comp']
        
        os.makedirs(f'{work_dir}/{source}_{target}', exist_ok=True)

        all_paths = list(nx.all_simple_paths(G, source=source, target=target, cutoff=3)) ##
        
        score = []
        score.clear()
        final_path = []
        final_path.clear()
        for path in range(len(all_paths)):
            #print(all_paths[path])
            pathGraph = nx.path_graph(all_paths[path])
            final_path.append(pathGraph)
            weights = []
            weights.clear()
            for ea in pathGraph.edges():
                weights.append(G.edges[ea[0], ea[1]]['weight'])
            
            score.append(round(sum(weights) / pathGraph.number_of_edges(), 4))
            #print('final score: ' + str(round(sum(weights) / pathGraph.number_of_edges(), 4)))

        print(f'{source}, {target}, score: {min(score, key=abs)}, {final_path[score.index(min(score, key=abs))]}')
        
        for ea in final_path[score.index(min(score, key=abs))].edges(data=True):
            if G.edges[ea[0], ea[1]]['lomap_score'] < 0.2:
                print(ea[0], '<-->', ea[1], '=> LOMAP score', G.edges[ea[0], ea[1]]['lomap_score'], '**poor edge**')
            else:
                print(ea[0], '<-->', ea[1], '=> LOMAP score', G.edges[ea[0], ea[1]]['lomap_score'])
            file_write = (f'cat {work_dir}/forOpenfe/{ea[0]} {work_dir}/forOpenfe/{ea[1]}  > {work_dir}/{source}_{target}/{ea[0]}_{ea[1]}')
            subprocess.check_call(file_write, shell=True)
        print('===================')
            

print(f'*****************')
print(f'Successfully Finished')

quit()


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
        #display(query, ref)
        #print(all_smiles.index(i))
        fp_q  = Chem.RDKFingerprint(query)
        fp_r  = Chem.RDKFingerprint(ref)
        sim   = DataStructs.TanimotoSimilarity(fp_r,fp_q)
        #print(sim)
        #print('========')
        if round(sim,1) > 0.9:
            #print(sim)
            #print('*************')
           ## save parent molecules with separate name
            sdfFile2d = open('tmp/parent_' + str(y) + '.sdf', 'w')
            writer = Chem.SDWriter(sdfFile2d)
            writer.write(mol)
            #sdfFile2d.close()
            ## save 3d sdf file
            sdfFile3d = open('1-tmp/parent_' + str(y) + '.sdf', 'w')
            confIds = Chem.AllChem.EmbedMultipleConfs(mol, 1)
            for confId in range(1):
                Chem.AllChem.UFFOptimizeMolecule(mol, confId=confId)
                writer = Chem.SDWriter(sdfFile3d)
                writer.write(mol, confId=confId)
        else:
            if round(sim,2) > cutoff:           
              ## save 2D sdf file fragments
                sdfFile2d = open('tmp/tmp_' +str(y) + '_' +str(x) + '.sdf', 'w')
                writer = Chem.SDWriter(sdfFile2d)
                writer.write(mol)
                #sdfFile2d.close()
                ## save 3d sdf file
                sdfFile3d = open('1-tmp/junk_' +str(y) + '_' +str(x) + '.sdf', 'w')
                confIds = Chem.AllChem.EmbedMultipleConfs(mol, 1)
                for confId in range(1):
                    Chem.AllChem.UFFOptimizeMolecule(mol, confId=confId)
                    writer = Chem.SDWriter(sdfFile3d)
                    writer.write(mol, confId=confId) 
                    #sdfFile3d.close()
                x += 1
            else:
                x += 0
    y +=1
    x = 1


## 4. -------------------------------------------------

db_mol = lomap.DBMolecules("tmp/", output=True)

strict, loose = db_mol.build_matrices()

nx_graph = db_mol.build_graph()


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
        
node1, node2


### check if the parent nodes are connected
#try:
#   test = nx.has_path(nx_graph, source=node1, target=node2)
#   print(test)
#   print('************')
#except nx.NetworkXNoPath:  
#    print('The parent molecules are not connected by lomap')

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
    quit()

print('Shortest path length: ')
print(nx.shortest_path_length(G, source=node1, target=node2, weight='weight', method='dijkstra'))

pathGraph = nx.path_graph(sp)

for ea in pathGraph.edges():
    print(ea, nx_graph.edges[ea[0], ea[1]])
    #print(nx_graph.nodes[ea[0]]['fname_comp'])
    #print(nx_graph.nodes[ea[1]]['fname_comp'])
    print('paths: ' +str(nx_graph.nodes[ea[0]]['fname_comp']) + ' -> ' + str(nx_graph.nodes[ea[1]]['fname_comp']))


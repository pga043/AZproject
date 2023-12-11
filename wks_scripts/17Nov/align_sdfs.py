import rdkit
from rdkit import Chem
from rdkit.Chem import AllChem
from rdkit.Chem import rdMolAlign
import numpy as np
import sys
p = AllChem.ETKDGv2()
p.verbose = False


## 1. open sdf files (one should be te reference file upon which to align)
mol1 = Chem.SDMolSupplier('junk_9_11/1-tmp/parent_1.sdf', removeHs=False)
mol2 =  Chem.SDMolSupplier('junk_9_11/1-tmp/junk_1_2.sdf', removeHs=False)

sdfs = [mol1, mol2]
molecules = []

for sdf in sdfs:
    smi = [Chem.MolToSmiles(x) for x in sdf]
    mol = [Chem.MolFromSmiles(x) for x in smi]
    molecules.append(Chem.AddHs(mol[0]))


print('Generating 100 conformers for each molecule')
for mol in molecules[0:]:
    AllChem.EmbedMultipleConfs(mol, 100, p)


## 2. get MMFF parameters for each molecule
mmff_params = [AllChem.MMFFGetMoleculeProperties(mol) for mol in molecules]
mmff_ref_param = mmff_params[0]
mmff_prob_params = mmff_params[1:]
ref_mol2 = molecules[0]
prob_mols_2 = molecules[1:]


## 3. perform alignment and save the best scored conformer
print('Performing alignment and saving the conformer with best score only.')

with Chem.SDWriter('aligned.sdf') as writer:
    writer.write(ref_mol2)
    pyO3A_score = []
    for idx, mol in enumerate(prob_mols_2):
        tempscore = []
        for cid in range(100):
            pyO3A = rdMolAlign.GetO3A(mol, ref_mol2, mmff_prob_params[idx], mmff_ref_param, cid, 0)
            pyO3A.Align()
            tempscore.append(pyO3A.Score())
        best = np.argmax(tempscore)
        ## save sdf file
        writer.write(mol, confId=int(best))
        pyO3A_score.append(tempscore[best])

print('pyO3A_score: ' + str(pyO3A_score))

## 4. Done
quit()



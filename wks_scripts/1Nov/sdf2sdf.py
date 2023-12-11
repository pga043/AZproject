from rdkit import Chem
from rdkit import DataStructs

# read the sdf files
file_name = '/projects/mai/users/kvnq006_parveen/data/fragmentation/literature/mobley_jctc_2023/MALT1/malt1.sdf'

supply = Chem.SDMolSupplier(file_name) # name the path of the SDF file 

mols = [x for x in  supply]

x = 0
for mol in mols:
    #print(type(mol))
    name = x
    sdf = Chem.SDWriter('malt1_' + str(name) + '.sdf')
    sdf.write(mol)
    sdf.close()
    x += 1



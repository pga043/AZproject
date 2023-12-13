import rdkit
from rdkit import Chem

suppl = Chem.SDMolSupplier('p38_vina_naive.sdf')

for mol in suppl:
     name = mol.GetProp("_Name")
     molh = Chem.AddHs(mol, explicitOnly=True, addCoords=True)
     with Chem.SDWriter(f'000/{name}.sdf') as writer:
         writer.write(molh)

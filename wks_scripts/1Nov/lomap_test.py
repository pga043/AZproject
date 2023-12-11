import lomap

# Generate the molecule database starting from a directory containing .mol2 files

input_dir = '.'

db_mol = lomap.DBMolecules(input_dir, output=True)

    #More graphing options:
    # Use the complete radial graph option. The ligand with the most structural similarity to all of the others will be picked as the 'lead compounds' and used as the central compound.
    #db_mol = lomap.DBMolecules(input_dir, output=True, radial=True)

    # Use a radial graph with a manually specified hub compound
#    db_mol = lomap.DBMolecules(input_dir, output=True, radial=True, hub=filename.mol2)

    # Use a radial graph with a manually specified hub compound and fast graphing option
    #the fast graphing option create the initial graph by connecting the hub ligand with the possible surrounding ligands and add surrounding edges based on the similarities accoss surrounding nodes
#    db_mol = lomap.DBMolecules("python string pointing to a directory with mol2 files", output=True, radial=True, hub=filename.mol2, fast=True)

# Calculate the similarity matrix betweeen the database molecules. Two molecules are generated
# related to the scrict rule and loose rule 

strict, loose = db_mol.build_matrices()

# Generate the NetworkX graph and output the results
nx_graph = db_mol.build_graph() 

# Calculate the Maximum Common Subgraph (MCS) between 
# the first two molecules in the molecule database 
# ignoring hydrogens and depicting the mapping in a file
    
MC = lomap.MCS.getMapping(db_mol[0].getMolecule(), db_mol[1].getMolecule(), hydrogens=False, fname='mcs.png')


# Alchemical transformation are usually performed between molecules with
# the same charges. However, it is possible to allow this transformation
# manually setting the electrostatic score for the whole set of molecules 
# producing a connected graph. The electrostatic scrore must be in the 
# range [0,1]


db_mol = lomap.DBMolecules(input_dir, output=True, ecrscore=0.1)
strict, loose = db_mol.build_matrices()
nx_graph = db_mol.build_graph() 

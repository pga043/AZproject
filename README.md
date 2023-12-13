# AZproject

#********************#
lomap_and_openfe.ipynb
#********************#

This jupyter notebook performs the following tasks:

1. generate lomap network for given inputs (should be more than two and preferablly sdfs, mol2 files are also fine)\
   (convery mol2 files to sdf files using openbabel. Rdkit can't handle mol2 files)

2. split the lomap netwok into good and bad edges based on lomap scores (bad < 0.4 < good)

3. fragment the molecules associated with bad edges using RDKit BRICS and \
   save fragments with tanimoto similarity (0.30 < cutoff > 0.9)

4. open a reference molecule (3D sdf file) upon which other molecules will be aligned

5. Take the lomap edges with good scores and align each edge molecule upon the reference (in step 4) \
   and save sdf files.

6. Now, take the molecules with bad lomap scores and all of their fragments and align them \
   to the reference structure in step 4 and save each molecule/fragment as sdf files.

7. Generate a maximal netwok using OpenFE for the molecules with bad lomap score and their fragments

8. For each molecule pair (with bad lomap score in step 1), find a shortest path within a cutoff (maximum number of edges)
   to connect them (via other molecules or fragments).

9. write out sdf files for each pair and it's corresponding shortest path.


Output directories or files:
a. frags : fragments from RDKit
b. fromLomap : Lomap edges with good scores, aligned to a given reference structure.
c. forOpenfe : Lomap edges with bad scores and their fragments, all aligned to a given reference structure.


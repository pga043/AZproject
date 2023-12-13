#!/bin/bash

## if you want to use the python script instead of the notebook
##
## usage: python lomap_and_openfe.py input_dir_containing_sdfs reference_mol_for_alignment output_dir

python lomap_and_openfe.py malt1 malt1/malt1_2.sdf malt1/

python lomap_and_openfe.py tnks tnks/mol_2.sdf tnks

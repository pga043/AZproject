import openfe
from rdkit import Chem
from rdkit.Chem import AllChem
from rdkit.Chem import Draw

ligand = Chem.SDMolSupplier("hybrid.sdf", removeHs=False)

ligands = [openfe.SmallMoleculeComponent.from_rdkit(mol) for mol in ligand]

# Create network from the two molecules
from openfe.setup.ligand_network_planning import generate_radial_network
from openfe.setup import LomapAtomMapper

network = generate_radial_network(ligands=ligands[1:],
                                  central_ligand=ligands[0],
                                  mappers=[LomapAtomMapper(threed=True, element_change=False),])

edges = [edge for edge in network.edges]

transform_edges = [edge for edge in network.edges]
print("molecule A smiles: ", transform_edges[0].componentA.smiles)
print("molecule B smiles: ", transform_edges[0].componentB.smiles)
print("map between molecule A and B: ", transform_edges[0].componentA_to_componentB)

edges = [edge for edge in network.edges]
print(f"edge 0 molecule 1: {edges[0].componentA.name}")
print(f"edge 0 molecule 2: {edges[0].componentB.name}")
print(f"edge 0 mapping: {edges[0].componentA_to_componentB}")
frag_name = edges[0].componentB.name


#######------------- solvation/neutralization ----------------------------

# First let's define the Protein and Solvent Components which we will be using
from openfe import SolventComponent, ProteinComponent
from openff.units import unit

#protein = ProteinComponent.from_pdb_file('inputs/181L_mod_capped_protonated.pdb')

# Note: the distance from the solute to add water is not defined here but in the
# the relevant RBFE solver method
solvent = SolventComponent(positive_ion='Na', negative_ion='Cl',
                           neutralize=False, ion_concentration=0*unit.molar)

####------------------------------------


# Extract the relevant edge for the benzene -> phenol transform in the radial graph
parent_2_frag = [edge for edge in network.edges if edge.componentB.name == frag_name][0]
print(parent_2_frag)

##-------------- alchemical system ----------------------------

# Let's create the four ChemicalSystems
from openfe import ChemicalSystem

#benzene_complex = ChemicalSystem({'ligand': benz_to_phenol.componentA,
#                                  'solvent': solvent,
#                                  'protein': protein,})
parent_solvent = ChemicalSystem({'ligand': edges[0].componentA,
                                  'solvent': solvent,})

#phenol_complex = ChemicalSystem({'ligand': benz_to_phenol.componentB,
#                                 'solvent': solvent,
#                                 'protein': protein,})
frag_solvent = ChemicalSystem({'ligand': edges[0].componentB,
                                 'solvent': solvent,})

##-----------------------------------------------------------------


#------------------- MD parameters ------------------------------
# Settings can be accessed from the various classes

from openfe.protocols.openmm_rfe.equil_rfe_settings import (
    SystemSettings, SolvationSettings, AlchemicalSettings,
    OpenMMEngineSettings, AlchemicalSamplerSettings,
    IntegratorSettings, SimulationSettings
)

system = SystemSettings(nonbonded_cutoff=1.2 * unit.nanometer)


from openfe.protocols.openmm_rfe import RelativeHybridTopologyProtocol

rbfe_settings = RelativeHybridTopologyProtocol.default_settings()

from pprint import pp

# Parameters can also be overriden after creation
# In this case, we'll reduce the equilibration length to 0.01 * nanosecond
# and the production to 0.05 * nanosecond in order to reduce the costs of
# running this notebook (in practice values of 2 and 5 nanoseconds
# respectively would be most appropriate)

rbfe_settings.alchemical_sampler_settings.sampler_method = 'sams'
rbfe_settings.alchemical_sampler_settings.n_replicas = 21
rbfe_settings.alchemical_sampler_settings.n_repeats = 3

rbfe_settings.alchemical_settings.lambda_windows = 21

rbfe_settings.engine_settings.compute_platform = 'CUDA'

rbfe_settings.integrator_settings.n_restart_attempts = 30

rbfe_settings.simulation_settings.equilibration_length = 1000 * unit.picosecond
rbfe_settings.simulation_settings.production_length = 5 * unit.nanosecond

pp(rbfe_settings.simulation_settings)

print(pp(rbfe_settings.simulation_settings))

##---------------------------------
##========== Print settings =============
print(SystemSettings())
print('-------------------------------------')
print(rbfe_settings.alchemical_settings)
print('-------------------------------------')
print(OpenMMEngineSettings)
print('-------------------------------------')
print(rbfe_settings.alchemical_sampler_settings)
print('-------------------------------------')
print(rbfe_settings.simulation_settings)
print('-------------------------------------')
print('-------------------------------------')

##=======================================
##---------------------------------

# Create RBFE Protocol class
rbfe_transform = RelativeHybridTopologyProtocol(
    settings=rbfe_settings
)

solvent_dag = rbfe_transform.create(
    stateA=parent_solvent, stateB=frag_solvent,
    mapping={'ligand': parent_2_frag},
)

list(solvent_dag.protocol_units)


# solvent dry-run
solvent_unit = list(solvent_dag.protocol_units)[0]

solvent_unit.run(dry=True, verbose=True)

from gufe.protocols import execute_DAG
import pathlib

# Next the solvent state transformation
solvent_path = pathlib.Path('./solvent')
solvent_path.mkdir()


### ------------- run and analysis ---------------------
solvent_dag_results = execute_DAG(solvent_dag, scratch_basedir=solvent_path, shared_basedir=solvent_path, keep_shared=True, keep_scratch=True)


# Get the complex and solvent results
#complex_results = rbfe_transform.gather([complex_dag_results])
solvent_results = rbfe_transform.gather([solvent_dag_results])

#print(f"Complex dG: {complex_results.get_estimate()}, err {complex_results.get_uncertainty()}")
print(f"Solvent dG: {solvent_results.get_estimate()}, err {solvent_results.get_uncertainty()}")

quit()



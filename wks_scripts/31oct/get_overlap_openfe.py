import netCDF4 as nc
from openmmtools.multistate import MultiStateReporter
from openfe.protocols.openmm_utils.multistate_analysis import MultistateEquilFEAnalysis
from openff.units import unit
import matplotlib.pyplot as plt
import sys


#reporter = MultiStateReporter("simulation.nc")

reporter = MultiStateReporter(sys.argv[1])

ana = MultistateEquilFEAnalysis(reporter=reporter, sampling_method="repex", result_units=unit.kilocalorie_per_mole)

ana.analyzer._unbiased_decorrelated_N_l

mat = ana.free_energy_overlaps["matrix"]

fig, axes = plt.subplots(2, 1, figsize=(4, 4.5), gridspec_kw={"height_ratios": [15, 1], "width_ratios": [3]})
ax, cax = axes
im = ax.matshow(mat, vmin=0.0, vmax=1.0, interpolation="nearest", cmap="plasma")
ax.tick_params(labelsize=16)
ax.set_title("Overlap", fontsize=16)
fig.colorbar(im, cax=cax, fraction=0.046, pad=0.04, orientation="horizontal")
cax.tick_params(labelsize=16)


#plt.show()

#plt.save('')

pl.show()



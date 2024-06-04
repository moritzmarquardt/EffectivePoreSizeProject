import MembraneAnalysisToolbox.EffectivePoreSizeAnalysis as EPSA
import numpy as np
import matplotlib.pyplot as plt
from scipy.stats import gaussian_kde

for i in [2]:

    # path = "/bigpool/users/ac130484/project/finished_sim/hex/poresize" + str(i) + "nm_NVT/"
    path = "/bigpool/data/projects/Carbon_pores_Sofia/correctdensity/hex_NVT/" + str(i) + "nm_NVT/"
    sim = 1

    Analysis = EPSA.EffectivePoreSizeAnalysis(
                topology_file = path + 'simulation_' + str(sim) + '/' + 'topol.tpr',
                trajectory_file = path + 'EPS_analysis/' + 'traj_simulation_' + str(sim) + '.xtc',
                membrane_resnames = ['C'],
                solvent_resnames = ['HEX', 'DOD'],
                y_middle = 35,
                y_range = 10,
                verbose = True
            )

    # Analysis.analyseDensity()
    Analysis.analyseDensityNormalised()
plt.show()
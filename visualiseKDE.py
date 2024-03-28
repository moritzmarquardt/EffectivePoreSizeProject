import MembraneAnalysisToolbox.EffectivePoreSizeAnalysis as EPSA
import numpy as np
import matplotlib.pyplot as plt
from scipy.stats import gaussian_kde

path = "/bigpool/data/projects/Carbon_pores_Sofia/correctdensity/hex_NVT/2nm_NVT/"
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
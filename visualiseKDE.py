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

membrane_pos = Analysis.membrane_atom_positions
solvent_pos = Analysis.solvent_atom_positions

filtered_indices_mem = np.where(
    (membrane_pos[:, :, 2] >= (Analysis.z_min + Analysis.z_max)/2 - 5) & 
    (membrane_pos[:, :, 2] <= (Analysis.z_min + Analysis.z_max)/2 + 5)
)
filtered_indices_sol = np.where(
    (solvent_pos[:, :, 2] >= (Analysis.z_min + Analysis.z_max)/2 - 5) & 
    (solvent_pos[:, :, 2] <= (Analysis.z_min + Analysis.z_max)/2 + 5)
)

filtered_positions_mem = membrane_pos[filtered_indices_mem[0], filtered_indices_mem[1], 0:2]
filtered_positions_sol = solvent_pos[filtered_indices_sol[0], filtered_indices_sol[1], 0:2]


kde_mem = gaussian_kde(filtered_positions_mem[::50,:].T)
kde_sol = gaussian_kde(filtered_positions_sol[::50,:].T)

# Evaluate the KDE on a grid
xmin = min(filtered_positions_mem[:, 0])
xmax = max(filtered_positions_mem[:, 0])
xn = 100
ymin = min(filtered_positions_mem[:, 1])
ymax = max(filtered_positions_mem[:, 1])
yn = 100
x = np.linspace(xmin, xmax, xn)
y = np.linspace(ymin, ymax, yn)
X, Y = np.meshgrid(x, y)
Z_mem = np.reshape(kde_mem([X.ravel(), Y.ravel()]), X.shape)
Z_sol = np.reshape(kde_sol([X.ravel(), Y.ravel()]), X.shape)
Z_sol_norm = Z_sol / np.max(Z_sol)
Z_mem_norm = Z_mem / np.max(Z_mem)
Z = Z_mem_norm + Z_sol_norm
plt.figure()
plt.pcolormesh(X, Y, Z, shading='gouraud')
plt.colorbar()
plt.axis('equal')
plt.show()
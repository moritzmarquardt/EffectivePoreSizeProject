import numpy as np
import matplotlib.pyplot as plt
import MDAnalysis as mda
import time

start = time.time()

# several time frames in the xtc file and tpr as topology file
# make xtc from gro : gmx trjconv -f /hugepool/data/projects/traj_dt10ns.gro -o /bigpool/users/st166545/TPS/240130kri_new_sim/traj.xtc -skip 1
# use one frame of gro or topology file as topology
# MDAnalysis uses Angstroms for distances, not nanometers. The GROMACS software, which produces .gro and .xtc files, uses nanometers. 
# So, when you read a .gro or .xtc file with MDAnalysis, the coordinates are automatically converted from nanometers to Angstroms. 

# u = mda.Universe('/bigpool/users/st166545/TPS/240130kri_new_sim/traj_firstframe.gro', '/bigpool/users/st166545/TPS/240130kri_new_sim/traj.xtc')
u = mda.Universe('/bigpool/users/st166545/TPS/4nm/topol.tpr', '/bigpool/users/st166545/TPS/4nm/traj.xtc')
# u = mda.Universe('/bigpool/data/projects/Carbon_pores_Sofia/kri_with_vacuum_test/hex_18_3_2_n_vacuum/simulation_1/topol.tpr', '/bigpool/users/st166545/TPS/kri_with_vacuum_test_skip10.xtc')
print(u)
print("timeframes analysed: " + str(u.trajectory.n_frames))


c_atoms = u.select_atoms('resname C')
hex_atoms = u.select_atoms('resname HEX')
dod_atoms = u.select_atoms('resname DOD')

c_atomic_positions = np.zeros((u.trajectory.n_frames,len(c_atoms),3))
dod_atomic_positions = np.zeros((u.trajectory.n_frames,len(dod_atoms),3))
hex_atomic_positions = np.zeros((u.trajectory.n_frames,len(hex_atoms),3))

for ts in u.trajectory:
    c_atomic_positions[ts.frame,:,:] = c_atoms.positions
    hex_atomic_positions[ts.frame,:,:] = hex_atoms.positions
    dod_atomic_positions[ts.frame,:,:] = dod_atoms.positions


print("pos aquire done after: " + str(time.time()-start))

# Set constraints on y and z directions
y_min, y_max = 40, 50  # Replace with your desired range
z_min, z_max = 125, 275  # Replace with your desired range


# Histograms for constraint analysis
'''plt.figure()
plt.hist(c_atomic_positions[:,:,0].flatten(),bins=50)
plt.title('Histogram for x-axis')
plt.figure()
plt.hist(c_atomic_positions[:,:,1].flatten(),bins=50)
plt.axvline(x=y_min, color='r', linestyle='--')  # y_min line
plt.axvline(x=y_max, color='r', linestyle='--')  # y_max line
plt.title('Histogram for y-axis')
plt.figure()
plt.hist(c_atomic_positions[:,:,2].flatten(),bins=50)
plt.axvline(x=z_min, color='r', linestyle='--')  # z_min line
plt.axvline(x=z_max, color='r', linestyle='--')  # z_max line
plt.title('Histogram for z-axis')
plt.show()'''




def filter_positions(positions, y_min, y_max, z_min, z_max):
    """
    Filter positions based on y and z constraints.

    Parameters:
    positions (ndarray): Array of positions with shape (timesteps, atom_amount, 3).
    y_min (float): lower y constraint.
    y_max (float): upper y constraint.
    z_min (float): lower z constraint.
    z_max (float): upper z constraint.

    Returns:
    ndarray: Array of filtered x positions.
    """
    filtered_indices = np.where(
        (positions[:, :, 1] >= y_min) & (positions[:, :, 1] <= y_max) &
        (positions[:, :, 2] >= z_min) & (positions[:, :, 2] <= z_max)
    )
    filtered_x_positions = positions[filtered_indices[0], filtered_indices[1], 0].flatten()
    return filtered_x_positions

c_filtered_x_positions = filter_positions(c_atomic_positions, y_min, y_max, z_min, z_max)
dod_filtered_x_positions = filter_positions(dod_atomic_positions, y_min, y_max, z_min, z_max)
hex_filtered_x_positions = filter_positions(hex_atomic_positions, y_min, y_max, z_min, z_max)
print("filtering done after: " + str(time.time()-start))


# Compute histograms
bins = 100
c_hist, c_bin_edges = np.histogram(c_filtered_x_positions, density=1, bins=bins)
dod_hex_hist, dod_hex_bin_edges= np.histogram(np.append(dod_filtered_x_positions,hex_filtered_x_positions), density=1, bins=bins)
print("hist calc done after: " + str(time.time()-start))


#compute true pore size:
def calculate_pore_size(c_hist, c_bin_edges, dod_hex_hist, dod_hex_bin_edges):
    """
    Calculate the lower and upper edges of the pore size distribution based on the averaging mathod.

    Parameters:
    c_hist (numpy.ndarray): Histogram of the concentration values.
    c_bin_edges (numpy.ndarray): Bin edges of the concentration histogram.
    dod_hex_hist (numpy.ndarray): Histogram of the DOD hex values.
    dod_hex_bin_edges (numpy.ndarray): Bin edges of the DOD hex histogram.

    Returns:
    tuple: A tuple containing the average lower edge and average upper edge of the pore size distribution.
    """

    first_zero_bin = np.where(c_hist == 0)[0][0]
    first_zero_middle = (c_bin_edges[first_zero_bin] + c_bin_edges[first_zero_bin + 1]) / 2
    first_non_zero_bin = np.where(dod_hex_hist > 0)[0][0]
    first_non_zero_middle = (dod_hex_bin_edges[first_non_zero_bin] + dod_hex_bin_edges[first_non_zero_bin + 1]) / 2
    avrg_lower_edge = np.abs(first_zero_middle + first_non_zero_middle) / 2

    last_zero_bin = np.where(c_hist == 0)[0][-1]
    last_zero_middle = (c_bin_edges[last_zero_bin] + c_bin_edges[last_zero_bin + 1]) / 2
    last_non_zero_bin = np.where(dod_hex_hist > 0)[0][-1]
    last_non_zero_middle = (dod_hex_bin_edges[last_non_zero_bin] + dod_hex_bin_edges[last_non_zero_bin + 1]) / 2
    avrg_upper_edge = np.abs(last_zero_middle + last_non_zero_middle) / 2

    return avrg_lower_edge, avrg_upper_edge

lower_edge, upper_edge = calculate_pore_size(c_hist, c_bin_edges, dod_hex_hist, dod_hex_bin_edges)
effective_pore_size = np.abs(lower_edge - upper_edge)  # in Angstroms
print("true pore size calculated after: " + str(time.time()-start))
print("true pore size in nm: " + str(effective_pore_size/10))

# Plot all histograms in one plot
plt.figure()
plt.plot(c_bin_edges[:-1], c_hist, label='C', linestyle='-', marker='o', markersize=3)
plt.plot(dod_hex_bin_edges[:-1], dod_hex_hist, label='DOD & HEX', linestyle='-', marker='o', markersize=3)
x1 = lower_edge
x2 = upper_edge
plt.axvline(x=x1, color='r', linestyle='--')
plt.axvline(x=x2, color='r', linestyle='--')
plt.xlabel('X-axis in Angstroms')
plt.ylabel('Frequency')
plt.title('Histogram Line Plot along the X-axis with Y and Z constraints for Atoms')
plt.grid(True)
plt.legend()
plt.show()

effective_pore_sizes = [2.942, 3.880]
init_pore_sizes = [3, 4]

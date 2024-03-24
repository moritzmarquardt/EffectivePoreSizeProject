import MembraneAnalysisToolbox.EffectivePoreSizeAnalysis as EPSA
import numpy as np
import matplotlib.pyplot as plt
import cv2
from scipy.stats import gaussian_kde
from skimage.measure import CircleModel, ransac

'''path = "/bigpool/data/projects/Carbon_pores_Sofia/correctdensity/hex_NVT/2nm_NVT/"

Analysis = EPSA.EffectivePoreSizeAnalysis(
            topology_file = path + 'simulation_' + str(1) + '/' + 'topol.tpr',
            trajectory_file = path + 'EPS_analysis/' + 'traj_simulation_' + str(1) + '.xtc',
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

np.save("filtered_positions_mem.npy", filtered_positions_mem)
np.save("filtered_positions_sol.npy", filtered_positions_sol)

print("safed")'''

filtered_positions_mem = np.load("filtered_positions_mem.npy")
filtered_positions_sol = np.load("filtered_positions_sol.npy")


kde_mem = gaussian_kde(filtered_positions_mem[::100,:].T)
kde_sol = gaussian_kde(filtered_positions_sol[::100,:].T)

# Evaluate the KDE on a grid
xmin = min(filtered_positions_mem[:, 0])
xmax = max(filtered_positions_mem[:, 0])
xn = int(np.ceil(xmax - xmin))
ymin = min(filtered_positions_mem[:, 1])
ymax = max(filtered_positions_mem[:, 1])
yn = int(np.ceil(xn * (ymax - ymin) / (xmax - xmin)))
x = np.linspace(xmin, xmax, xn)
y = np.linspace(ymin, ymax, yn)
X, Y = np.meshgrid(x, y)
Z_mem = np.reshape(kde_mem([X.ravel(), Y.ravel()]), X.shape)
Z_sol = np.reshape(kde_sol([X.ravel(), Y.ravel()]), X.shape)
Z_sol_norm = Z_sol / np.max(Z_sol)
Z_mem_norm = Z_mem / np.max(Z_mem)
Z = Z_mem_norm + Z_sol_norm
Z2 = np.abs(Z_mem_norm - Z_sol_norm)

'''plt.figure()
plt.pcolormesh(X, Y, Z, shading='gouraud')
plt.colorbar()
plt.axis('equal')'''


plt.figure()
Z_8bit = cv2.normalize(Z2, None, 0, 255, cv2.NORM_MINMAX).astype(np.uint8)
canny = cv2.Canny(Z_8bit,10,200)
plt.imshow(canny)
circle = None
approx_radius = 10
circles = cv2.HoughCircles(canny, cv2.HOUGH_GRADIENT, dp=1, minDist=50, param1=10, param2 = 10, minRadius=approx_radius-3, maxRadius=approx_radius+3)
print(circles)
circle = circles[0][0]
if circle is not None:
    c = plt.Circle((circle[0], circle[1]), circle[2], color='red', fill=False)
    plt.gca().add_artist(c)



plt.figure()
plt.pcolormesh(X, Y, Z2, shading='gouraud')
plt.colorbar()
plt.axis('equal')
if circle is not None:
    c = plt.Circle((circle[0], circle[1]), circle[2], color='red', fill=False)
    plt.gca().add_artist(c)

plt.show()


"""# Calculate the range of the x and y data
x_range = np.ptp(filtered_positions[:, 0])
y_range = np.ptp(filtered_positions[:, 1])

# Calculate the number of bins to use to ensure square bins
bins = int(np.ceil(max(x_range, y_range) / min(x_range, y_range) * 100))

h, xedges, yedges, image =  plt.hist2d(filtered_positions[:, 0], filtered_positions[:, 1], bins = bins)
plt.colorbar()
# plt.show()
image_8bit = cv2.normalize(h, None, 0, 255, cv2.NORM_MINMAX).astype(np.uint8)
# cv2.imshow("Image", image_8bit)
# cv2.waitKey(0)
# cv2.destroyAllWindows()

circles = cv2.HoughCircles(image_8bit, cv2.HOUGH_GRADIENT, dp=1, minDist=10, param1=50, param2 = 50)
print(circles, )

# If at least one circle was found
if circles is not None:
    # Convert the (x, y) coordinates and radius of the circles to integers
    circles = np.round(circles[0, :]).astype("int")


    # Loop over the circles and draw them on the image
    for (x, y, r) in circles:
        # cv2.circle(image_8bit, (x, y), r, (0, 255, 0), 4)
        
        # Draw the circle on the histogram plot
        circle = plt.Circle((49.5, 34.5), 11.8, color='red', fill=False)
        plt.gca().add_artist(circle)

# Display the image with circles
#plt.figure()
#plt.imshow(image_8bit, cmap='gray')
plt.axis('equal')
plt.show()"""
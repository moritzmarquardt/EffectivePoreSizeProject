import MembraneAnalysisToolbox.EffectivePoreSizeAnalysis as mat

Analysis3nm = mat.EffectivePoreSizeAnalysis( 
    topology_file = '/bigpool/data/projects/Carbon_pores_Sofia/correctdensity/hex/3nm_1musonly/simulation_1/topol.tpr', 
    trajectory_file = '/bigpool/data/projects/Carbon_pores_Sofia/correctdensity/hex/3nm_1musonly/simulation_1/traj.xtc',
    membrane_resnames = ['C'],
    solvent_resnames = ['HEX', 'DOD'],
    y_middle = 40, 
    y_range = 10, 
    z_min = 150, 
    z_max = 300,
    verbose = True
)
print("calculated effective pore size (3nm) is: " + str(Analysis3nm.calculate_effective_pore_size() / 10))
Analysis3nm.plot()
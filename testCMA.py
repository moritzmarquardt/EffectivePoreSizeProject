import MembraneAnalysisToolbox.EffectivePoreSizeAnalysis as mat

Analysis2nm = mat.EffectivePoreSizeAnalysis( 
    topology_file = '/bigpool/data/projects/Carbon_pores_Sofia/correctdensity/hex/2nm_1musonly/simulation_1/topol.tpr', 
    trajectory_file = '/bigpool/data/projects/Carbon_pores_Sofia/correctdensity/hex/2nm_1musonly/simulation_1/traj.xtc',
    membrane_resnames = ['C'],
    solvent_resnames = ['HEX', 'DOD'],
    y_middle = 35, 
    y_range = 10, 
    z_min = 250, 
    z_max = 400
)
print("calculated effective pore size (2nm) is: " + str(Analysis2nm.calculate_effective_pore_size() / 10))


Analysis3nm = mat.EffectivePoreSizeAnalysis( 
    topology_file = '/bigpool/data/projects/Carbon_pores_Sofia/correctdensity/hex/3nm_1musonly/simulation_1/topol.tpr', 
    trajectory_file = '/bigpool/data/projects/Carbon_pores_Sofia/correctdensity/hex/3nm_1musonly/simulation_1/traj.xtc',
    membrane_resnames = ['C'],
    solvent_resnames = ['HEX', 'DOD'],
    y_middle = 40, 
    y_range = 10, 
    z_min = 150, 
    z_max = 300
)
print("calculated effective pore size (3nm) is: " + str(Analysis3nm.calculate_effective_pore_size() / 10))


Analysis4nm = mat.EffectivePoreSizeAnalysis( 
    topology_file = '/bigpool/data/projects/Carbon_pores_Sofia/correctdensity/hex/4nm_1musonly/simulation_1/topol.tpr', 
    trajectory_file = '/bigpool/data/projects/Carbon_pores_Sofia/correctdensity/hex/4nm_1musonly/simulation_1/traj.xtc',
    membrane_resnames = ['C'],
    solvent_resnames = ['HEX', 'DOD'],
    y_middle = 45, 
    y_range = 10, 
    z_min = 125, 
    z_max = 275
)
print("calculated effective pore size (4nm) is: " + str(Analysis4nm.calculate_effective_pore_size() / 10))


Analysis6nm = mat.EffectivePoreSizeAnalysis( 
    topology_file = '/bigpool/data/projects/Carbon_pores_Sofia/correctdensity/hex/6nm_1musonly/simulation_1/topol.tpr', 
    trajectory_file = '/bigpool/data/projects/Carbon_pores_Sofia/correctdensity/hex/6nm_1musonly/simulation_1/traj.xtc',
    membrane_resnames = ['C'],
    solvent_resnames = ['HEX', 'DOD'],
    y_middle = 50, 
    y_range = 10, 
    z_min = 200, 
    z_max = 325
)
print("calculated effective pore size (6nm) is: " + str(Analysis6nm.calculate_effective_pore_size() / 10))


Analysis2nm.plot()
Analysis3nm.plot()
Analysis4nm.plot()
Analysis6nm.plot()
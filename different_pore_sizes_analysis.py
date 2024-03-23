import MembraneAnalysisToolbox.EffectivePoreSizeAnalysis as mat

Analysis2nm = mat.EffectivePoreSizeAnalysis( 
    topology_file = '/bigpool/data/projects/Carbon_pores_Sofia/correctdensity/hex/2nm_1musonly/simulation_1/topol.tpr', 
    trajectory_file = '/bigpool/data/projects/Carbon_pores_Sofia/correctdensity/hex/2nm_1musonly/simulation_1/traj.xtc',
    membrane_resnames = ['C'],
    solvent_resnames = ['HEX', 'DOD'],
    y_middle = 35, 
    y_range = 10, 
)
print("calculated effective pore size (2nm) is: " + str(Analysis2nm.calculate_effective_pore_size() / 10))
print("relation between pore size and effective pore size: " + str(Analysis2nm.calculate_effective_pore_size()/10/2*100) + "%")
Analysis2nm.analyseConstraints()


Analysis3nm = mat.EffectivePoreSizeAnalysis( 
    topology_file = '/bigpool/data/projects/Carbon_pores_Sofia/correctdensity/hex/3nm_1musonly/simulation_1/topol.tpr', 
    trajectory_file = '/bigpool/data/projects/Carbon_pores_Sofia/correctdensity/hex/3nm_1musonly/simulation_1/traj.xtc',
    membrane_resnames = ['C'],
    solvent_resnames = ['HEX', 'DOD'],
    y_middle = 40, 
    y_range = 10, 
)
print("calculated effective pore size (3nm) is: " + str(Analysis3nm.calculate_effective_pore_size() / 10))
print("relation between pore size and effective pore size: " + str(Analysis3nm.calculate_effective_pore_size()/10/3*100) + "%")


Analysis4nm = mat.EffectivePoreSizeAnalysis( 
    topology_file = '/bigpool/data/projects/Carbon_pores_Sofia/correctdensity/hex/4nm_1musonly/simulation_1/topol.tpr', 
    trajectory_file = '/bigpool/data/projects/Carbon_pores_Sofia/correctdensity/hex/4nm_1musonly/simulation_1/traj.xtc',
    membrane_resnames = ['C'],
    solvent_resnames = ['HEX', 'DOD'],
    y_middle = 45, 
    y_range = 10, 
)
print("calculated effective pore size (4nm) is: " + str(Analysis4nm.calculate_effective_pore_size() / 10))
print("relation between pore size and effective pore size: " + str(Analysis4nm.calculate_effective_pore_size()/10/4*100) + "%")


Analysis6nm = mat.EffectivePoreSizeAnalysis( 
    topology_file = '/bigpool/data/projects/Carbon_pores_Sofia/correctdensity/hex/6nm_1musonly/simulation_1/topol.tpr', 
    trajectory_file = '/bigpool/data/projects/Carbon_pores_Sofia/correctdensity/hex/6nm_1musonly/simulation_1/traj.xtc',
    membrane_resnames = ['C'],
    solvent_resnames = ['HEX', 'DOD'],
    y_middle = 50, 
    y_range = 10, 
)
print("calculated effective pore size (6nm) is: " + str(Analysis6nm.calculate_effective_pore_size() / 10))
print("relation between pore size and effective pore size: " + str(Analysis6nm.calculate_effective_pore_size()/10/6*100) + "%")


path = "/bigpool/data/projects/Carbon_pores_Sofia/correctdensity/hex_NVT/3nm/hex_18_3_2_n_NVT/" # new nvt sim , z middle 30
Analysis3nm_newNVT = mat.EffectivePoreSizeAnalysis( 
    topology_file = path + 'topol.tpr', 
    trajectory_file = '/bigpool/users/st166545/EffectivePoreSizeProject/simdata/new3nmNVThex_traj.xtc',
    membrane_resnames = ['C'],
    solvent_resnames = ['HEX', 'DOD'],
    y_middle = 30, 
    y_range = 10, 
)
# Analysis3nm_newNVT.analyseConstraints()
print("calculated effective pore size (3nm) is: " + str(Analysis3nm_newNVT.calculate_effective_pore_size(strategy = "intersection") / 10) + "nm")
print("relation between pore size and effective pore size: " + str(Analysis3nm_newNVT.calculate_effective_pore_size(strategy = "intersection")/10/3*100) + "%")


Analysis2nm.plot()
Analysis3nm.plot()
Analysis4nm.plot()
Analysis6nm.plot()
Analysis3nm_newNVT.plot()
import numpy as np
import opinf
import os
from matplotlib import pyplot as plt
import normaopinf
import normaopinf.readers
import normaopinf.calculus
import romtools

if __name__ == "__main__":
    # Load in snapshots
    cur_dir = os.getcwd()
    solution_id = 1
    displacement_snapshots,times = normaopinf.readers.load_displacement_csv_files(solution_directory=cur_dir,solution_id=solution_id,skip_files=1)

    # Identify which DOFs are free
    free_dofs = normaopinf.readers.get_free_dofs(solution_directory=cur_dir,solution_id=solution_id)

    # Set values = 0 if DOFs are fixed
    displacement_snapshots[free_dofs[:,:]==False] = 0.

    my_energy_truncater = romtools.vector_space.utils.EnergyBasedTruncater(1. - 1.e-5)
    #my_energy_truncater = romtools.vector_space.utils.BasisSizeTruncater(30)

    trial_space = romtools.VectorSpaceFromPOD(snapshots=displacement_snapshots[:],
                                              truncater=my_energy_truncater,
                                              shifter = None,
                                              orthogonalizer=romtools.vector_space.utils.EuclideanL2Orthogonalizer(),
                                              scaler = romtools.vector_space.utils.NoOpScaler())

    uhat = romtools.rom.optimal_l2_projection(displacement_snapshots,trial_space)
    u_ddots,times_dummy = normaopinf.readers.load_acceleration_csv_files(solution_directory=cur_dir,solution_id=solution_id,skip_files=1)
    # Set values = 0 if DOFs are fixed.
    u_ddots[free_dofs[:,:]==False] = 0.
    uhat_ddots = romtools.rom.optimal_l2_projection(u_ddots*1.,trial_space)
    l2solver = opinf.lstsq.L2Solver(regularizer=1e-11)
    opinf_model = opinf.models.ContinuousModel("A",solver=l2solver)

    #opinf_model.fit(states=uhat, ddts=uhat_ddots,inputs=reduced_stacked_sideset_snapshots)
    opinf_model.fit(states=uhat, ddts=uhat_ddots)
    K = -opinf_model.A_.entries
    f = np.zeros(K.shape[0])

    vals_to_save = {}
    vals_to_save["basis"] = trial_space.get_basis()
    vals_to_save["K"] = K
    vals_to_save["f"] = f

    np.savez('opinf-operator',**vals_to_save)

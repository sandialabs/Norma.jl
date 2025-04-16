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
    solution_id = 3
    displacement_snapshots,times = normaopinf.readers.load_displacement_csv_files(solution_directory=cur_dir,solution_id=solution_id,skip_files=1)
    
    # Identify which DOFs are free
    free_dofs = normaopinf.readers.get_free_dofs(solution_directory=cur_dir,solution_id=solution_id)
    
    # Set values = 0 if DOFs are fixed
    displacement_snapshots[free_dofs[:,:]==False] = 0.
    
    #Get sideset snapshots
    sidesets = ["nodeset_positive_y-y","sideset_negative_y"]
    sideset_snapshots = normaopinf.readers.load_sideset_displacement_csv_files(solution_directory=cur_dir,sidesets=sidesets,solution_id=solution_id,skip_files=1)
    
    # Create an energy-based truncater 
    #tolerance = 1.e-5
    #my_energy_truncater = romtools.vector_space.utils.EnergyBasedTruncater(1. - tolerance)
    my_energy_truncater = romtools.vector_space.utils.BasisSizeTruncater(4)

    # Now load in sidesets and create reduced spaces
    # Note that I construct a separate basis for each x,y,z component. This isn't necessary 
    ss_tspace = {}
    reduced_sideset_snapshots = {}
    for sideset in sidesets:
        if sideset_snapshots[sideset].shape[0] ==  1:
          ss_tspace[sideset] = romtools.VectorSpaceFromPOD(snapshots=sideset_snapshots[sideset],
                                              truncater=my_energy_truncater,
                                              shifter = None,
                                              orthogonalizer=romtools.vector_space.utils.EuclideanL2Orthogonalizer(),
                                              scaler = romtools.vector_space.utils.NoOpScaler())
          # Compute L2 orthogonal projection onto trial spaces
          reduced_sideset_snapshots[sideset] = romtools.rom.optimal_l2_projection(sideset_snapshots[sideset],ss_tspace[sideset]) 
        else:
          comp_trial_space = []
          for i in range(0,3):
            tspace = romtools.VectorSpaceFromPOD(snapshots=sideset_snapshots[sideset][i:i+1],
                                              truncater=my_energy_truncater,
                                              shifter = None,
                                              orthogonalizer=romtools.vector_space.utils.EuclideanL2Orthogonalizer(),
                                              scaler = romtools.vector_space.utils.NoOpScaler())
            comp_trial_space.append(tspace)
          ss_tspace[sideset] = romtools.CompositeVectorSpace(comp_trial_space)
          reduced_sideset_snapshots[sideset] = romtools.rom.optimal_l2_projection(sideset_snapshots[sideset],ss_tspace[sideset]) 


    ## Stack sidesets into one matrix for OpInf 
    reduced_stacked_sideset_snapshots = None
    for sideset in sidesets:
        if reduced_stacked_sideset_snapshots is None:
            reduced_stacked_sideset_snapshots = reduced_sideset_snapshots[sideset]*1.
        else: 
            reduced_stacked_sideset_snapshots = np.append(reduced_stacked_sideset_snapshots,reduced_sideset_snapshots[sideset],axis=0)
    
    my_energy_truncater = romtools.vector_space.utils.BasisSizeTruncater(20)
    
    # Create trial space for displacement vector
    # Note again that I construct a separate basis for each x,y,z component. This isn't necessary 
    trial_spaces = []
    for i in range(0,3):
      trial_space = romtools.VectorSpaceFromPOD(snapshots=displacement_snapshots[i:i+1],
                                              truncater=my_energy_truncater,
                                              shifter = None,
                                              orthogonalizer=romtools.vector_space.utils.EuclideanL2Orthogonalizer(),
                                              scaler = romtools.vector_space.utils.NoOpScaler())
      trial_spaces.append(trial_space)
    
    trial_space = romtools.CompositeVectorSpace(trial_spaces)
    
    # Compute L2 orthogonal projection onto trial spaces 
    uhat = romtools.rom.optimal_l2_projection(displacement_snapshots,trial_space)
    #u_ddots = normaopinf.calculus.d2dx2(displacement_snapshots,times)
    u_ddots,times_dummy = normaopinf.readers.load_acceleration_csv_files(solution_directory=cur_dir,solution_id=solution_id,skip_files=1)
    # Set values = 0 if DOFs are fixed
    u_ddots[free_dofs[:,:]==False] = 0.
    uhat_ddots = romtools.rom.optimal_l2_projection(u_ddots*1.,trial_space)
   
    # Construct an opinf "AB" model (linear in the state and linear in the exogenous inputs)
    #   Note: I don't construct a cAB ROM in this example since I know there is no forcing vector
    l2solver = opinf.lstsq.L2Solver(regularizer=1e-7)
    opinf_model = opinf.models.ContinuousModel("AB",solver=l2solver)
    opinf_model.fit(states=uhat, ddts=uhat_ddots,inputs=reduced_stacked_sideset_snapshots)

    # Flip signs to match convention of K on the LHS
    K = -opinf_model.A_.entries
    B = opinf_model.B_.entries

    ## Now extract boundary operators and create dictionary to save
    col_start = 0
    sideset_operators = {}
    for sideset in sidesets:
      num_dofs = reduced_sideset_snapshots[sideset].shape[0]
      val = np.einsum('kr,vnr->vkn',B[:,col_start:col_start + num_dofs] , ss_tspace[sideset].get_basis() )
      shape2 = B[:,col_start:col_start + num_dofs] @ ss_tspace[sideset].get_basis()[0].transpose()
      sideset_operators["B_" + sideset] = val#
      col_start += num_dofs
    
    f = np.zeros(K.shape[0])
    vals_to_save = sideset_operators 
    vals_to_save["basis"] = trial_space.get_basis() 
    vals_to_save["K"] = K
    vals_to_save["f"] = f
    
    np.savez('opinf-operator-holder1',**vals_to_save)




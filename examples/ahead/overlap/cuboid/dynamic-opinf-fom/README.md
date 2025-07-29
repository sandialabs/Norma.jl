#Instructions for running OpInf ROMs in Norma



###Requirements:
- [noma-opinf](https://gitlab-ex.sandia.gov/ejparis/norma-opinf) python package for data processing and opinf ROM construction
- Above package links to UT's (Shane's) [OpInf](https://willcox-research-group.github.io/rom-operator-inference-Python3/source/index.html) package along with [rom-tools-and-workflows ](https://github.com/Pressio/rom-tools-and-workflows)

##OpInf ROM for Schwarz
This example details how to construct an OpInf ROM for Schrwarz of the form

$$\ddot{x} + Kx = Bu + f$$
where $x$ is the state, $u$ are inputs arising, e.g., from BCs, and $K,B$ are the learned matrices. We will work with the cuboid  example.

###1 - Training data generation
First, we are going to update the cuboids.yaml yaml file to write out to CSV and writeout the CSV sidesets. I plan on adding exodus support, but this is not done yet.

```yaml
type: single
input mesh file: ../laser-weld.g
output mesh file: laser-weld.e
Exodus output interval: 1
CSV output interval: 1
CSV write sidesets: true
```
CSV output interval: 1 will result in nodal displacements, accelerations, and velocities being written out to a csv file at each step, while CSV write sidesets: true will write out the the side-set conditions for each side-set listed in the yaml. These are needed, e.g., if you want to incorporate boundary conditions. An additional file will be written for what DOFs are free.

Once this is done, run norma:

```bash
norma cuboids.yaml
```

###2 - Build OpInf ROM
The next step is to build the OpInf ROM. Here, we will build the model for domain 2. The output of this process needs to be a .npz file with keys:

* K: reduced stiffness matrix
* f: reduced forcing vector
* basis: basis for *all* DOFs (free and fixed). Fixed DOFs will be overwritten within Norma
* B\_sideset\_name for all sidesets: sideset matrices where sideset name is:
	*  Standard Dirichlet BCs: The the name of the sideset + x,y,z component (e.g., 'B\_surface\_negative\_y-y').
	*  Schwarz DBCs: The name of the sideset. There is no x,y,z component since Schwarz assigns BCs for all components

This operator can be constructed in any manner, and here we will use the norma-opinf tools. In the following snippet, we will load in the displacement csv files, identify which DOFs are free, and build an OpInf ROM for those DOFs. We will additionally load in the sideset DOFs, which are presecribed via BCs, and treat these as exogenous inputs to the model. We perform dimension reduction on them as well for scalibility. rom-tools-and-workflows is used to perform dimension reduction.

```py
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
    solution_id = 2
    displacement_snapshots,times = normaopinf.readers.load_displacement_csv_files(solution_directory=cur_dir,solution_id=solution_id,skip_files=1)

    # Identify which DOFs are free
    free_dofs = normaopinf.readers.get_free_dofs(solution_directory=cur_dir,solution_id=solution_id)

    # Set values = 0 if DOFs are fixed
    displacement_snapshots[free_dofs[:,:]==False] = 0.

    #Get sideset snapshots
    # Note the different convention for ssz-, which is the Schwarz DBC
    sidesets = ["nsx--x","nsy--y","ssz-","nsz+-z"]
    sideset_snapshots = normaopinf.readers.load_sideset_displacement_csv_files(solution_directory=cur_dir,sidesets=sidesets,solution_id=solution_id,skip_files=1)

    # Create an energy-based truncater
    tolerance = 1.e-5
    my_energy_truncater = romtools.vector_space.utils.EnergyBasedTruncater(1. - tolerance)

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
    u_ddots = normaopinf.calculus.d2dx2(displacement_snapshots,times)
    uhat_ddots = romtools.rom.optimal_l2_projection(u_ddots*1.,trial_space)

    # Construct an opinf "AB" model (linear in the state and linear in the exogenous inputs)
    #   Note: I don't construct a cAB ROM in this example since I know there is no forcing vector
    l2solver = opinf.lstsq.L2Solver(regularizer=5e-3)
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

    np.savez('opinf-operator',**vals_to_save)

```

###3 Run OpInf ROM
Now, we can run an OpInf ROM. We simply need to change the model type to "linear opinf rom" and specify the model-file that we saved above (e.g., "opinf-operator.npz"). The yaml for cuboids-2.yaml should look like this:

```yaml
type: single
input mesh file: cuboid-2.g
output mesh file: cuboid-2.e
model:
  type: linear opinf rom
  model-file: opinf-operator.npz
  material:
    blocks:
      coarse: hyperelastic
    hyperelastic:
      model: neohookean
      elastic modulus: 1.0e+09
      Poisson's ratio: 0.25
      density: 1000.0
time integrator:
  type: Newmark
  β: 0.25
  γ: 0.5
boundary conditions:
  Dirichlet:
    - node set: nsx-
      component: x
      function: "0.0"
    - node set: nsy-
      component: y
      function: "0.0"
    - node set: nsz+
      component: z
      function: "1.0 * t"
  Schwarz overlap:
    - side set: ssz-
      source: cuboid-1.yaml
      source block: fine
      source side set: ssz+
solver:
  type: Hessian minimizer
  step: full Newton
  minimum iterations: 1
  maximum iterations: 16
  relative tolerance: 1.0e-10
  absolute tolerance: 1.0e-06
```


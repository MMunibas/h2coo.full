# h2coo.full


### Data set

The data at CASPT2/aug-cc-pVTZ level of theory is saved in a compressed numpy binary file (H2COO_13877.npz), which contains a python dictionary with seven numpy arrays:

N: Number of atoms in structure (num_data,)
R: Cartesian Coordinates of atoms (in Angstrom [A]), (num_data, N, 3)
Q: Total charge (in elementary charges [e]), (num_data,)
D: Dipole moment vector (in elementary charges times Angstrom [eA]), (num_data, 3)
E: Potential energy with respect to free atoms (in electronvolt [eV]), (num_data,)
F: Forces acting on atoms (in electronvolt per Angstrom [eV/A]), (num_data, N, 3)
Z: Atomic number of atoms, (num_data, N)

The data sets can be accessed using python:
>>> data = np.load("H2COO_13877.npz")

The different keywords of the python dictionary can be listed using
>>> data.files
>>>['N', 'E', 'Q', 'D', 'Z', 'R', 'F']

and the individual entries can be loaded using the appropriate keyword, e.g. for the energy
>>> energies = data["E"]


### PhysNet PES

PhysNet model is stored in best_model.ckpt-245800. Note that, the model was trained using PhysNet based on Tensorflow 1.


### Inputs for MD simualtions

Input for MD simulations with the initial coordinates and velocities are stored in "VibEx_CH_COO". 

The file "test.str" defined the initial excitaions, i.e., the excitation of CH stretch and the COO bending mode.

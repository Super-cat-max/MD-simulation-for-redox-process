# Catalyst Reconstruction Simulation

This program simulates the reconstruction phenomenon of alloy catalysts during oxidation-reduction processes to obtain basic information about active sites on the catalyst surface.

## Features

- Simulates oxidation and reduction steps
- Uses a cyclic calculation strategy under different ensembles
- Generates and analyzes structures at each step
- Calculates adsorption energies
- Performs Molecular Dynamics (MD) simulations

## Main Components

1. `main.py`: The main script that orchestrates the simulation process
2. `Gen_structures.py`: Generates structures for each step
3. `Read_energy.py`: Reads and processes energy calculations
4. `add_O2_first.py`: Handles the initial oxygen addition step
5. `clean_finder.py`: Finds adsorption sites on the surface

## Usage

1. Ensure all dependencies are installed (ASE, NumPy, SciPy, Pymatgen)
2. Set up your initial structure and simulation parameters
3. Run the main script:

```
python main.py
```

## Output

The program generates various files and folders for each oxidation and reduction step, including:

- Structure files (CIF, ARC)
- Energy calculations
- MD simulation results

The final reconstructed structure can be found in the `final_structure.cif` file.

## Note

This simulation is designed for specific alloy catalyst systems. Adjust parameters in the code as needed for your particular system.

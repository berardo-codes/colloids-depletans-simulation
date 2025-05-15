# Monte Carlo Simulation of Hard Sphere Gas (2D)

This project simulates a system of hard colloidal particles and depletants in two dimensions using a Monte Carlo algorithm.

## Features

- Particles are constrained in a periodic box.
- Colloids interact with each other and with depletants via excluded volume.
- Simulation parameters: number of particles, interaction distance, random step size.
- Outputs include particle distance statistics and simulation parameters.

## File Structure

- `montecarlo_hard_spheres_simulation.c` – the main simulation code.
- `Dati Simulazione_v3.dat`, `Simulation##_v3.dat` – output files containing simulation data.
- `output/` – recommended folder for output files.

## Compilation

You can compile the code with GCC:

```bash
gcc -O2 montecarlo_hard_spheres_simulation.c -o mc_sim -lm
```

## Execution

```bash
./mc_sim
```

Output files will be created in the same directory.

## License

This project is licensed under the MIT License. See the `LICENSE` file for more information.

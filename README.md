# GR Compact Objects Toolkit

This project aims to develop a general-purpose toolkit for simulating compact objects in General Relativity. The toolkit is designed to handle a wide range of scenarios, including the construction of compact stars with general metric expressions, accommodating features such as magnetic fields and rotation.

## Project Structure

- `src/`: Core library for GR simulations
  - `lib.rs`: Library root
  - `rk4.rs`: Runge-Kutta 4th order integrator
  - `stellar_structure.rs`: General stellar structure simulation module

- `research_project/`: Application of the core library to specific scenarios
  - `src/`: 
    - `free_fall.rs`: Free fall simulations
    - `main.rs`: Stellar structure simulations of 4 different types of compact objects using the 
    polytropic equation of state.
  - `data/`: Output data from simulations
  - `plots/`: Generated plots and visualizations

## Features

- General-purpose toolkit for GR simulations of compact objects
- Flexible framework to accommodate various metric expressions
- Support for magnetic and rotating star models (planned)
- Currently implements Newtonian simulations as a starting point
- Runge-Kutta 4th order integration method
- Handles both relativistic and non-relativistic equations of state

## Current Capabilities

- Simulates stellar structure for various types of degenerate matter:
  - Non-relativistic and relativistic electron gas
  - Non-relativistic and relativistic neutron gas
- Generates individual stellar profiles
- Produces mass-radius relationships for different central densities

## Installation

1. Run the Newtonian stellar structure simulations in `main.rs`. Under the `main()`, set the variable `gas_type`. After running, this file will generate data files inside `data/`.

```bash
git clone https://github.com/kkin1995/gr-compact-objects.git

cd research_project/

cargo build

cd src/

cargo run
```

2. Generate plots from the simulation data.

```bash
python plot_data.py
```
# WaterSim

## Overview

WaterSim is a small molecular dynamics (MD) program that models a system of rigid water molecules in two dimensions.  The code implements a simple velocity Verlet integration scheme and includes Lennard–Jones and Coulomb interactions.  Molecules are placed on a square lattice and interact with periodic boundary conditions.

The project was originally written for experimentation and is intentionally minimal.  All computation is done in a single process with plain C++ classes.

## Building

The project uses CMake.  A typical build looks like:

```bash
mkdir build
cd build
cmake ..
make
```

This produces the `water_main` executable as defined in `CMakeLists.txt`.

## Usage

After building, run the simulation from the build directory:

```bash
./water_main
```

The default `main` function constructs an `MD` instance with a 19Å box, 225 water molecules, a starting temperature of 3 K and 100 time steps.  The program prints progress information and writes final velocities to `VELout.csv`.

## Core Files

### `Water`

`Water` represents a rigid water molecule and stores its position, velocity and moment of inertia.  Key fields include `m_x` for the centre of mass, `m_v` for translational and rotational velocity and `m_ol`, `m_al` and `m_bl` for the locations of the oxygen and two hydrogens.  Methods allow updating position/velocity and computing derived coordinates.

Relevant declarations appear in `Water.h`:

```cpp
class Water{
    ...
    Force m_f; //The sum of fx, fy, t on the molecule
    Point m_x; //COM location
    Point m_v; //Velocity of molecule
    Point m_ol; //Oxygen location (th = 0)
    Point m_al; //Ha location (th = 0)
    Point m_bl; //Hb location (th = 0)
    double m_I; //Moment of intertia
    void update_pos(Point newpos);
    void update_vel(Point newvel);
    void init();
};
```

### `MD`

`MD` drives the simulation.  It creates the molecule array, performs neighbour searches and advances the system each time step.  `MD::evolve()` runs the main loop and calls `single_timestep()` which implements velocity Verlet integration and periodic boundaries.

Part of the interface in `MD.h`:

```cpp
class MD{
    ...
    int m_N;      //num particles
    double m_L;   //box length angstroms
    double m_T;   //initial temp
    int m_simL;   //number of fs to run the sim
    void evolve();
    ...
};
```

### `Library.h` constants

`Library.h` defines a collection of physical constants and helper structures.  Important constants include:

```cpp
#define HOH_A (109.47*M_PI/180.0) // Degrees
#define ROH 1.0                   // Angstrom
#define EPS 0.650                 // KJ/mol
#define SIGMA 3.166               // Angstrom
#define QO -0.820                 // Charge of Oxygen
#define QH 0.410                  // Charge of Hydrogen
#define TS 1                      // timestep fs
#define NEIGHTS 20                // time between updating neighbors fs
#define RC 9.0                    // cutoff distance in angstroms
```

`Force` and `Point` structures are also declared here and provide basic vector math utilities.

## Notes

* The code assumes the number of molecules forms a perfect square so that molecules can be arranged on a lattice.  Supplying a different value may lead to warnings or incorrect geometry.
* The implementation is two dimensional and meant for experimentation only.  Energy conservation and long term stability have not been extensively tested.
* Planned improvements include exporting trajectory data and supporting three‐dimensional simulations.


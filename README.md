# AMO-TDSE

## What is this?

AMO-TDSE is a Time-Independent and Time-Dependent Schrodinger Equation solver I wrote for my graduate research at JILA. The purpose of the simulation is to model how atoms interact with ultrashort, intense laser pulses and study the phenomena that occurs as a result of that interaction. The simulation is broken into three main steps:

The simulation is broken into three main steps:
1. Solve the TISE to get the bound states for an atom.
2. Solve the TDSE for a laser interacting with the atom, using a bound state as the initial condition.
3. Analyze the final state after the TDSE to study phenomena (Ionization, Excited State Populations etc)

## How does it work?

The simulation relies on PETSc and SLEPc for sparse, parallelized linear algebra in two key places: Solving an eigenvalue problem for the bound states (SLEPc) and solving a linear system to propagate the wavefunction in time (PETSc). This allows the problems to be solved much faster, for example, on a HPC cluster such as the one I used at the university. 

## Numerical Techniques

There are a few numerical techniques that were employed that make this code stand out from other similar projects one can find online.

### BSpline Basis Functions

When dealing with the radial coordinate a natural representation is the finite difference basis. This basis effectively represents your radial coordinate with a discrete set of position eigenstates (usually with a fixed step size). There are a few downsides to this:

1. For atoms whose potential wells are deeper than Hydrogen, the bound states oscillate more quickly near the core. These rapid oscillations make accurate representation more difficult compared to Hydrogen, requiring more grid points in that region in order to converge the correct eigenvalues/eigenvectors.

2. Compounding with the last point, most finite difference schemes assume a fixed spacial step size. The problem with this is that while we may need a 10x smaller step size near the core to converge the eigenvalue problem, we are then forced to use the smaller step size everywhere. This dramatically increases dimensionality, that increased resolution is wasted in regions where it's not required. While where are some finite difference schemes for varying step size, they are more difficult to implement and manage. 

A BSpline basis addresses this in two ways. First, by adjusting the knot vector it becomes trivial to increase the density of BSPlines in a particular region. In fact, you have full control of how these basis functions are distributed. For example,you can have quadratically increasing spacing as you get further from the core, which transitions to a constant linear spacing at some sufficient distance, providing high resolution where needed and not wasting dimensionality where it is not. Second, with a BSpline basis you effectively trade sparsity for dimensionality. This is because a BSpline basis function is able to represent more structure than a position eigenstate (which in position space is a dirac delta). So you need far less of them to accurately represent your wavefunction. 

### Exterior Complex Scaling

There are at least two concessions that must be made when solving the Schrodinger Equation numerically. 





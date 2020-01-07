# MinimumEnergyPoints
Distribute points on hyperspheres such that user-defined energy is minimized.
MATLAB GUI and command line.

Start the GUI by running `Start_ME_Points.m`.
A command line script template is provided in the subfolder `examples/`.

User-defined kernel functions have to be implemented in
`source/PointPotential.m` and in `source/dPointPotential.m`.
See the paper for more information or contact the authors.

This software package is related to the research article

> Oliver Kunc and Felix Fritzen: *Generation of energy-minimizing point sets on
> spheres and their application in mesh-free interpolation and differentiation*,
>
> Advances in Computational Mathematics **45**(5–6), 3021–3056, 2019
>
> DOI   [10.1007/s10444-019-09726-5](https://doi.org/10.1007/s10444-019-09726-5 "Paper in Advances in Computational Mathematics")

For licensing information, see the `LICENSE` file.

# RBFKit

> A Wolfram Language toolkit for radial basis function interpolation and meshless numerical methods.

[![Wolfram Language](https://img.shields.io/badge/Wolfram%20Language-14.2%2B-DD1100)](https://www.wolfram.com/language/)
[![License: MIT](https://img.shields.io/badge/License-MIT-blue.svg)](LICENSE)
![Status](https://img.shields.io/badge/status-early%20development-yellow)

**RBFKit** is a Wolfram Language paclet for experimenting with **radial basis function (RBF)** methods in scientific computing.

The current implementation focuses on RBF interpolation with arbitrary user-defined kernels and optional polynomial augmentation. The package is also intended to provide a foundation for **global RBF collocation methods for partial differential equations (PDEs)**.

## Features

* RBF interpolation in one or more spatial dimensions
* User-defined radial basis functions
* Optional polynomial augmentation
* Callable `RBFInterpolant` objects
* Metadata access for generated interpolants
* Utilities for constructing pairwise RBF/distance matrices
* Paclet-based Wolfram Language project structure
* Global RBF PDE solver infrastructure under development

## Mathematical form

Given centers

$$
X = \{x_1,\ldots,x_N\},
$$

an RBF interpolant has the form

$$
s(x) = \sum_{j=1}^{N} \lambda_j \, \phi\left(\left\lVert x-x_j \right\rVert\right) + p(x).
$$

where:

- $\phi(r)$ is a radial basis function,
- $\lambda_j$ are interpolation coefficients,
- $p(x)$ is an optional polynomial term.

RBFKit constructs and solves the corresponding augmented interpolation system automatically.

## Requirements

* **Wolfram Language 14.2 or newer**
* Mathematica or another environment capable of evaluating Wolfram Language paclets

## Installation

Clone the repository:

```bash
git clone https://github.com/LDBreton/RBFKit.git
```

The paclet itself is located in the inner `RBFKit` directory. From Wolfram Language, install the local checkout with:

```wl
PacletInstall["/absolute/path/to/RBFKit/RBFKit"]
```

Then load the interpolation package:

```wl
Needs["RBFKit`RBFInterpolation`"]
```

During development, you can reinstall the paclet after making changes:

```wl
PacletUninstall["RBFKit"];
PacletInstall["/absolute/path/to/RBFKit/RBFKit"];
```

## Quick start

### 1D interpolation

Create sample data:

```wl
data = Table[
    {x, Sin[2 Pi x]},
    {x, 0., 1., 0.1}
];
```

Define any radial kernel as a function of the radial distance `r`. For example, a Gaussian RBF:

```wl
gaussian[r_] := Exp[-(3 r)^2];
```

Construct the interpolant:

```wl
interp = RBFInterpolation[data, gaussian];
```

The returned object can be evaluated directly:

```wl
interp[0.35]
```

Compare the interpolation with the original function:

```wl
Plot[
    {
        Sin[2 Pi x],
        interp[x]
    },
    {x, 0, 1},
    PlotLegends -> {"Exact", "RBF interpolation"}
]
```

## Polynomial augmentation

Polynomial terms can be added with the `"PolynomialDegree"` option:

```wl
interp = RBFInterpolation[
    data,
    gaussian,
    "PolynomialDegree" -> 1
];
```

The default is:

```wl
"PolynomialDegree" -> 0
```

which uses only the RBF part of the interpolant.

## Multidimensional interpolation

Input data are stored row-wise in the form

```text
{x1, x2, ..., xd, f}
```

where the first `d` entries are coordinates and the last entry is the function value.

For example:

```wl
data2D = Flatten[
    Table[
        {x, y, Sin[Pi x] Cos[Pi y]},
        {x, 0., 1., 0.2},
        {y, 0., 1., 0.2}
    ],
    1
];

imq[r_] := 1/Sqrt[1 + (2 r)^2];

interp2D = RBFInterpolation[data2D, imq];

interp2D[0.25, 0.4]
```

## Interpolant metadata

`RBFInterpolation` returns an `RBFInterpolant` object. In addition to being callable, it stores useful information about the interpolation:

```wl
interp["Variables"]
interp["Degree"]
interp["NCenters"]
```

For example:

```wl
interp["NCenters"]
```

returns the number of RBF centers used to construct the interpolant.

## Utility functions

RBFKit also provides a utility for constructing matrices from two sets of points.

Load the utilities package:

```wl
Needs["RBFKit`RBFutils`"]
```

Then use:

```wl
distancematrix[f, nodes, centers]
```

where `f` is a two-argument function applied to every pair of points.

For a Euclidean distance matrix:

```wl
D = distancematrix[
    Norm[#1 - #2] &,
    nodes,
    centers
];
```

For an RBF matrix:

```wl
A = distancematrix[
    gaussian[Norm[#1 - #2]] &,
    nodes,
    centers
];
```

## API overview

| Symbol                                                 | Description                                                     |
| ------------------------------------------------------ | --------------------------------------------------------------- |
| `RBFInterpolation[data, rbf]`                          | Construct an RBF interpolant from scattered data                |
| `RBFInterpolation[data, rbf, "PolynomialDegree" -> p]` | Construct an RBF interpolant with polynomial augmentation       |
| `RBFInterpolant[...]`                                  | Object returned by `RBFInterpolation`; callable like a function |
| `distancematrix[f, nodes, centers]`                    | Apply a two-point function to all node-center pairs             |

## Project structure

```text
RBFKit/
├── LICENSE
├── README.md
└── RBFKit/
    ├── PacletInfo.wl
    └── Kernel/
        ├── init.m
        ├── RBFInterpolation.wl
        ├── RBFutils.wl
        └── RBFGlobalSolve.wl
```

The package is intentionally split into modules so interpolation, utilities, and PDE-solving functionality can evolve independently.

## Development status

RBFKit is currently under active development.

### Implemented

* RBF interpolation
* Custom radial kernels
* Polynomial augmentation
* `RBFInterpolant` wrapper
* Pairwise matrix construction utilities

### In development

* Global RBF collocation for PDEs
* Differential-operator assembly
* Boundary-condition handling
* Higher-level PDE solver interface
* Expanded examples and documentation

The API may change as these components are developed.

## Planned PDE workflow

The long-term global RBF solver is intended to support workflows of the form

```text
Geometry / nodes
       ↓
RBF basis
       ↓
Differential and boundary operators
       ↓
Global collocation matrix
       ↓
Linear or nonlinear algebraic system
       ↓
RBF approximation of the PDE solution
```

This is aimed at meshless discretizations of problems such as Poisson, diffusion, advection-diffusion, and related PDEs.

## Contributing

Contributions, bug reports, numerical examples, and suggestions are welcome.

If you find an issue, please open a GitHub issue and include:

* Wolfram Language version
* minimal reproducible example
* RBF/kernel used
* node set or data set
* expected and observed behavior

## Author

**Louis Breton**

RBFKit is developed for research and experimentation with radial basis functions, meshless approximation, and numerical PDE methods in the Wolfram Language.

## License

RBFKit is released under the [MIT License](LICENSE).

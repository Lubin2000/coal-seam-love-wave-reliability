# QEDispInv Windows Build and Validation

## Source version

- Repository: https://github.com/pan3rock/QEDispInv
- Commit: `ed1b5dd7`
- License at repository root: LGPL-3.0
- Source headers in the inspected copy: GPL-3.0

Pinned submodules:

- fmt: `40626af88bd7df9a5fb80be7b25ac85b122d6c21`
- toml11: `be08ba2be2a964edcdb3d3e3ea8d100abc26f286`
- HighFive: `1a5ef376f6ceb5ca111a02c936f9918699a7a0b0`
- LBFGSpp: `c524a407fb85b74807f53de5a3ca2ddbcc164e54`

## Windows build environment

A portable UCRT64 MSYS2 environment was used outside the system path.

- MSYS2 base archive SHA-256:
  `421ee2cafd704b4c27d1c8ccb40febbb6334732284ef1c546d8454eaeda7cd3c`
- GCC: 16.2.0
- gfortran: 16.2.0
- CMake: 4.4.5
- Ninja: 1.13.2
- Eigen: 5.0.1
- HDF5: 2.2.0
- OpenMP: disabled for the initial Windows build

Build command:

```bash
cmake -G Ninja -S . -B build-win -DUSE_OPENMP=OFF -DCMAKE_BUILD_TYPE=Release
cmake --build build-win --target forward -j 2
```

## Local compatibility patches

These changes are limited to the build copy:

1. Removed QEDispInv's strict `Eigen3 3.3` version check because the current
   package exports Eigen 5.0.1 as a newer ABI-compatible header package.
2. Removed the same strict version check in the pinned LBFGSpp submodule.
3. Added fallback propagation of `Eigen3::Eigen` include directories because
   Eigen 5 exports the target without defining the legacy
   `EIGEN3_INCLUDE_DIR` variable.
4. Added `#undef FILE_CREATE` in HighFive's `H5PropertyList.hpp` because
   Windows headers define `FILE_CREATE` as a macro and collide with the
   HighFive enum of the same name.

No dispersion or Fortran numerical code was changed.

## Coal-seam comparison model

- Coal thickness: 2.4 m.
- Coal: \(V_S=1.12\) km/s, density 1.38 g/cm3, \(V_P=2.0\) km/s.
- Roof: \(V_S=2.4\) km/s, density 2.45 g/cm3, \(V_P=4.0\) km/s.
- Floor: \(V_S=2.65\) km/s, density 2.6 g/cm3, \(V_P=4.3\) km/s.
- Roof depth interval: 0-100 m.
- Coal depth interval: 100-102.4 m.
- QEDispInv uses a free surface at 0 km and a floor half-space below 102.4 m.

The 100-m roof is an approximation to the analytical upper half-space.
QEDispInv also returns free-surface roof modes. Each analytical coal mode is
therefore matched to its nearest root among all QEDispInv modes.

Command:

```bash
../../bin/forward.exe --sh -m 100 -o qedisp_love_all.txt
```

## Results

- Frequency range: 100-560 Hz at 20-Hz spacing.
- Coal-mode roots compared: 34.
- Within 1 m/s: 27/34.
- Within 10 m/s: 34/34.
- Median absolute difference: 0.477 m/s.
- Mean absolute difference: 0.972 m/s.
- P95 absolute difference: 2.582 m/s.
- Maximum absolute difference: 9.104 m/s.

The largest discrepancy is at 100 Hz, where the fundamental analytical branch
approaches the roof velocity and is most sensitive to replacing the roof
half-space by a finite layer.

Compiled executable SHA-256:

`AC142210A190A3EE8DD467CC66E941B351431479240DC9E7C19AEABBC4394750`

Raw QEDispInv output SHA-256:

`0071861A7DC39FD41FF46755047073CD88395051BA07AC2198EF0D98C2094C35`

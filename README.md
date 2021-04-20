# Cobiveco
**Co**nsistent **bive**ntricular **co**ordinates (https://arxiv.org/abs/2102.02898)

![Coordinates](https://user-images.githubusercontent.com/31965820/103575830-8bc56400-4ed2-11eb-913d-1bc1c338622d.png)

## Overview

This is a MATLAB implementation to compute local coordinates on tetrahedral meshes of biventricular cardiac geometries. Also provided are functions to utilize the results for transferring data, standardized visualization of data and alignment of the heart with the global coordinate axes.

## Dependencies

Before using Cobiveco, its dependencies [[1](#1),[2](#2),[3](#3)] need to be installed by running [dependencies/install_cobiveco.sh](dependencies/install_cobiveco.sh).

## Computing coordinates

Coordinates can be computed with:

```matlab
% Create a cobiveco object, providing a config struct with input and output path prefixes
c = cobiveco(struct('inPrefix','geo1/heart', 'outPrefix','result_geo1/'));

% Run computation of all coordinates
c.computeAll;

% Optional: Retrieve the result and a config struct with all parameters
result = c.result;
config = c.cfg;
```

Run ``help cobiveco`` for more information.

## Utilities

The following [utilities](utilities) are provided:
* [cobiveco_computeMappingMatrix.m](utilities/cobiveco_computeMappingMatrix.m): Computes a matrix to map point data from a source to a target mesh.
* [cobiveco_createPolarProjection.m](utilities/cobiveco_createPolarProjection.m): Creates a standardized visualization by projecting scalar data onto polar plots.
* [cobiveco_applyAlignmentMatrix.m](utilities/cobiveco_applyAlignmentMatrix.m): Uses the matrix ``R`` computed by Cobiveco to align the heart axes with the global coordinate axes.

## Examples

The following [examples](examples) illustrate the use of Cobiveco:
* [example_geo1.m](examples/example_geo1.m): Computes coordinates on [geo1](examples/geo1) &ndash; the mean shape of a statistical shape model [[4](#4),[5](#5)].
* [example_geo2.m](examples/example_geo2.m): Computes coordinates on [geo2](examples/geo2) &ndash; a patient geometry.
* [example_mapping.m](examples/example_mapping.m): Uses the coordinates to map data between geo1 and geo2.
* [example_polarProjection.m](examples/example_polarProjection.m): Uses the coordinates to project data onto polar plots.

## Hints

To speed up the computation for fine meshes (several millions of nodes), we recommend to increase the sizing parameters of Mmg. This way, all coordinates except the binary transventricular coordinate are automatically computed on coarser meshes and interpolated to the original mesh. As the non-binary coordinates are spatially low-frequent, an edge length of slightly below 1 mm (or at maximum one third of the smallest wall thickness) is sufficient.
If you want to compute coordinates for a mesh with an average edge length of 0.3 mm, for example, scale the default value for the input 'mmgSizingParam' by a factor of 0.9/0.3 = 3 to use an effective edge length of about 0.9 mm for the computations:

```matlab
c = cobiveco(struct('inPrefix','heart', 'outPrefix','result/', 'mmgSizingParam',3*[0.1 0.9 1.1]));
```

## License

All source code is subject to the terms of the Apache License 2.0.  
Copyright 2021 Steffen Schuler, Karlsruhe Institute of Technology.

## Contact

Steffen Schuler  
Institute of Biomedical Engineering  
Karlsruhe Institute of Technology  
www.ibt.kit.edu

## References

<a id="1">[1]</a> [Dapogny, C. et al., 2014. Three-dimensional adaptive domain remeshing, implicit domain meshing, and applications to free and moving boundary problems. J Comput Phys.](https://github.com/MmgTools/mmg)  
<a id="2">[2]</a> [Jacobson, A., 2018. gptoolbox: Geometry processing toolbox.](https://github.com/alecjacobson/gptoolbox)  
<a id="3">[3]</a> [Schuler, S., 2020. vtkToolbox: A MEX interface to the VTK library.](https://github.com/KIT-IBT/vtkToolbox)  
<a id="4">[4]</a> [Bai, W. et al., 2015. A bi-ventricular cardiac atlas built from 1000+ high resolution MR images of healthy subjects and an analysis of shape and motion. Med Image Anal 26, 133–45.](https://github.com/UK-Digital-Heart-Project/Statistical-Shape-Model)  
<a id="5">[5]</a> [Schuler, S. et al., 2021. Biventricular statistical shape model of the human heart adapted for computer simulations. Zenodo.](https://doi.org/10.5281/zenodo.4419783)

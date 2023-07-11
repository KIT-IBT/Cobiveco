# CobivecoX

**Co**nsistent **bive**ntricular **co**ordinates E**X**tended (https://arxiv.org/abs/2102.02898, submitted)

### CobivecoX on geometries including the four valve annuli
![abstract_medima2023_whitebackground](https://user-images.githubusercontent.com/51908398/235209955-ebe255e9-03d2-462b-a67d-7cf24192612a.png)

### Cobiveco 1.0 on geometries with a flat base
![Coordinates](https://user-images.githubusercontent.com/31965820/103575830-8bc56400-4ed2-11eb-913d-1bc1c338622d.png)

## Overview

This is a MATLAB implementation to compute local coordinates on tetrahedral meshes of biventricular cardiac geometries. Also provided are functions to utilize the results for transferring data, standardized visualization of data and alignment of the heart with the global coordinate axes.
CobivecoX is backwards compatible with Cobiveco and can be applied to both meshes with and without including the valve annuli.

## How to install CobivecoX (macOS Monterey)

1. A MATLAB installation version 2022a (or above) is required to run CobivecoX. Check [this link](https://se.mathworks.com/help/install/install-products.html) on how to install it if required.

Here is the list of all Add-ons needed to add in MATLAB:

* Image Processing Toolbox
* Mapping Toolbox
* Signal Processing Toolbox
* Statistics and Machine Learning Toolbox
* Curve Fitting Toolbox
* Parallel Computing Toolbox

If you already have MATLAB and need to change licences, see [this link](https://se.mathworks.com/matlabcentral/answers/99457-how-do-i-activate-matlab-or-other-mathworks-products) to activate the app. Ones MATLAB is running, you can check if a spesific add-on is anabled by using the command `matlab.addons.isAddonEnabled('NameOfToolbox')` (it returns the logical value 1 if it is enabled).
Alternatively, `addons = matlab.addons.installedAddons` shows a list of all installed add-ons.

2. Install vtk (For MacOS you can follow [this guide](https://stackoverflow.com/questions/32853082/how-to-install-vtk-on-a-mac)).
3. Install cmake using the command

   ```
   $ brew install cmake 
   ```
4. Install the dependencies [[1](#1),[2](#2),[3](#3)] needed running the [dependencies/install_cobiveco.sh](dependencies/install_cobiveco.sh) using the command

   ```
   $ sh install_cobiveco.sh 
   ```

   **Note**: once mmg is installed, check that the ``bin`` directory is non-empty using the command

   ```
   $ ls dependencies/mmg/build/bin
   ```

   If the files ``mmg2d_O3, mmg3d_O3, mmgs_O3`` are present, you can go on to the next step. Else perform download the mmg binary from [here](https://www.mmgtools.org/mmg-remesher-downloads) and move ``mmg2d_O3, mmg3d_O3, mmgs_O3`` to ``dependencies/mmg/build/bin``.

   When running any example later, a warning about unidentified developer may occur. This can be solved by editing the setting as descibed [here](https://support.apple.com/en-us/HT202491).
5. Paraview will be used during the calculation of the coordinates and is highly recommended for visualization. Install Paraview. You can download Paraview [here](https://www.paraview.org/download/) (recommended on Mac: ParaView-5.11.0-RC1-MPI-OSX10.13-Python3.9-x86_64.dmg). The installation was tested was Paraview 5.10 and 5.11 .
6. Check that ``pvpython`` can be evoked from the commandline. If not, add the path to your ``~/.zshrc``.
   First, check where ``pvpython`` is located.

   ```$
   $ which pvpython
   ```

   The output will look something like this:

   ```
   /Applications/ParaView-5.10.1.app/Contents/bin/pvpython
   ```

   Modify the following commands to represent your path. Then, add pvpython to your `` ~/.zshrc`` using your favorite editor or just by doing:

   ```
   $ echo "export PATH="$PATH:/Applications/ParaView-5.10.1.app/Contents/bin/"" >> ~/.zshrc
   $ ln -s /Applications/ParaView-5.10.0.app/Contents/bin/pvpython ~/bin
   ```
7. To switch between environments in Matlab, we recommend using conda. However, one can also create a python environment using the ``requirements.txt`` file, which is more light-weight, than choosing anaconda. However, one should be aware, that you can switch between calling ``pvpython`` and ``python`` and might need to modify the execution of condalab in the matlab code. See step 9 below.

To create a python environment one can simply execute the following:
   ```
   pip install -r requirements.txt
   ```
or use conda to create a minimal environment
   ```
   conda install --file requirements.txt
   ```

The tested and recommended version of this code is currently using anaconda. You can install Anaconda using [this guide](https://docs.anaconda.com/anaconda/ install/mac-os/) (commandline-installation recommended).
   This [guide for conda](https://conda.io/projects/conda/en/latest/user-guide/getting-started.html) might also be useful if the is user not familiar with it.

8. Create a conda environment called base using the  ``cobivecox.yml. ``

   You can use the following command to create the python environment with the dependencies needed:

   ```
   $ conda env create -f cobivecox.yml
   ```
9. The conda path has to be put in manually in the code.Use the command

   ```
   $ which conda
   ```

   or

   ```
   $ conda env list
   ```

   in the terminal to receive the path.

   Copy the path you got in line 19 of the function ``createBoundarySurfacesAB.m`` (you can find this function in the directory ``functions``)

   ```
   $ condaPath = '{*the path from which conda*}';
   ```
10. To avoid communication errors, open MATLAB from your terminal.
    One can add the shortcut to their ``~/.zshrc``

```$
   $ alias matlab="/Applications/MATLAB_R2022b.app/bin/matlab"
```

   and then use the following command in the directory ``examples``

```$
   $ matlab
```

   to open the program.

## Dependencies

Before using Cobiveco, its dependencies [[1](#1),[2](#2),[3](#3)] need to be installed by running [dependencies/install_cobiveco.sh](dependencies/install_cobiveco.sh).
Note that this script was tested on MacOS.

## Usage

1. Open MATLAB from the commandline.

If you set the ``matlab``alias as described above you can do that in the directory ``examples`` using the follwoing command:

```$
$ matlab
```

to open the program.

## Input Files To Compute Coordinates

1. To use CobivecoX, the user needs to provide the following input files in a directory named `{filename}` containing the following files:

* {filename}_av.ply
* {filename}_endo_lv.ply
* {filename}_endo_rv.ply
* {filename}_epi_base.ply
* {filename}_epi_no_base.ply
* {filename}_mv.ply
* {filename}_pv.ply
* {filename}_tv.ply
* {filename}.vtu

Note that the folder `{filename}` needs to be in `examples`.

2. In case of needing the original Cobiveco, the following input files have to be provided:

* {filename}_base.ply
* {filename}_endo_lv.ply
* {filename}_endo_rv.ply
* {filename}_epi.ply
* {filename}.vtu

## Computing coordinates

Coordinates can be computed with:

```matlab
% Create a cobiveco object, providing a config struct with input and output path prefixes
c = cobiveco(struct('inPrefix','{filename}/{filename}', 'outPrefix','result_{filename}/'));

% Run computation of all coordinates
c.prepareMesh0;
if c.cfg.CobivecoX == true
    c.computeAllCobivecoX;
else 
    c.computeAllCobiveco;
end

% Optional: Retrieve the result and a config struct with all parameters
result = c.result;
config = c.cfg;
```

Run ``help cobiveco`` for more information.

An example can be found in the directory `examples`  named ``examples_single_cobiveco_run.m.``
*Please note that at line 18 all `{filename}` placeholders have to be substituted with the name of the file and directory names.*

Finally, once MATLAB has finished exporting the result, it can be visualized opening ParaView and selecting the ``result.vtu`` file in the ``result_{filename}`` directory.

## Utilities

The following [utilities](utilities) are provided:

* [cobiveco_computeMappingMatrix.m](utilities/cobiveco_computeMappingMatrix.m): Computes a matrix to map point data from a source to a target mesh.
* [cobiveco_createPolarProjection.m](utilities/cobiveco_createPolarProjection.m): Creates a standardized visualization by projecting scalar data onto polar plots.
* [cobiveco_applyAlignmentMatrix.m](utilities/cobiveco_applyAlignmentMatrix.m): Uses the matrix ``R`` computed by Cobiveco to align the heart axes with the global coordinate axes.

## Examples

The following [examples](examples) illustrate the use of Cobiveco:

* [example_geo1.m](examples/example_geo1.m): Computes coordinates on [geo1](examples/geo1) &ndash; the mean shape of a statistical shape model [[4](#4),[5](#5)].
* [example_geo2.m](examples/example_geo2.m): Computes coordinates on [geo2](examples/geo2) &ndash; a patient geometry.
* [example_mapping_Cobiveco.m](examples/example_mapping_Cobiveco.m): Uses the coordinates to map data between geo1 and geo2 defaulting to Cobiveco1.0.
* [example_polarProjection.m](examples/example_polarProjection.m): Uses the coordinates to project data onto polar plots.
* [example_mapping_batch_CobivecoX.m](examples/example_mapping_batch_CobivecoX.m): calculates two way mapping error as described by [Bayer et al., 2018](https://www.sciencedirect.com/science/article/pii/S1361841518300203?via%3Dihub) by mapping coordinates between reference and cohort.
* [example_CobivecoX.m](examples/example_CobivecoX.m): Computes coordinates on the mean shape of a Tetralogy of Fallot atlas including valve annuli.
* [example_parcellation.m](examples/example_parcellation.m) and [example_extended_parcellation.m](examples/example_extended_parcellation.m): Add a parcellation in four different areas (left anterior, left posterior, right anterior and right posterior) to the given result form cobiveco or cobivecoX respectively.
* [example_aha_parcellation.m](examples/example_aha_parcellation.m) and [example_extended_aha_parcellation.m](examples/example_extended_aha_parcellation.m): Add a parcellation according to the [AHA 17 segmentation model](https://www.ahajournals.org/doi/10.1161/hc0402.102975) and extending it also to the right ventricle to the given result form cobiveco or cobivecoX respectively.

## Hints

To speed up the computation for fine meshes (several millions of nodes), we recommend to increase the sizing parameters of Mmg. This way, all coordinates except the binary transventricular coordinate are automatically computed on coarser meshes and interpolated to the original mesh. As the non-binary coordinates are spatially low-frequent, an edge length of slightly below 1 mm (or at maximum one third of the smallest wall thickness) is sufficient.
If you want to compute coordinates for a mesh with an average edge length of 0.3 mm, for example, scale the default value for the input 'mmgSizingParam' by a factor of 0.9/0.3 = 3 to use an effective edge length of about 0.9 mm for the computations:

```matlab
c = cobiveco(struct('inPrefix','heart', 'outPrefix','result/', 'mmgSizingParam',3*[0.1 0.9 1.1]));
```

## License

All source code is subject to the terms of the Apache License 2.0.
Copyright 2021 Steffen Schuler, Karlsruhe Institute of Technology.

## Contacts

Lisa Pankewitz
Simula Research Laboratory
www.simula.no

Steffen Schuler & Axel Loewe
Institute of Biomedical Engineering
Karlsruhe Institute of Technology
www.ibt.kit.edu

## References

[1] [Dapogny, C. et al., 2014. Three-dimensional adaptive domain remeshing, implicit domain meshing, and applications to free and moving boundary problems. J Comput Phys.](https://github.com/MmgTools/mmg)

[2] [Jacobson, A., 2018. gptoolbox: Geometry processing toolbox.](https://github.com/alecjacobson/gptoolbox)

[3] [Schuler, S., 2020. vtkToolbox: A MEX interface to the VTK library.](https://github.com/KIT-IBT/vtkToolbox)

[4] [Bai, W. et al., 2015. A bi-ventricular cardiac atlas built from 1000+ high resolution MR images of healthy subjects and an analysis of shape and motion. Med Image Anal 26, 133–45.](https://github.com/UK-Digital-Heart-Project/Statistical-Shape-Model)

[5] [Schuler, S. et al., 2021. Biventricular statistical shape model of the human heart adapted for computer simulations. Zenodo.](https://doi.org/10.5281/zenodo.4419783)

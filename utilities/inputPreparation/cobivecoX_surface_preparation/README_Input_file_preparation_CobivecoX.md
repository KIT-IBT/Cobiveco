# README_Input_file_preparation_CobivecoX.md

Example directory structure befor running the input preparation for cobivecoX:

```
.
├── README_Input_file_preparation_CobivecoX.md
├── data
│   ├── README.md
│   ├── input
│   │   ├── filename.msh
│   │   ├── filename_vol.msh
│   │   ├── README.md
│   │   └── example.xlsx
│   └── output
│       └── README.md
└── scripts
    ├── cp2cobiveco.py
    ├── extract.py
    ├── extract_midmyocard.py
    ├── extractfunc.py
    ├── extractfunc_orig.py
    ├── msh2ply.py
    ├── msh2surfacesply.py
    └── prepare_input_surface.py
```

Directory structure after running the input preparation for cobivecoX:

```
.
├── README_Input_file_preparation_CobivecoX.md
├── data
│   ├── README.md
│   ├── input
│   │   ├── filename.msh
│   │   ├── filename_vol.msh
│   │   ├── README.md
│   │   └── example.xlsx
│   └── output
│       ├── filename
│       │   ├── filename.ply
|       |   ├── filename.config
│       │   └── filename.vtu
│       ├── README.md
│       ├── ids_created.csv
│       └── ids_failed.csv
└── scripts
    ├── cp2cobiveco.py
    ├── extract.py
    ├── extract_midmyocard.py
    ├── extractfunc.py
    ├── extractfunc_orig.py
    ├── msh2ply.py
    ├── msh2surfacesply.py
    └── prepare_input_surface.py
```

Please note that the scripts work for gmsh files (in ``input``: ``filename.msh`` and ``filename_vol.msh``) version 2.2 and in ascii format.

Make sure that the directory structure is illustrated as in the example above.
In the ``example.xlxs``, write down all your file names in the first column ('Instance ID'). You might need to change the header name if needed.
Run the script ``prepare_input_surfaces.py`` in the ``scripts`` directory:

```
$ python prepare_input_surfaces.py
```

**Warning:** The start indices in the created ``.config`` file have to be input manually.
If you get the message:

The configuration file does not contain the required start indices, please modify it according to the``README.md file``.


Follow this step by step guide:

1. Open the newly created ``{filename}.ply``file in ParaView.
2. Select the option "Hover Points On."
3. Go with your pointer on any point of the mitral valve ('mv') and mark down the "Id" number in ``{filename}.config`` after the '=' sign. Preferentially, select points in the middle of a surface!
4. Repeat for every surface in the file. Please be sure to take points in the endocardium for both the right and the left ventricle.
5. Save the file and run the script again.

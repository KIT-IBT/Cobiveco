#!/bin/bash

echo -e "\n========== Installing mmg ==========\n"
if [[ -d "mmg" ]]
then
	echo "The directory 'mmg' already exists. Remove it to reinstall."
else
	git clone https://github.com/MmgTools/mmg.git
	mkdir mmg/build
	cd mmg/build
	cmake ..
	make -j
	cd ../..
fi

echo -e "\n========== Installing gptoolbox ==========\n"
if [[ -d "gptoolbox" ]]
then
	echo "The directory 'gptoolbox' already exists. Remove it to reinstall."
else
	git clone https://github.com/alecjacobson/gptoolbox.git
fi

echo -e "\n========== Installing vtkToolbox ==========\n"
if [[ -d "vtkToolbox" ]]
then
	echo "The directory 'vtkToolbox' already exists. Remove it to reinstall."
else
	git clone https://github.com/KIT-IBT/vtkToolbox.git
	mkdir vtkToolbox/build
	cd vtkToolbox/build
	cmake .. # you might want to specify the path to VTK here, e.g.: cmake .. -D VTK_DIR=../../vtk-8.2.0/build
	make -j
	cd ../..
fi

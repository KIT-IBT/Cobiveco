#!/bin/bash

echo -e "\n========== Installing numerical-tours ==========\n"
if [[ -d "numerical-tours" ]]
then
	echo "The directory 'numerical-tours' already exists. Remove it to reinstall."
else
	git clone --filter=blob:none --no-checkout https://github.com/gpeyre/numerical-tours.git
	cd numerical-tours
	git config core.sparseCheckout true
	echo "matlab/toolbox_general/" >> .git/info/sparse-checkout
	echo "matlab/toolbox_graph/" >> .git/info/sparse-checkout
	git checkout master
	cd ..
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

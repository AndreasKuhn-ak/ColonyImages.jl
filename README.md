# ColonyImages

ColonyImages is a local Julia package for quantification of Trypanosoma colony images. It provides a set of functions to manipulate and analyze colony images. Additionally, it allows for the creation of artificial colony images that can be used as test data.

[![Stable](https://img.shields.io/badge/docs-stable-blue.svg)](https://andreaskuhn-ak.github.io/ColonyImages.jl/)
[![Build Status](https://github.com/AndreasKuhn-ak/ColonyImages.jl/actions/workflows/CI.yml/badge.svg?branch=master)](https://github.com/AndreasKuhn-ak/ColonyImages.jl/actions/workflows/CI.yml?query=branch%3Amaster)
[![Coverage](https://codecov.io/gh/AndreasKuhn-ak/ColonyImages.jl/branch/master/graph/badge.svg)](https://codecov.io/gh/AndreasKuhn-ak/ColonyImages.jl)

## How to Use

To get started with ColonyImages.jl, follow these steps:

1. **Download the Repository**: Clone or download the ColonyImages.jl repository from GitHub to your local machine.

   ```bash
   git clone https://github.com/AndreasKuhn-ak/ColonyImages.jl.git  
    ``` 
2. **Navigate to the Package Directory**:  Open a terminal and change directory to the ColonyImages.jl folder.
   ```bash
   cd ColonyImages
    ``` 
3. **Activate the Julia Environment**: Start Julia in the terminal and activate the package environment.
   ```bash
    julia  
      ``` 
Within the Julia REPL, activate and instantiate the project:
```julia
using Pkg
Pkg.activate(".")
Pkg.instantiate()
```
This sets up the environment with all necessary dependencies.

3. **Using ColonyImages.jl**: Now, you can start using ColonyImages.jl in your Julia scripts or REPL.
```julia
using ColonyImages
```
You're now ready to use the functions provided by ColonyImages.jl for analyzing and manipulating colony images.
## Features

### Image Functions

The `image_functions.jl` file contains a variety of functions for processing and analyzing stacked tif images of Trypanosoma colonies. These functions provide the tools necessary to quantify and understand the growth patterns of these colonies.

One of the key functions in this file is `angular_metric`, which is inspired by the paper "Quantifying Two-Dimensional Filamentous and Invasive Growth Spatial Patterns in Yeast Colonies" by Binder et al. This function measures the angular distribution of a colony's growth, providing a quantitative measure of the colony's spatial pattern. It calculates an angular metric for a given image, where each element in the resulting vector represents the number of pixels in a certain angular sector of the image.

Here is an example of how to use this function:

```julia
using ColonyImages

# Load a single image
img = load("path_to_your_image.tif")
# Make it binary
img_int = b_w(img)

# Define the center of the circle
center = centroid(img)

# Calculate the angular metric
angular_metric = angular_metric(img, center)
```

This will return a vector where each element represents the number of pixels in a certain angular sector of the image. You can then use this vector to analyze the spatial distribution of the colony in the image.



### Artificial Colony Creation

The `artificial_colony_creation.jl` file contains functions for creating artificial colony images. These images can be used as test data for the image functions.



### Tutorial 

To assist users in understanding how to utilize these functions, we have provided a Jupyter notebook titled `Tutorial.ipynb`. This notebook guides you through the process of analyzing test data using the functions in `image_functions.jl` and simulating colony growth using the functions from `artificial_colony_creation.jl`.

It offers a practical example of our workflow to analyze and compare preprocessed colony images with simulated colonies.

The image preprocessing of fluorescent colony images is automated in `Fiji` with a custom-written macro. We have also provided a tutorial (`Fiji_tutorial.pdf`) on how to use it, along with all necessary files in the `fiji` folder.

### Documentation
You are not confined to using `ColonyImages.jl` in a specific way. All functions are documented and can be used independently or in any combination of your choice. To access the documentation, use the `Julia` macros `@doc` or our custom macro `@h` for a more nicely formatted output of the docstrings. Alternatively, you can find all the [documentation](https://andreaskuhn-ak.github.io/ColonyImages.jl/) on our GitHub pages.

### Further Content 
We have also provided all Jupyter notebooks (image import, colony simulation, analysis) together with their outputs, that were used to create our paper `Quantification of Trypanosoma brucei social motility indicates different colony growth phases`. These notebooks can be used if you wish to reproduce our results. However, they are not documented, so we highly recommend working through the `tutorial.ipynb` first to gain an understanding of the workflow and the package.

### Contributing
Contributions to ColonyImages are very welcome! If you have a feature request, bug report, or proposal, please open an issue on the GitHub repository.

### License
ColonyImages is licensed under the MIT license.

### Software Requirements
All software & Julia packages and their respective versions are listed in the `requirements.md` file.
# ColonyImages

ColonyImages is a Julia package for quantification of Trypansoma colony images. It provides a set of functions to manipulate and analyze colony images. Additionally, it allows for the creation of artificial colony images that can be used as test data.

[![Stable](https://img.shields.io/badge/docs-stable-blue.svg)](https://andreaskuhn-ak.github.io/ColonyImages.jl/)
[![Build Status](https://github.com/AndreasKuhn-ak/ColonyImages.jl/actions/workflows/CI.yml/badge.svg?branch=master)](https://github.com/AndreasKuhn-ak/ColonyImages.jl/actions/workflows/CI.yml?query=branch%3Amaster)
[![Coverage](https://codecov.io/gh/AndreasKuhn-ak/ColonyImages.jl/branch/master/graph/badge.svg)](https://codecov.io/gh/AndreasKuhn-ak/ColonyImages.jl)

## Features

### Image Functions

The `image_functions.jl` file contains a variety of functions for processing and analyzing stacked tif images of Trypanosoma colonies. These functions provide the tools necessary to quantify and understand the growth patterns of these colonies.

One of the key functions in this file is `angular_metric`, which is inspired by the paper "Quantifying Two-Dimensional Filamentous and Invasive Growth Spatial Patterns in Yeast Colonies" by Binder et al. This function measures the angular distribution of a colony's growth, providing a quantitative measure of the colony's spatial pattern. It calculates an angular metric for a given image, where each element in the resulting vector represents the number of pixels in a certain angular sector of the image.

Here is an example of how to use this function:

```julia
using ColonyImages

# Load a single image
img = load("path_to_your_image.tif")
#make it binary
img_int = b_w(img)


# Define the center of the circle
center = centroid(img)

# Calculate the angular metric
angular_metric = angular_metric(img, center)
```

This will return a vector where each element represents the number of pixels in a certain angular sector of the image. You can then use this vector to analyze the spatial distribution of the colony in the image.

To help users understand how to use these functions, we have provided a Jupyter notebook titled "image_pipeline_tutorial_notebook.ipynb". This notebook walks through the process of analyzing some given test data using the functions in image_functions.jl. It provides a practical example of how these functions can be used to extract meaningful information from colony images. 

### Artificial Colony Creation

The `artificial_colony_creation.jl` file contains functions for creating artificial colony images. These images can be used as test data for the image functions.

### Usage
Here is a basic example of how to use the expand_colony_radom_cov! function:

```julia
using ColonyImages

img = zeros(100, 100)
img[50:55, 50:55] .= 1
expand_colony_radom_cov!(img, 100)
```
### Documentation
For more detailed information about the functions and their usage, please refer to the [documentation](https://andreaskuhn-ak.github.io/ColonyImages.jl/).

### Contributing
Contributions to ColonyImages are very welcome! If you have a feature request, bug report, or proposal, please open an issue on the GitHub repository.

### License
ColonyImages is licensed under the MIT license.
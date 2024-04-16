var documenterSearchIndex = {"docs":
[{"location":"#ColonyImages","page":"Home","title":"ColonyImages","text":"","category":"section"},{"location":"","page":"Home","title":"Home","text":"Documentation for ColonyImages.","category":"page"},{"location":"","page":"Home","title":"Home","text":"    Modules = [ColonyImages]","category":"page"},{"location":"#ColonyImages.analysis_parameters","page":"Home","title":"ColonyImages.analysis_parameters","text":"analysis_parameters\n\nA struct for holding parameters related to the analysis.\n\nFields\n\nplot_theme::Attributes: An Attributes object for setting the theme of the plots. Default is a Theme object with a fontsize of 25, a size of (1000,800), a markersize of 15 for Scatter plots, and a linewidth of 4 for Lines plots.\n\nExample\n\nparams = analysis_parameters()\n\n\n\n\n\n","category":"type"},{"location":"#ColonyImages.parameters","page":"Home","title":"ColonyImages.parameters","text":"parameters\n\nA struct that holds various parameters for colony image processing, creation and analysis.\n\nFields\n\ntime_points::Vector{Float64}: The time points for the analysis. Default is [0,48,96,144].\nthreshold_conv::Float64: The threshold for the convolution. Default is 0.8.\nthreshold_c::Float64: Another threshold parameter. Default is 0.8.\nkernel_ratio::Float64: The ratio for the kernel. Default is 0.4.\nsteps_angular::Int: The number of angular steps. Default is 360.\nsamples_pair::Int: The number of sample pairs. Default is 2000000.\ntimesteps::Int64: The number of time steps. Default is 300.\nim_size::Vector{Int}: The size of the image. Default is [600,600].\nstacks::Int: The number of stacks, which is the length of time_points.\nradius_colony::Int: The radius of the colony. Default is round(Int,(im_size[1]*0.05)).\nCenter::Vector{Int}: The center of the image. Default is round.(Int,im_size./2).\ngrowth_rate::Float64: The growth rate of the colony. Default is 0.02971700864000873.\ncolony_size::Function: A function to calculate the size of the colony. Default is t-> (1+growth_rate).^t.\nrelative_size_filles_holes::Float64: The relative size of filled holes. Default is 0.01.\nlaplac_kernel::Matrix{Int}: The Laplacian kernel. Default is [0 1 0; 1 -4 1; 0 1 0].\ncolony_nr::Int: The number of colonies. Default is 4.\ncolonies::Vector{String}: The names of the colonies. Default is [\"Colony x artifical\" for x in 1:colony_nr].\nplot_factor::AbstractFloat: The factor for plotting. Default is 2.0.\nPoints::Vector{Vector{Vector{Int}}}: The lattice points. Default is lattice_points(Int(maximum(im_size)÷2)).\ncol_size_add::Vector{Float64}: The additional size of the colony. Default is colony_size.(time_points).-1.\ncol_size_add_diff::Vector{Float64}: The difference in the additional size of the colony. Default is col_size_add[2:end]-col_size_add[1:end-1].\n\n\n\n\n\n","category":"type"},{"location":"#ColonyImages.angular_metric-Tuple{Union{BitMatrix, Matrix{<:Real}}, Vector{Int64}}","page":"Home","title":"ColonyImages.angular_metric","text":"angular_metric(img::Union{Matrix{<:Real}, BitMatrix}, center::Vector{Int}; steps::Int = 360)\n\nCalculates an angular metric for a given image. The metric is a vector where each element represents the number of pixels in a certain angular sector of the image.  The sectors are determined by dividing a circle centered at a given point into a certain number of equal parts.\n\nArguments\n\nimg::Union{Matrix{<:Real}, BitMatrix}: The input image.\ncenter::Vector{Int}: The center of the circle.\nsteps::Int: The number of sectors into which the circle is divided. Default is 360.\n\nReturns\n\nA vector where each element represents the number of pixels in a certain angular sector of the image.\n\n\n\n\n\n","category":"method"},{"location":"#ColonyImages.approx_radi_colo-Tuple{Union{BitMatrix, Matrix{<:Real}}}","page":"Home","title":"ColonyImages.approx_radi_colo","text":"approx_radi_colo(img::Union{Matrix{<:Real}, BitMatrix})\n\nCalculates the approximate diameter of a colony by summing up all the pixel values and taking the square root of the sum.\n\nThis function assumes that the pixel values represent the area of the colony. The diameter is then approximated using the formula for the diameter of a circle given its area.\n\nArguments\n\nimg::Union{Matrix{<:Real}, BitMatrix}: A 2D array representing the image. The pixel values are assumed to represent the area of the colony.\n\nReturns\n\nA float representing the approximate diameter of the colony.\n\n\n\n\n\n","category":"method"},{"location":"#ColonyImages.b_w-Tuple{Any}","page":"Home","title":"ColonyImages.b_w","text":"b_w(img)\n\nConverts a grayscale colony image into a binary image/BitArray. If more than half of the image is black, it inverts the image.  This ensures that in the output image, the pixels inside the colony are always set to 1 and the background pixels to 0,  regardless of the inversion status of the input image.\n\nArguments\n\nimg: The input image.\n\nReturns\n\nA binary image where the colony pixels are set to 1 and the background pixels are set to 0.\n\n\n\n\n\n","category":"method"},{"location":"#ColonyImages.build_artifical_colony!-Tuple{Vector{Int64}, AbstractArray, Int64, Vector{Vector{Vector{Int64}}}}","page":"Home","title":"ColonyImages.build_artifical_colony!","text":"build_artifical_colony!(center::Vector{Int}, img::AbstractArray}, radius::Int, points::Vector{Vector{Vector{Int}}})\n\nThis function constructs an artificial spherical colony within a given image. The colony is represented as a circle with a specified center and radius.  The function directly modifies the input image.\n\nArguments\n\ncenter::Vector{Int}: A vector representing the center coordinates of the colony.\nimg::AbstractArray: The input image where the colony will be built.\nradius::Int: The radius of the colony.\npoints::Vector{Vector{Vector{Int}}}: A nested vector containing the points used to construct the colony.\n\nReturns\n\nThe image with the built colony.\n\n\n\n\n\n","category":"method"},{"location":"#ColonyImages.build_circle-Tuple{Vector{Int64}, Union{BitMatrix, Matrix{<:Real}}, Vector{Vector{Vector{Int64}}}}","page":"Home","title":"ColonyImages.build_circle","text":"build_circle(center::Vector{Int}, img::Union{Matrix{<:Real}, BitMatrix}, points::Vector{Vector{Vector{Int}}}; threshold = 0.8::Float64)\n\nCreates a binary image with the same size as the input image. The binary image is a circle with a given center.  The circle is built by iterating over a set of points and setting the corresponding pixel in the binary image to 1 if the point is within the circle. The occupation in the outermost circle band is calculated in each iteration and stored in the occupation vector.  The function stops building the circle when the mean of the occupation vector is less than a given threshold.\n\nArguments\n\ncenter::Vector{Int}: The center of the circle.\nimg::Union{Matrix{<:Real}, BitMatrix}: The input image.\npoints::Vector{Vector{Vector{Int}}}: A set of points used to build the circle.\nthreshold::Float64: The threshold for the mean of the occupation vector. Defaults to 0.8.\n\nReturns\n\ncircle_kernel: A binary image representing the circle.\n\n\n\n\n\n","category":"method"},{"location":"#ColonyImages.centroid-Tuple{Union{BitMatrix, Matrix{<:Real}}}","page":"Home","title":"ColonyImages.centroid","text":"centroid(img::Union{Matrix{<:Real}, BitMatrix})\n\nCalculates the centroid of a given image img.\n\nThe centroid is calculated as the average of the x and y coordinates of all non-zero pixels in the image,  weighted by their intensity. The coordinates are then rounded to the nearest integer.\n\nArguments\n\nimg::Union{Matrix{<:Real}, BitMatrix}: A 2D array representing the image. Non-zero values are considered as part of the object to find the centroid of.\n\nReturns\n\ncentroid_norm::Vector{Int64}: A vector containing the x and y coordinates of the centroid.\n\n\n\n\n\n","category":"method"},{"location":"#ColonyImages.conv-Tuple{AbstractArray, AbstractArray}","page":"Home","title":"ColonyImages.conv","text":"conv(img::Union{Matrix{<:Real}, BitMatrix}, kernel::Union{Matrix{<:Real}, BitMatrix})\n\nPerforms a convolution operation on an image using a given kernel. The input image and kernel are 2D arrays of Int or Float or Bool.  The function returns a 2D Float64 array of the same size as the input image.\n\nArguments\n\nimg::Union{Matrix{<:Real}, BitMatrix}: The input image, a 2D array of Int, Float or Bool.\nkernel::Union{Matrix{<:Real}, BitMatrix}: The kernel used for the convolution, a smaller 2D array of Int,Float or Bool.\n\nReturns\n\nA 2D Float64 array representing the convolved image.\n\n\n\n\n\n","category":"method"},{"location":"#ColonyImages.create_kernel-Tuple{Int64}","page":"Home","title":"ColonyImages.create_kernel","text":"create_kernel(rad::Int;  geometry::String = \"circle\")\n\nCreates a binary image kernel with a given radius. The kernel can be either a circle or a square.\n\nArguments\n\nrad::Int: The radius of the kernel.\ngeometry::String: The shape of the kernel. Can be either \"circle\" or \"square\". Defaults to \"circle\".\n\nReturns\n\nkernel: A 2D array representing the kernel.\n\n\n\n\n\n","category":"method"},{"location":"#ColonyImages.expand_colony_circular!-Tuple{AbstractArray, Vector{Vector{Vector{Int64}}}, Vector{Int64}, Int64}","page":"Home","title":"ColonyImages.expand_colony_circular!","text":"expand_colony_circular!(img::AbstractArray, points::Vector{Vector{Vector{Int}}}, center::Vector{Int}, pixels_to_add::Int)\n\nThis function expands an artifical colony in a circular pattern within an image. The expansion starts from the center of the colony and proceeds outward.  The function directly modifies the input image.\n\nArguments\n\nimg::AbstractArray: The input image where the colony will be expanded.\npoints::Vector{Vector{Vector{Int}}}: A nested vector containing the points used to expand the colony.\ncenter::Vector{Int}: A vector representing the center coordinates of the colony.\npixels_to_add::Int: The number of pixels to add to the colony.\n\n\n\n\n\n","category":"method"},{"location":"#ColonyImages.expand_colony_finger_radom_cov!-Tuple{AbstractArray, Int64, Any}","page":"Home","title":"ColonyImages.expand_colony_finger_radom_cov!","text":"expand_colony_finger_radom_cov!(img::AbstractArray, pixels_to_add::Int, dir::Vector{Int}; dir_match_rate::AbstractFloat = 0.999, still_spawn_rate::AbstractFloat = 0.99, min_neigbour::Int = 2)\n\nExpand the colony in the image img by adding pixels_to_add pixels. The expansion is done randomly  at the border of the colony, but it is more likely to expand in the direction specified by dir.  The border is determined by convolving the image with a Laplacian kernel and finding points where the convolution  is greater than 0.1. The function modifies the input image in-place.\n\nArguments\n\nimg::AbstractArray: A 2D array representing the image of the colony. The colony is represented by 1s and the background by 0s.\npixels_to_add::Int: The number of pixels to add to the colony.\ndir::Vector{Int}: A vector representing the preferred direction of expansion.\ndir_match_rate::AbstractFloat: A float representing the threshold for the dot product between the preferred direction and the point. Default is 0.999.\nstill_spawn_rate::AbstractFloat: A float representing the probability of expanding in a direction opposite to dir. Default is 0.99.\nmin_neigbour::Int: The minimum number of neighboring points that must be occupied for a point to be added to the colony. Default is 2.\n\nExample\n\nimg = zeros(100, 100)\nimg[50:55, 50:55] .= 1\nexpand_colony_finger_radom_cov!(img, 100, [1, 0])\n\n\n\n\n\n","category":"method"},{"location":"#ColonyImages.expand_colony_point!-Tuple{AbstractArray, AbstractArray, CartesianIndex{2}}","page":"Home","title":"ColonyImages.expand_colony_point!","text":"expand_colony_point!(img::AbstractArray, cov_img::AbstractArray, point::Vector{Int})\n\nExpands a colony at a given point in an image. The function uses a Laplacian kernel to convolve the image and find border points.  If the given point is in the border points and is in the background of the image, it is added to the colony and the convolution image is updated.\n\nArguments\n\nimg::AbstractArray: The input image.\ncov_img::AbstractArray: The convolution image.\npoint::Vector{Int}: The point at which to expand the colony.\n\nReturns\n\nThe updated convolution image.\n\n\n\n\n\n","category":"method"},{"location":"#ColonyImages.expand_colony_radom!-Tuple{AbstractArray, Int64}","page":"Home","title":"ColonyImages.expand_colony_radom!","text":"expand_colony_radom!(img::AbstractArray, pixels_to_add::Int)\n\nThis function expands a colony in a random pattern within an image. The expansion is performed by adding pixels to the border of the colony.  The function directly modifies the input image. \n\nArguments\n\nimg::AbstractArray: The input image where the colony will be expanded.\npixels_to_add::Int: The number of pixels to add to the colony.\n\nCompared to expand_colony_radom_cov!, this function is slower for large images and many pixels to add,  but faster for small images and and fewer pixels to add .This is due to the fact that the computational heavy convolution only needs to be      calculated once for the whole image, whereas the distance transform in expand_colony_radom needs to be calculated for each iteration of the loop.\n\n\n\n\n\n","category":"method"},{"location":"#ColonyImages.expand_colony_radom_cov!-Tuple{AbstractArray, Int64}","page":"Home","title":"ColonyImages.expand_colony_radom_cov!","text":"expand_colony_radom_cov!(img::AbstractArray, pixels_to_add::Int)\n\nExpand the colony in the image img by adding pixels_to_add pixels. The expansion is done randomly  at the border of the colony. The border is determined by convolving the image with a Laplacian kernel  and finding points where the convolution is greater than 0.1. The function modifies the input image in-place.\n\nCompared to expand_colony_radom!, this function is faster for large images and many pixels to add,  but slower for small images and fewer pixels to add. This is due to the fact that the computationally heavy convolution only needs to be  calculated once for the whole image, whereas the distance transform in expand_colony_radom needs to be calculated for each iteration of the loop.\n\nArguments\n\nimg::AbstractArray: A 2D array representing the image of the colony. The colony is represented by 1s and the background by 0s.\npixels_to_add::Int: The number of pixels to add to the colony.\n\nExample\n\nimg = zeros(100, 100)\nimg[50:55, 50:55] .= 1\nexpand_colony_radom_cov!(img, 100)\n\nCompared to expand_colony_radom!, this function is faster for large images and many pixels to add,  but slower for small images and fewer pixels to add. This is due to the fact that the computational heavy convolution only needs to be      calculated once for the whole image, whereas the distance transform in expand_colony_radom needs to be calculated for each iteration of the loop.\n\n\n\n\n\n","category":"method"},{"location":"#ColonyImages.expand_colony_un_random_cov!-Tuple{AbstractArray, Int64, Vector{Int64}}","page":"Home","title":"ColonyImages.expand_colony_un_random_cov!","text":"expand_colony_un_radom_cov!(img::AbstractArray, pixels_to_add::Int, dir::Vector{Int}; still_spawn_rate::AbstractFloat = 0.5)\n\nThis function expands a colony in a given direction by adding pixels to the image.\n\nArguments\n\nimg::AbstractArray: A 2D array representing the image of the colony.\npixels_to_add::Int: The number of pixels to add to the colony.\ndir::Vector{Int}: A vector representing the direction in which to expand the colony.\nstill_spawn_rate::AbstractFloat: The probability of adding a pixel even if it is not in the desired direction. Default is 0.5.\n\nExample\n\nimg = zeros(100, 100)\nimg[50:55, 50:55] .= 1\nexpand_colony_un_radom_cov!(img, 100, [1, 0])\n\n\n\n\n\n","category":"method"},{"location":"#ColonyImages.expand_matrix-Tuple{Union{BitMatrix, Matrix{<:Real}}}","page":"Home","title":"ColonyImages.expand_matrix","text":"expand_matrix(mat::Union{Matrix{<:Real}, BitMatrix}; annuli::Int = 2 )\n\nExpands a matrix by adding zero-filled columns to its outermost right column and rotating it 90 degrees counterclockwise. This process is repeated 4 times per annulus. The number of additional bands (annuli) added to the matrix is determined by the annuli parameter.\n\nArguments mat::Union{Matrix{<:Real}, BitMatrix}: The input matrix. annuli::Int: The number of additional bands to be added to the matrix. Default is 2. Returns The expanded matrix. Examples\n\nmat = [1 2; 3 4]\nexpand_matrix(mat, annuli = 1)\n4×4 Matrix{Int64}:\n 0  0  0  0\n 0  1  2  0\n 0  3  4  0\n 0  0  0  0\n\n\n\n\n\n","category":"method"},{"location":"#ColonyImages.fill_holes-Tuple{Any, Real}","page":"Home","title":"ColonyImages.fill_holes","text":"fill_holes(img, size_holes::Real)\n\nFills holes in a binary image. The size of the holes to be filled is determined by the size_holes parameter.\n\nArguments\n\n`img: The input binary image.\nsize_holes::Real: The relative size of the holes to be filled. This is a fraction of the total number of pixels in the image.\n\nReturns\n\nA binary image with the holes filled.\n\n\n\n\n\n","category":"method"},{"location":"#ColonyImages.filter_fourier_alpha-Tuple{Vector{<:Real}}","page":"Home","title":"ColonyImages.filter_fourier_alpha","text":"filter_fourier_alpha(vec; a = 5)\n\nFilters a vector by setting all elements after the a-th element to zero. This function is useful for filtering the results of a Fourier transform, where the elements of the vector represent the amplitudes of the frequencies, and the index of the element represents the frequency. By setting all elements after the a-th element to zero, we effectively remove all frequencies higher than a.\n\nArguments\n\nvec::Vector{<:Real}: The input vector.\na::Int: The cutoff frequency. All frequencies higher than a are removed. Default is 5.\n\nReturns\n\nA new vector where all elements after the a-th element are zero.\n\n\n\n\n\n","category":"method"},{"location":"#ColonyImages.filter_fourier_beta-Tuple{Vector{<:Real}}","page":"Home","title":"ColonyImages.filter_fourier_beta","text":"filter_fourier_beta(vec::Vector{<:Real}; b::AbstractFloat = 0.5)\n\nFilters a vector by setting all elements whose absolute value is less than or equal to b times the maximum absolute value of the elements to zero. This function is useful for filtering the results of a Fourier transform, where the elements of the vector represent the amplitudes of the frequencies. By setting all elements whose amplitude is less than or equal to b times the maximum amplitude to zero, we effectively remove all frequencies with low amplitudes.\n\nArguments\n\nvec::Vector{<:Real}: The input vector.\nb::AbstractFloat: The threshold for the amplitude. All frequencies with amplitudes less than or equal to b times the maximum amplitude are removed. Default is 0.5.\n\nReturns\n\nA new vector where all elements whose absolute value is less than or equal to b times the maximum absolute value are zero.\n\n\n\n\n\n","category":"method"},{"location":"#ColonyImages.find_freq-Tuple{Vector{<:Real}}","page":"Home","title":"ColonyImages.find_freq","text":"find_freq(vec::Vector{<:Real}; ignore_latter_half = true )\n\nFinds the frequencies in a vector that have non-zero amplitudes. This function is useful for analyzing the results of a Fourier transform, where the elements of the vector represent the amplitudes of the frequencies, and the index of the element represents the frequency.\n\nArguments\n\nvec::Vector{<:Real}: The input vector.\nignore_latter_half::Bool: If true, the function only considers the first half of the vector. This is useful when the input vector is the result of a Fourier transform, where the second half of the vector contains the same information as the first half but in reverse order. Default is true.\n\nReturns\n\nA vector of the frequencies that have non-zero amplitudes.\n\n\n\n\n\n","category":"method"},{"location":"#ColonyImages.find_peaks-Tuple{Vector{<:Real}}","page":"Home","title":"ColonyImages.find_peaks","text":"find_peaks(signal::Vector{<:Real}; threshold::AbstractFloat = 1.0)\n\nFinds the local maxima inside a 1D signal, in areas where the signal exceeds its mean value by a given factor. \n\nArguments\n\nsignal::Vector{<:Real}: The input 1D signal.\nthreshold::AbstractFloat: The factor by which the signal needs to exceed its mean value to be considered a peak. Default is 1.0.\n\nReturns\n\nposition_peaks: A vector of the positions of the peaks in the signal.\nnr_peaks: The number of peaks found in the signal.\n\n\n\n\n\n","category":"method"},{"location":"#ColonyImages.generate_dir_vec-Tuple{parameters}","page":"Home","title":"ColonyImages.generate_dir_vec","text":"generate_dir_vec(para::parameters)\n\nGenerates vectors of x and y coordinates that span vectors ranging from 0 to 2π.  This function is used to generate the vectors used in expand_colony_finger_radom_cov!.\n\nArguments\n\npara::parameters: A parameters object containing various parameters for the analysis. \nnumber_finger: The number of vectors to generate.\nfinger_dist: The random distance to add to the vectors.\n\nReturns\n\ndir::Array{Array{Float64,1},1}: A vector of vectors, each containing the y and x coordinates of a vector.\n\nExample\n\nusing CairoMakie\npara = parameters(number_finger = 20, finger_dist = 0.1)\ndir = generate_dir_vec(para)\nyy = zeros(para.number_finger)\narrows(yy,yy, [d[1] for d in dir], [d[2] for d in dir])\n\n\n\n\n\n","category":"method"},{"location":"#ColonyImages.lattice_points-Tuple{Int64}","page":"Home","title":"ColonyImages.lattice_points","text":"lattice_points(r::Int)\n\nGenerates a lattice of points within a circle of radius r.\n\nThe function returns a nested vector of points, each represented by its x and y coordinates.  The points are sorted by their distance from the origin and grouped into bands,  each band containing points that have a similar distance to the origin.\n\nArguments\n\nr::Int: The radius of the circle within which the lattice points are generated.\n\nReturns\n\npoints3::Vector{Vector{Vector{Int}}}: A nested vector of points.\n\npoints3 is a vector of vectors of Lattice vectors in ascending order by their length . Whereas the first entry of Points contains a vector of all lattice vectors which lengths are lower than 2:\n\nPoints[1] = [[0, 1], [0, 1], [0, -1], [0, -1], [1, 0], [-1, 0], [1, 0], [-1, 0], [1, 1], [-1, 1], [1, -1], [-1, -1]]\n\nPoints[2] does the same for length lower than 3:\n\nPoints[2] = [[0, 2], [0, 2], [0, -2], [0, -2], [2, 0], [-2, 0], [2, 0], [-2, 0], [1, 2], [-1, 2], [1, -2], [-1, -2], [2, 1], [-2, 1], [2, -1], [-2, -1], [2, 2], [-2, 2], [2, -2], [-2, -2]]\n\nand so on...\n\n\n\n\n\n","category":"method"},{"location":"#ColonyImages.occupied_points-Tuple{Union{BitMatrix, Matrix{<:Real}}}","page":"Home","title":"ColonyImages.occupied_points","text":"occupied_points(img::Union{Matrix{<:Real}, BitMatrix})\n\nCalculates the proportion of occupied points in a given binary image img.\n\nThe function sums up all the pixel values in the image and divides by the total number of pixels.  This gives the proportion of occupied points, assuming that non-zero pixel values represent occupied points.\n\nArguments\n\nimg::Union{Matrix{<:Real}, BitMatrix}: A 2D array representing the image. Non-zero values are considered as occupied points.\n\nReturns\n\nA float representing the proportion of occupied points in the image.\n\n\n\n\n\n","category":"method"},{"location":"#ColonyImages.pair_cor_metric3-Tuple{Union{BitMatrix, Matrix{<:Real}}, Vector{Int64}}","page":"Home","title":"ColonyImages.pair_cor_metric3","text":"pair_cor_metric3(img::Union{Matrix{<:Real}, BitMatrix}, center::Vector{Int}; samples::Int = 10000, steps::Int = 360)\n\nCalculates a pair correlation metric for a given image. The metric is a vector where each element represents the number of pairs of pixels that have a certain relative angle.  The angles are determined by dividing a half-circle into a certain number of equal parts.\n\nArguments\n\nimg::Union{Matrix{<:Real}, BitMatrix}: The input image.\ncenter::Vector{Int}: The reference point for calculating the relative angles.\nsamples::Int: The number of pairs of pixels to sample. Default is 10000.\nsteps::Int: The number of sectors into which the half-circle is divided. Default is 360.\n\nReturns\n\nA vector where each element represents the number of pairs of pixels that have a certain relative angle.\n\n\n\n\n\n","category":"method"},{"location":"#ColonyImages.plot_and_save_time_series_metrics!-Tuple{AbstractArray, parameters, DataFrames.DataFrame, Any}","page":"Home","title":"ColonyImages.plot_and_save_time_series_metrics!","text":"plot_and_save_time_series_metrics!(img_vec::AbstractArray, para::parameters, df::DataFrame; name_data =\"random_growth\")\n\nThis function creates a series of plots for each image in a stack, showing the original image, the angular metric, and the pair correlation metric. It also updates a DataFrame with the metrics.\n\nArguments\n\nimg_vec::AbstractArray: A 4D array where the first two dimensions are the image dimensions, the third dimension is the time point, and the fourth dimension is the image stack.\npara::parameters: A parameters object containing various parameters for the analysis. \nplot_factor: The scaling factor for the plot size.\nkernel_ratio::Float64: The ratio for the kernel operation.\nthreshold_conv: The threshold for the convolution operation.\nsteps_angular: The number of steps for the angular metric calculation.\nsamples_pair: The number of samples for the pair correlation metric calculation.\nPoints: The number of points for the circle fitting.\nthreshold_c: The threshold for the circle fitting.\ncolonies: The names of the colonies.\ndf::DataFrame: A DataFrame to update with the metrics.\nname_data::String: A string to prepend to the date for the data set name in the DataFrame. Default is \"random_growth\".\n\nReturns\n\nfig_big::Figure: A Figure object containing the plots.\n\nDetails\n\nThe function first creates a Figure object with a size determined by the res_scaling function. It then loops over each image stack in img_vec. For each stack, it calculates the centroid of the first image and creates a kernel based on the approximate radius of the colony in the image.\n\nThe function then loops over each time point in the image stack. For each time point, it calculates the centroid of the convoluted image and calculates the angular metric and pair correlation metric for the original image and the image minus the first image in the stack.\n\nThe function then fits a circle to the image and calculates the angular metric and pair correlation metric for the image minus the fitted circle.\n\nThe function then creates three Axis objects for the original image, the angular metric, and the pair correlation metric, and adds plots to each axis. If the time point is not the first, it also plots the angular metric and pair correlation metric for the image minus the fitted circle.\n\nFinally, the function increments a counter, updates the DataFrame with the metrics, and returns the Figure object after looping over all image stacks and time points.\n\n\n\n\n\n","category":"method"},{"location":"#ColonyImages.plot_convolution_schematic2-Tuple{AbstractArray, AbstractArray, parameters}","page":"Home","title":"ColonyImages.plot_convolution_schematic2","text":"plot_convolution_schematic2(colony::AbstractArray, colony_cov::AbstractArray, para::parameters; name = \"cov_fig1\")\n\nThis function creates a plot of a colony image and its convolution. \n\nArguments\n\ncolony::AbstractArray: A 3D array representing the colony image stack or a 2D array representing a single colony image.\ncolony_cov::AbstractArray: A 3D array representing the convoluted colony image stack or a 2D array representing a single convoluted colony image.\npara::parameters: A parameters object containing various parameters for the analysis.\nname::String: An optional name for the output plot. Default is \"cov_fig1\".\n\nReturns\n\ncov_fig::Figure: A Figure object containing the plots.\n\nDetails\n\nThe function first determines the type of the colony and colony_cov inputs and extracts the data for the last time point if they are 3D arrays. \n\nIt then creates a Figure object and two Axis objects for the original colony and the convoluted colony. It adds heatmaps to each of these axes using the respective data.\n\nThe function then iterates over each point in the images and adds a text object to each point in the heatmaps, displaying the value at that point.\n\nFinally, the function saves the figure as a PNG image and returns the Figure object.\n\n\n\n\n\n","category":"method"},{"location":"#ColonyImages.plot_convolution_schematic3-Tuple{AbstractArray, AbstractArray, parameters}","page":"Home","title":"ColonyImages.plot_convolution_schematic3","text":"plot_convolution_schematic3(colony, colony_cov, para; name = \"cov_fig_kernel\")\n\nThis function creates a schematic plot of a convolution operation on a colony image. \n\nArguments\n\ncolony::Array: A 3D array representing the colony image stack.\ncolony_cov::Array: A 3D array representing the convoluted colony image stack.\npara::parameters: A parameters object containing various parameters for the analysis.\nname::String: An optional name for the output plot. Default is \"covfigkernel\".\n\nReturns\n\ncov_fig_k::Figure: A Figure object containing the plots.\n\nDetails\n\nThe function first sets up a theme for the plot with a specific font size. It then extracts the data for the last time point from the colony and convoluted images. It creates a kernel and places it in a larger matrix of zeros.\n\nThe function then creates a Figure object and three Axis objects for the original colony, the kernel, and the convoluted colony. It adds heatmaps to each of these axes using the respective data.\n\nThe function then iterates over each point in the images and adds a text object to each point in the heatmaps, displaying the value at that point.\n\nFinally, the function saves the figure as a PNG image and returns the Figure object.\n\n\n\n\n\n","category":"method"},{"location":"#ColonyImages.plot_metric_schematic-Tuple{AbstractArray, parameters, Any}","page":"Home","title":"ColonyImages.plot_metric_schematic","text":"plot_metric_schematic(colony::AbstractArray, para::parameters, colormap ; name = \"metric_fig1\")\n\nThis function creates a schematic plot of a colony image and how it is cut into angle section to be analysed later. \n\nArguments\n\ncolony::AbstractArray: A  3D array representing the colony image stack.\npara::parameters: A parameters object containing various parameters for the analysis.\ncolormap: A colormap to use for the heatmap.\nname::String: An optional name for the output plot. Default is \"metric_fig1\".\n\nReturns\n\nfig::Figure: A Figure object containing the plots.\n\nDetails\n\nThe function first creates a Figure object and extracts the first and last images from the colony stack. It creates a kernel based on the approximate radius of the colony and the kernel ratio parameter.\n\nThe function then cuts a portion of the last image and rotates it 90 degrees. It creates two Axis objects on the Figure and hides their decorations.\n\nThe function then cuts a portion of the first image and rotates it 90 degrees. It adds a heatmap to the first Axis, showing the difference between the last and first images.\n\nThe function then creates a pie chart on the second Axis, with each slice having an equal value, representing the angle sections of the colony.\n\nFinally, the function returns the Figure object.\n\n\n\n\n\n","category":"method"},{"location":"#ColonyImages.plot_time_series_cov_centroid-Tuple{AbstractArray, parameters}","page":"Home","title":"ColonyImages.plot_time_series_cov_centroid","text":"plot_time_series_cov_centroid(img_vec::AbstractArray, para::parameters)\n\nGenerates a time series plot of the centroids of the images in img_vec. The centroids are calculated in three ways: \n\nThe original centroid of the first timestep of the timeseries.\nThe centroid after applying a convolution operation.\nThe current centroid of the image.\n\nArguments\n\nimg_vec::AbstractArray: A vector of 3D image stacks.\npara::parameters: A parameters object containing various parameters for the analysis. \nkernel_ratio::Float64: The ratio for the kernel operation. Default is 0.4.\nplot_factor: The scaling factor for the plot size.\nthreshold_conv: The threshold for the convolution operation.\ncolonies: The names of the colonies.\n\nReturns\n\nfig_big::Figure: A Figure object with the time series plot of the centroids.\n\nExample\n\nimg_vec = [rand(100, 100, 10) for _ in 1:5]\npara = parameters(kernel_ratio = 0.4, plot_factor = 1, threshold_conv = 0.5, colonies = [\"Colony i\" for i in 1:5])\nfig = plot_time_series_cov_centroid(img_vec, para)\n\n\n\n\n\n","category":"method"},{"location":"#ColonyImages.plot_time_series_metrics-Tuple{AbstractArray, parameters}","page":"Home","title":"ColonyImages.plot_time_series_metrics","text":"plot_time_series_metrics(img_vec, para)\n\nThis function creates a series of plots for each image in a stack, showing the original image, the angular metric, and the pair correlation metric. \n\nArguments\n\nimg_vec::Array: A 4D array where the first two dimensions are the image dimensions, the third dimension is the time point, and the fourth dimension is the image stack.\npara::parameters: A parameters object containing various parameters for the analysis.\n\nReturns\n\nfig_big::Figure: A Figure object containing the plots.\n\nDetails\n\nThe function first creates a Figure object with a size determined by the res_scaling function. It then loops over each image stack in img_vec. For each stack, it calculates the centroid of the first image and creates a kernel based on the approximate radius of the colony in the image.\n\nThe function then loops over each time point in the image stack. For each time point, it calculates the centroid of the convoluted image and calculates the angular metric and pair correlation metric for the original image and the image minus the first image in the stack.\n\nThe function then fits a circle to the image and calculates the angular metric and pair correlation metric for the image minus the fitted circle.\n\nThe function then creates three Axis objects for the original image, the angular metric, and the pair correlation metric, and adds plots to each axis. If the time point is not the first, it also plots the angular metric and pair correlation metric for the image minus the fitted circle.\n\nFinally, the function increments a counter and returns the Figure object after looping over all image stacks and time points.\n\n\n\n\n\n","category":"method"},{"location":"#ColonyImages.plot_timeseries_heatmap-Tuple{AbstractArray, parameters}","page":"Home","title":"ColonyImages.plot_timeseries_heatmap","text":"plot_timeseries_heatmap(colony::AbstractArray, para::parameters; name = \"Eden Growth Model\", colormap  = reverse(co.Blues))\n\nThis function creates a heatmap plot of a colony image stack over time.\n\nArguments\n\ncolony::AbstractArray: A 3D array representing the colony image stack.\npara::parameters: A parameters object containing various parameters for the analysis.\nname::String: An optional name for the output plot. Default is \"Eden Growth Model\".\ncolormap: An optional colormap to use for the heatmap. Default is the reverse of Blues.\n\nReturns\n\nfig_eden::Figure: A Figure object containing the heatmap plot.\n\nDetails\n\nThe function first creates a Figure object and an Axis object with the title set to the provided name. It then initializes an intensity image with zeros.\n\nThe function then iterates over the z-dimension of the colony stack, adding each image to the intensity image.\n\nIt then creates a heatmap on the Axis, showing the intensity image with the maximum intensity subtracted and the sign reversed.\n\nThe function then hides the decorations of the Axis and adds a colorbar to the Figure, with the ticks set to the time points and the label set to \"time [h]\".\n\nFinally, the function saves the Figure as a PDF in the plots folder and returns the Figure object.\n\n\n\n\n\n","category":"method"},{"location":"#ColonyImages.rad2deg_discrete-Tuple{AbstractFloat}","page":"Home","title":"ColonyImages.rad2deg_discrete","text":"rad2deg_discrete(ϕ::AbstractFloat; steps::Int =360)\n\nConverts an angle in radians to a discrete angle in degrees. The output is the number of a circular sector on a unit circle divided into a given number of sectors.  For example, if the angle is 0.0 and the number of steps is 360, the output is 1, which corresponds to the first circular sector on the unit circle.  If the angle is (2pi - 0.01) and the number of steps is 360, the output is 360, which corresponds to the last circular sector on the unit circle.\n\nArguments\n\nϕ::AbstractFloat: The input angle in radians.\nsteps::Int: The number of circular sectors into which the unit circle is divided. Default is 360.\n\nReturns\n\nThe number of the circular sector to which the angle corresponds.\n\n\n\n\n\n","category":"method"},{"location":"#ColonyImages.res_scaling-Tuple{AbstractArray}","page":"Home","title":"ColonyImages.res_scaling","text":"res_scaling(img_int_vec; factor = 3, plots = 1)\n\nScales the resolution of a plot based on a given factor. The function counts the number of images in the given image_vec and and scales the resolution  of plot containg all these images accordingly.\n\nArguments\n\nimg_int_vec: A vector of images.\nfactor: The scaling factor. Defaults to 3.\nplots: The number of plots per images. Defaults to 1.\n\nReturns\n\nA tuple containing the scaled width and height of the image.\n\n\n\n\n\n","category":"method"},{"location":"#ColonyImages.save_time_series_metrics!-Tuple{AbstractArray, parameters, DataFrames.DataFrame, Any}","page":"Home","title":"ColonyImages.save_time_series_metrics!","text":"save_time_series_metrics!(img_vec::AbstractArray, para::parameters, df::DataFrame; name_data=\"random_growth\")\n\nThis function calculates and saves time series metrics for a given set of images.\n\nArguments\n\nimg_vec::AbstractArray: An array of image stacks.\npara::parameters: An object containing various parameters for the analysis.\ndf::DataFrame: The DataFrame to which the calculated metrics will be appended.\nname_data::String: (optional) A string representing the name of the data set. Default is \"random_growth\".\n\nReturns\n\ndata_set::String: The name of the data set.\n\n\n\n\n\n","category":"method"},{"location":"#ColonyImages.@h-Tuple{Any}","page":"Home","title":"ColonyImages.@h","text":"@h methodname\n\nOutputs documentations in jupyternotenbooks in VScode as markdown without bugs.\n\nExample of how to use the @h macro:\n\n@h res_scaling(img_int_vec; factor = 3, plots = 1)\n\nOutputs documentations in jupyternotenbooks in VScode as markdown without bugs.\n\n\n\n\n\n","category":"macro"}]
}

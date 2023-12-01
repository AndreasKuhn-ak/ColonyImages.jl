var documenterSearchIndex = {"docs":
[{"location":"#ColonyImages","page":"Home","title":"ColonyImages","text":"","category":"section"},{"location":"","page":"Home","title":"Home","text":"Documentation for ColonyImages.","category":"page"},{"location":"","page":"Home","title":"Home","text":"    Modules = [ColonyImages]","category":"page"},{"location":"#ColonyImages.angular_metric-Tuple{Union{BitMatrix, Matrix{<:Real}}, Vector{Int64}}","page":"Home","title":"ColonyImages.angular_metric","text":"angular_metric(img::Union{Matrix{<:Real}, BitMatrix}, center::Vector{Int}; steps::Int = 360)\n\nCalculates an angular metric for a given image. The metric is a vector where each element represents the number of pixels in a certain angular sector of the image.  The sectors are determined by dividing a circle centered at a given point into a certain number of equal parts.\n\nArguments\n\nimg::Union{Matrix{<:Real}, BitMatrix}: The input image.\ncenter::Vector{Int}: The center of the circle.\nsteps::Int: The number of sectors into which the circle is divided. Default is 360.\n\nReturns\n\nA vector where each element represents the number of pixels in a certain angular sector of the image.\n\n\n\n\n\n","category":"method"},{"location":"#ColonyImages.approx_radi_colo-Tuple{Union{BitMatrix, Matrix{<:Real}}}","page":"Home","title":"ColonyImages.approx_radi_colo","text":"approx_radi_colo(img::Union{Matrix{<:Real}, BitMatrix})\n\nCalculates the approximate diameter of a colony by summing up all the pixel values and taking the square root of the sum.\n\nThis function assumes that the pixel values represent the area of the colony. The diameter is then approximated using the formula for the diameter of a circle given its area.\n\nArguments\n\nimg::Union{Matrix{<:Real}, BitMatrix}: A 2D array representing the image. The pixel values are assumed to represent the area of the colony.\n\nReturns\n\nA float representing the approximate diameter of the colony.\n\n\n\n\n\n","category":"method"},{"location":"#ColonyImages.b_w-Tuple{Any}","page":"Home","title":"ColonyImages.b_w","text":"b_w(img)\n\nConverts a grayscale colony image into a binary image/BitArray. If more than half of the image is black, it inverts the image.  This ensures that in the output image, the pixels inside the colony are always set to 1 and the background pixels to 0,  regardless of the inversion status of the input image.\n\nArguments\n\nimg: The input image.\n\nReturns\n\nA binary image where the colony pixels are set to 1 and the background pixels are set to 0.\n\n\n\n\n\n","category":"method"},{"location":"#ColonyImages.build_artifical_colony!-Tuple{Vector{Int64}, AbstractArray, Int64, Vector{Vector{Vector{Int64}}}}","page":"Home","title":"ColonyImages.build_artifical_colony!","text":"build_artifical_colony!(center::Vector{Int}, img::AbstractArray}, radius::Int, points::Vector{Vector{Vector{Int}}})\n\nThis function constructs an artificial spherical colony within a given image. The colony is represented as a circle with a specified center and radius.  The function directly modifies the input image.\n\nArguments\n\ncenter::Vector{Int}: A vector representing the center coordinates of the colony.\nimg::AbstractArray: The input image where the colony will be built.\nradius::Int: The radius of the colony.\npoints::Vector{Vector{Vector{Int}}}: A nested vector containing the points used to construct the colony.\n\nReturns\n\nThe image with the built colony.\n\n\n\n\n\n","category":"method"},{"location":"#ColonyImages.build_circle-Tuple{Vector{Int64}, Union{BitMatrix, Matrix{<:Real}}, Vector{Vector{Vector{Int64}}}}","page":"Home","title":"ColonyImages.build_circle","text":"build_circle(center::Vector{Int}, img::Union{Matrix{<:Real}, BitMatrix}, points::Vector{Vector{Vector{Int}}}; threshold = 0.8::Float64)\n\nCreates a binary image with the same size as the input image. The binary image is a circle with a given center.  The circle is built by iterating over a set of points and setting the corresponding pixel in the binary image to 1 if the point is within the circle. The occupation in the outermost circle band is calculated in each iteration and stored in the occupation vector.  The function stops building the circle when the mean of the occupation vector is less than a given threshold.\n\nArguments\n\ncenter::Vector{Int}: The center of the circle.\nimg::Union{Matrix{<:Real}, BitMatrix}: The input image.\npoints::Vector{Vector{Vector{Int}}}: A set of points used to build the circle.\nthreshold::Float64: The threshold for the mean of the occupation vector. Defaults to 0.8.\n\nReturns\n\ncircle_kernel: A binary image representing the circle.\n\n\n\n\n\n","category":"method"},{"location":"#ColonyImages.centroid-Tuple{Union{BitMatrix, Matrix{<:Real}}}","page":"Home","title":"ColonyImages.centroid","text":"centroid(img::Union{Matrix{<:Real}, BitMatrix})\n\nCalculates the centroid of a given image img.\n\nThe centroid is calculated as the average of the x and y coordinates of all non-zero pixels in the image,  weighted by their intensity. The coordinates are then rounded to the nearest integer.\n\nArguments\n\nimg::Union{Matrix{<:Real}, BitMatrix}: A 2D array representing the image. Non-zero values are considered as part of the object to find the centroid of.\n\nReturns\n\ncentroid_norm::Vector{Int64}: A vector containing the x and y coordinates of the centroid.\n\n\n\n\n\n","category":"method"},{"location":"#ColonyImages.conv-Tuple{AbstractArray, AbstractArray}","page":"Home","title":"ColonyImages.conv","text":"conv(img::Union{Matrix{<:Real}, BitMatrix}, kernel::Union{Matrix{<:Real}, BitMatrix})\n\nPerforms a convolution operation on an image using a given kernel. The input image and kernel are 2D arrays of Int or Float or Bool.  The function returns a 2D Float64 array of the same size as the input image.\n\nArguments\n\nimg::Union{Matrix{<:Real}, BitMatrix}: The input image, a 2D array of Int, Float or Bool.\nkernel::Union{Matrix{<:Real}, BitMatrix}: The kernel used for the convolution, a smaller 2D array of Int,Float or Bool.\n\nReturns\n\nA 2D Float64 array representing the convolved image.\n\n\n\n\n\n","category":"method"},{"location":"#ColonyImages.create_kernel-Tuple{Int64}","page":"Home","title":"ColonyImages.create_kernel","text":"create_kernel(rad::Int;  geometry::String = \"circle\")\n\nCreates a binary image kernel with a given radius. The kernel can be either a circle or a square.\n\nArguments\n\nrad::Int: The radius of the kernel.\ngeometry::String: The shape of the kernel. Can be either \"circle\" or \"square\". Defaults to \"circle\".\n\nReturns\n\nkernel: A 2D array representing the kernel.\n\n\n\n\n\n","category":"method"},{"location":"#ColonyImages.expand_colony_circular!-Tuple{AbstractArray, Vector{Vector{Vector{Int64}}}, Vector{Int64}, Int64}","page":"Home","title":"ColonyImages.expand_colony_circular!","text":"expand_colony_circular!(img::AbstractArray, points::Vector{Vector{Vector{Int}}}, center::Vector{Int}, pixels_to_add::Int)\n\nThis function expands an artifical colony in a circular pattern within an image. The expansion starts from the center of the colony and proceeds outward.  The function directly modifies the input image.\n\nArguments\n\nimg::AbstractArray: The input image where the colony will be expanded.\npoints::Vector{Vector{Vector{Int}}}: A nested vector containing the points used to expand the colony.\ncenter::Vector{Int}: A vector representing the center coordinates of the colony.\npixels_to_add::Int: The number of pixels to add to the colony.\n\n\n\n\n\n","category":"method"},{"location":"#ColonyImages.expand_colony_finger_radom_cov!-Tuple{AbstractArray, Int64, Any}","page":"Home","title":"ColonyImages.expand_colony_finger_radom_cov!","text":"expand_colony_finger_radom_cov!(img::AbstractArray, pixels_to_add::Int, dir::Vector{Int}; dir_match_rate::AbstractFloat = 0.999, still_spawn_rate::AbstractFloat = 0.99, min_neigbour::Int = 2)\n\nExpand the colony in the image img by adding pixels_to_add pixels. The expansion is done randomly  at the border of the colony, but it is more likely to expand in the direction specified by dir.  The border is determined by convolving the image with a Laplacian kernel and finding points where the convolution  is greater than 0.1. The function modifies the input image in-place.\n\nArguments\n\nimg::AbstractArray: A 2D array representing the image of the colony. The colony is represented by 1s and the background by 0s.\npixels_to_add::Int: The number of pixels to add to the colony.\ndir::Vector{Int}: A vector representing the preferred direction of expansion.\ndir_match_rate::AbstractFloat: A float representing the threshold for the dot product between the preferred direction and the point. Default is 0.999.\nstill_spawn_rate::AbstractFloat: A float representing the probability of expanding in a direction opposite to dir. Default is 0.99.\nmin_neigbour::Int: The minimum number of neighboring points that must be occupied for a point to be added to the colony. Default is 2.\n\nExample\n\nimg = zeros(100, 100)\nimg[50:55, 50:55] .= 1\nexpand_colony_finger_radom_cov!(img, 100, [1, 0])\n\n\n\n\n\n","category":"method"},{"location":"#ColonyImages.expand_colony_fractal_radom_cov!-Tuple{AbstractArray, Int64}","page":"Home","title":"ColonyImages.expand_colony_fractal_radom_cov!","text":"expand_colony_fractal_radom_cov!(img::AbstractArray, pixels_to_add::Int; max_neigbours::Int, = 1 prob_neigbour::AbstractFloat = 0.3)\n\nExpand the colony in the image img by adding pixels_to_add pixels. The expansion is done randomly  at the border of the colony, but it is more likely to expand to points that only have max_neigbours or less.  The border is determined by convolving the image with a Laplacian kernel and finding points where the convolution  is greater than 0.1. The function modifies the input image in-place.\n\nArguments\n\nimg::AbstractArray: A 2D array of type Matrix with elements of any subtype of Real or a BitMatrix, representing the image of the colony. The colony is represented by 1s and the background by 0s.\npixels_to_add::Int: An integer representing the number of pixels to add to the colony.\nmax_neigbours: An integer representing the maximum number of neighbours a point can have to be considered for expansion. Default is 1.\nprob_neigbour::AbstractFloat: A float representing the probability of expanding to a point with more than max_neigbours. Default is 0.3.\n\nExample\n\nimg = zeros(100, 100)\nimg[50:55, 50:55] .= 1\nexpand_colony_fractal_radom_cov!(img, 100)\n\n\n\n\n\n","category":"method"},{"location":"#ColonyImages.expand_colony_point!-Tuple{AbstractArray, AbstractArray, Vector{Int64}}","page":"Home","title":"ColonyImages.expand_colony_point!","text":"expand_colony_point!(img::AbstractArray, cov_img::AbstractArray, point::Vector{Int})\n\nExpands a colony at a given point in an image. The function uses a Laplacian kernel to convolve the image and find border points.  If the given point is in the border points and is in the background of the image, it is added to the colony and the convolution image is updated.\n\nArguments\n\nimg::AbstractArray: The input image.\ncov_img::AbstractArray: The convolution image.\npoint::Vector{Int}: The point at which to expand the colony.\n\nReturns\n\nThe updated convolution image.\n\n\n\n\n\n","category":"method"},{"location":"#ColonyImages.expand_colony_radom!-Tuple{AbstractArray, Int64}","page":"Home","title":"ColonyImages.expand_colony_radom!","text":"expand_colony_radom!(img::AbstractArray, pixels_to_add::Int)\n\nThis function expands a colony in a random pattern within an image. The expansion is performed by adding pixels to the border of the colony.  The function directly modifies the input image. \n\nArguments\n\nimg::AbstractArray: The input image where the colony will be expanded.\npixels_to_add::Int: The number of pixels to add to the colony.\n\nCompared to expand_colony_radom_cov!, this function is slower for large images and many pixels to add,  but faster for small images and and fewer pixels to add .This is due to the fact that the computational heavy convolution only needs to be      calculated once for the whole image, whereas the distance transform in expand_colony_radom needs to be calculated for each iteration of the loop.\n\n\n\n\n\n","category":"method"},{"location":"#ColonyImages.expand_colony_radom_cov!-Tuple{AbstractArray, Int64}","page":"Home","title":"ColonyImages.expand_colony_radom_cov!","text":"expand_colony_radom_cov!(img::AbstractArray, pixels_to_add::Int)\n\nExpand the colony in the image img by adding pixels_to_add pixels. The expansion is done randomly  at the border of the colony. The border is determined by convolving the image with a Laplacian kernel  and finding points where the convolution is greater than 0.1. The function modifies the input image in-place.\n\nCompared to expand_colony_radom!, this function is faster for large images and many pixels to add,  but slower for small images and fewer pixels to add. This is due to the fact that the computationally heavy convolution only needs to be  calculated once for the whole image, whereas the distance transform in expand_colony_radom needs to be calculated for each iteration of the loop.\n\nArguments\n\nimg::AbstractArray: A 2D array representing the image of the colony. The colony is represented by 1s and the background by 0s.\npixels_to_add::Int: The number of pixels to add to the colony.\n\nExample\n\nimg = zeros(100, 100)\nimg[50:55, 50:55] .= 1\nexpand_colony_radom_cov!(img, 100)\n\nCompared to expand_colony_radom!, this function is faster for large images and many pixels to add,  but slower for small images and fewer pixels to add. This is due to the fact that the computational heavy convolution only needs to be      calculated once for the whole image, whereas the distance transform in expand_colony_radom needs to be calculated for each iteration of the loop.\n\n\n\n\n\n","category":"method"},{"location":"#ColonyImages.expand_colony_radom_cov_show!-Tuple{AbstractArray, Int64}","page":"Home","title":"ColonyImages.expand_colony_radom_cov_show!","text":"expand_colony_radom_cov_show!(img::AbstractArray, pixels_to_add::Int)\n\nExpand the colony in the image img by adding pixels_to_add pixels. The expansion is done randomly  at the border of the colony. The border is determined by convolving the image with a Laplacian kernel  and finding points where the convolution is greater than 0.1. The function modifies the input image in-place.\n\nCompared to expand_colony_radom!, this function is faster for large images and many pixels to add,  but slower for small images and fewer pixels to add. This is due to the fact that the computationally heavy convolution only needs to be  calculated once for the whole image, whereas the distance transform in expand_colony_radom needs to be calculated for each iteration of the loop.\n\nArguments\n\nimg::AbstractArray: A 2D array representing the image of the colony. The colony is represented by 1s and the background by 0s.\npixels_to_add::Int: The number of pixels to add to the colony.\n\nExample\n\nimg = zeros(100, 100)\nimg[50:55, 50:55] .= 1\nexpand_colony_radom_cov!(img, 100)\n\nCompared to expand_colony_radom!, this function is faster for large images and many pixels to add,  but slower for small images and fewer pixels to add. This is due to the fact that the computational heavy convolution only needs to be      calculated once for the whole image, whereas the distance transform in expand_colony_radom needs to be calculated for each iteration of the loop.\n\n\n\n\n\n","category":"method"},{"location":"#ColonyImages.expand_colony_un_radom_cov!-Tuple{AbstractArray, Int64, Vector{Int64}}","page":"Home","title":"ColonyImages.expand_colony_un_radom_cov!","text":"expand_colony_un_radom_cov!(img::AbstractArray, pixels_to_add::Int, dir::Vector{Int}; still_spawn_rate::AbstractFloat = 0.5)\n\nExpand the colony in the image img by adding pixels_to_add pixels. The expansion is done randomly  at the border of the colony, but it is more likely to expand in the direction specified by dir.  The border is determined by convolving the image with a Laplacian kernel and finding points where the convolution  is greater than 0.1. The function modifies the input image in-place.\n\nArguments\n\nimg::AbstractArray: A 2D array representing the image of the colony. The colony is represented by 1s and the background by 0s.\npixels_to_add::Int: The number of pixels to add to the colony.\ndir::Vector{Int}: A vector representing the preferred direction of expansion.\nstill_spawn_rate::AbstractFloat: A float representing the probability of expanding in a direction opposite to dir. Default is 0.5.\n\nExample\n\nimg = zeros(100, 100)\nimg[50:55, 50:55] .= 1\nexpand_colony_un_radom_cov!(img, 100, [1, 0])\n\n\n\n\n\n","category":"method"},{"location":"#ColonyImages.expand_matrix-Tuple{Union{BitMatrix, Matrix{<:Real}}}","page":"Home","title":"ColonyImages.expand_matrix","text":"expand_matrix(mat::Union{Matrix{<:Real}, BitMatrix}; annuli::Int = 2 )\n\nExpands a matrix by adding zero-filled columns to its outermost right column and rotating it 90 degrees counterclockwise. This process is repeated 4 times per annulus. The number of additional bands (annuli) added to the matrix is determined by the annuli parameter.\n\nArguments mat::Union{Matrix{<:Real}, BitMatrix}: The input matrix. annuli::Int: The number of additional bands to be added to the matrix. Default is 2. Returns The expanded matrix. Examples\n\nmat = [1 2; 3 4]\nexpand_matrix(mat, annuli = 1)\n4×4 Matrix{Int64}:\n 0  0  0  0\n 0  1  2  0\n 0  3  4  0\n 0  0  0  0\n\n\n\n\n\n","category":"method"},{"location":"#ColonyImages.fill_holes-Tuple{Any, Real}","page":"Home","title":"ColonyImages.fill_holes","text":"fill_holes(img, size_holes::Real)\n\nFills holes in a binary image. The size of the holes to be filled is determined by the size_holes parameter.\n\nArguments\n\n`img: The input binary image.\nsize_holes::Real: The relative size of the holes to be filled. This is a fraction of the total number of pixels in the image.\n\nReturns\n\nA binary image with the holes filled.\n\n\n\n\n\n","category":"method"},{"location":"#ColonyImages.filter_fourier_alpha-Tuple{Vector{<:Real}}","page":"Home","title":"ColonyImages.filter_fourier_alpha","text":"filter_fourier_alpha(vec; a = 5)\n\nFilters a vector by setting all elements after the a-th element to zero. This function is useful for filtering the results of a Fourier transform, where the elements of the vector represent the amplitudes of the frequencies, and the index of the element represents the frequency. By setting all elements after the a-th element to zero, we effectively remove all frequencies higher than a.\n\nArguments\n\nvec::Vector{<:Real}: The input vector.\na::Int: The cutoff frequency. All frequencies higher than a are removed. Default is 5.\n\nReturns\n\nA new vector where all elements after the a-th element are zero.\n\n\n\n\n\n","category":"method"},{"location":"#ColonyImages.filter_fourier_beta-Tuple{Vector{<:Real}}","page":"Home","title":"ColonyImages.filter_fourier_beta","text":"filter_fourier_beta(vec::Vector{<:Real}; b::AbstractFloat = 0.5)\n\nFilters a vector by setting all elements whose absolute value is less than or equal to b times the maximum absolute value of the elements to zero. This function is useful for filtering the results of a Fourier transform, where the elements of the vector represent the amplitudes of the frequencies. By setting all elements whose amplitude is less than or equal to b times the maximum amplitude to zero, we effectively remove all frequencies with low amplitudes.\n\nArguments\n\nvec::Vector{<:Real}: The input vector.\nb::AbstractFloat: The threshold for the amplitude. All frequencies with amplitudes less than or equal to b times the maximum amplitude are removed. Default is 0.5.\n\nReturns\n\nA new vector where all elements whose absolute value is less than or equal to b times the maximum absolute value are zero.\n\n\n\n\n\n","category":"method"},{"location":"#ColonyImages.find_freq-Tuple{Vector{<:Real}}","page":"Home","title":"ColonyImages.find_freq","text":"find_freq(vec::Vector{<:Real}; ignore_latter_half = true )\n\nFinds the frequencies in a vector that have non-zero amplitudes. This function is useful for analyzing the results of a Fourier transform, where the elements of the vector represent the amplitudes of the frequencies, and the index of the element represents the frequency.\n\nArguments\n\nvec::Vector{<:Real}: The input vector.\nignore_latter_half::Bool: If true, the function only considers the first half of the vector. This is useful when the input vector is the result of a Fourier transform, where the second half of the vector contains the same information as the first half but in reverse order. Default is true.\n\nReturns\n\nA vector of the frequencies that have non-zero amplitudes.\n\n\n\n\n\n","category":"method"},{"location":"#ColonyImages.find_peaks-Tuple{Vector{<:Real}}","page":"Home","title":"ColonyImages.find_peaks","text":"find_peaks(signal::Vector{<:Real}; threshold::AbstractFloat = 1.0)\n\nFinds the local maxima inside a 1D signal, in areas where the signal exceeds its mean value by a given factor. \n\nArguments\n\nsignal::Vector{<:Real}: The input 1D signal.\nthreshold::AbstractFloat: The factor by which the signal needs to exceed its mean value to be considered a peak. Default is 1.0.\n\nReturns\n\nposition_peaks: A vector of the positions of the peaks in the signal.\nnr_peaks: The number of peaks found in the signal.\n\n\n\n\n\n","category":"method"},{"location":"#ColonyImages.generate_dir_vec-Tuple{Any, Any}","page":"Home","title":"ColonyImages.generate_dir_vec","text":"generate_dir_vec(number_fingers, rand_dist)\n\nGenerates vectors of x and y coordinates that span vectors ranging from 0 to 2π.  This function is used to generate the vectors used in expand_colony_finger_radom_cov!.\n\nArguments\n\nnumber_fingers: The number of vectors to generate.\nrand_dist: The random distance to add to the vectors.\n\nReturns\n\ndir: A vector of vectors, each containing the y and x coordinates of a vector.\n\nExample\n\nusing CairoMakie\nnumber_vecs = 20\ny,x = generate_dir_vec(20,0.1)\nyy = zeros(number_vecs)\narrows(yy,yy, y, x)\n\n\n\n\n\n","category":"method"},{"location":"#ColonyImages.lattice_points-Tuple{Int64}","page":"Home","title":"ColonyImages.lattice_points","text":"lattice_points(r::Int)\n\nGenerates a lattice of points within a circle of radius r.\n\nThe function returns a nested vector of points, each represented by its x and y coordinates.  The points are sorted by their distance from the origin and grouped into bands,  each band containing points that have a similar distance to the origin.\n\nArguments\n\nr::Int: The radius of the circle within which the lattice points are generated.\n\nReturns\n\npoints3::Vector{Vector{Vector{Int}}}: A nested vector of points.\n\npoints3 is a vector of vectors of Lattice vectors in ascending order by their length . Whereas the first entry of Points contains a vector of all lattice vectors which lengths are lower than 2:\n\nPoints[1] = [[0, 1], [0, 1], [0, -1], [0, -1], [1, 0], [-1, 0], [1, 0], [-1, 0], [1, 1], [-1, 1], [1, -1], [-1, -1]]\n\nPoints[2] does the same for length lower than 3:\n\nPoints[2] = [[0, 2], [0, 2], [0, -2], [0, -2], [2, 0], [-2, 0], [2, 0], [-2, 0], [1, 2], [-1, 2], [1, -2], [-1, -2], [2, 1], [-2, 1], [2, -1], [-2, -1], [2, 2], [-2, 2], [2, -2], [-2, -2]]\n\nand so on...\n\n\n\n\n\n","category":"method"},{"location":"#ColonyImages.occupied_points-Tuple{Union{BitMatrix, Matrix{<:Real}}}","page":"Home","title":"ColonyImages.occupied_points","text":"occupied_points(img::Union{Matrix{<:Real}, BitMatrix})\n\nCalculates the proportion of occupied points in a given binary image img.\n\nThe function sums up all the pixel values in the image and divides by the total number of pixels.  This gives the proportion of occupied points, assuming that non-zero pixel values represent occupied points.\n\nArguments\n\nimg::Union{Matrix{<:Real}, BitMatrix}: A 2D array representing the image. Non-zero values are considered as occupied points.\n\nReturns\n\nA float representing the proportion of occupied points in the image.\n\n\n\n\n\n","category":"method"},{"location":"#ColonyImages.pair_cor_metric3-Tuple{Union{BitMatrix, Matrix{<:Real}}, Vector{Int64}}","page":"Home","title":"ColonyImages.pair_cor_metric3","text":"pair_cor_metric3(img::Union{Matrix{<:Real}, BitMatrix}, center::Vector{Int}; samples::Int = 10000, steps::Int = 360)\n\nCalculates a pair correlation metric for a given image. The metric is a vector where each element represents the number of pairs of pixels that have a certain relative angle.  The angles are determined by dividing a half-circle into a certain number of equal parts.\n\nArguments\n\nimg::Union{Matrix{<:Real}, BitMatrix}: The input image.\ncenter::Vector{Int}: The reference point for calculating the relative angles.\nsamples::Int: The number of pairs of pixels to sample. Default is 10000.\nsteps::Int: The number of sectors into which the half-circle is divided. Default is 360.\n\nReturns\n\nA vector where each element represents the number of pairs of pixels that have a certain relative angle.\n\n\n\n\n\n","category":"method"},{"location":"#ColonyImages.rad2deg_discrete-Tuple{AbstractFloat}","page":"Home","title":"ColonyImages.rad2deg_discrete","text":"rad2deg_discrete(ϕ::AbstractFloat; steps::Int =360)\n\nConverts an angle in radians to a discrete angle in degrees. The output is the number of a circular sector on a unit circle divided into a given number of sectors.  For example, if the angle is 0.0 and the number of steps is 360, the output is 1, which corresponds to the first circular sector on the unit circle.  If the angle is (2pi - 0.01) and the number of steps is 360, the output is 360, which corresponds to the last circular sector on the unit circle.\n\nArguments\n\nϕ::AbstractFloat: The input angle in radians.\nsteps::Int: The number of circular sectors into which the unit circle is divided. Default is 360.\n\nReturns\n\nThe number of the circular sector to which the angle corresponds.\n\n\n\n\n\n","category":"method"},{"location":"#ColonyImages.res_scaling-Tuple{Any}","page":"Home","title":"ColonyImages.res_scaling","text":"res_scaling(img_int_vec; factor = 3, plots = 1)\n\nScales the resolution of a plot based on a given factor. The function counts the number of images in the given image_vec and and scales the resolution  of plot containg all these images accordingly.\n\nArguments\n\nimg_int_vec: A vector of images.\nfactor: The scaling factor. Defaults to 3.\nplots: The number of plots per images. Defaults to 1.\n\nReturns\n\nA tuple containing the scaled width and height of the image.\n\n\n\n\n\n","category":"method"}]
}
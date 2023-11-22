"""
    centroid(img::Union{Matrix{<:Real}, BitMatrix})

Calculates the centroid of a given image `img`.

The centroid is calculated as the average of the x and y coordinates of all non-zero pixels in the image, 
weighted by their intensity. The coordinates are then rounded to the nearest integer.

# Arguments
- `img::Union{Matrix{<:Real}, BitMatrix}`: A 2D array representing the image. Non-zero values are considered as part of the object to find the centroid of.

# Returns
- `centroid_norm::Vector{Int64}`: A vector containing the x and y coordinates of the centroid.

"""
function centroid(img::Union{Matrix{<:Real}, BitMatrix})
    # Get the size of the image
    Y,X = size(img)

    # Initialize the centroid and counter
    centroid = Float64[0.0,0.0]
    counter = 0

    # Iterate over all pixels in the image
    for x in 1:X
        for y in 1:Y
            # If the pixel is part of the object (non-zero), update the centroid and counter
            if img[y,x] != 0
                counter += img[y,x]
                centroid += img[y,x].*[y,x]
            end
        end
    end

    # Normalize the centroid by the total intensity of the object
    centroid_norm = round.(Int64,centroid./counter)

    return centroid_norm
end


"""
    lattice_points(r::Int)

Generates a lattice of points within a circle of radius `r`.

The function returns a nested vector of points, each represented by its x and y coordinates. 
The points are sorted by their distance from the origin and grouped into bands, 
each band containing points that have a similar distance to the origin.

# Arguments
- `r::Int`: The radius of the circle within which the lattice points are generated.

# Returns
- `points3::Vector{Vector{Vector{Int}}}`: A nested vector of points.

`points3` is a vector of vectors of `Lattice vectors` in ascending order by their length . Whereas the first entry of `Points` contains a vector of all lattice vectors which lengths are lower than 2:
```julia
Points[1] = [[0, 1], [0, 1], [0, -1], [0, -1], [1, 0], [-1, 0], [1, 0], [-1, 0], [1, 1], [-1, 1], [1, -1], [-1, -1]]
```
`Points[2]` does the same for length lower than 3:
```julia
Points[2] = [[0, 2], [0, 2], [0, -2], [0, -2], [2, 0], [-2, 0], [2, 0], [-2, 0], [1, 2], [-1, 2], [1, -2], [-1, -2], [2, 1], [-2, 1], [2, -1], [-2, -1], [2, 2], [-2, 2], [2, -2], [-2, -2]]
```
and so on...

"""
function lattice_points(r::Int)
    # Initialize the vector of points
    points = Vector{Vector{Real}}(undef,0)

    # Generate the points within the circle
    for x in 0:r , y in 0:r
        if !(x == 0 && y== 0)
            if sqrt(x^2 +y^2)< r
                # Add the point and its reflections in the other quadrants
                push!(points,[x,y,sqrt(x^2 +y^2)],[-x,y,sqrt(x^2 +y^2)],[-x,-y,sqrt(x^2 +y^2)],[x,-y,sqrt(x^2 +y^2)])
            end
        end
    end

    # Sort the points by their distance from the origin
    sort!(points, by = x-> x[3])

    # Round the coordinates of the points to the nearest integer
    points_r = [round.(Int,x[1:2]) for x in points]
    
    # Initialize the vector of grouped points
    points3 = Vector{Vector{Vector{Int}}}(undef,0)
    j = 1 
    val = 2
    
    # Group the points by their distance from the origin
    for i in 1:length(points)
        if points[i][3] >= val
            append!(points3,[points_r[j:i-1]])
            j = i
            val += 1
        end
    end

    return points3
end



"""
    occupied_points(img::Union{Matrix{<:Real}, BitMatrix})

Calculates the proportion of occupied points in a given binary image `img`.

The function sums up all the pixel values in the image and divides by the total number of pixels. 
This gives the proportion of occupied points, assuming that non-zero pixel values represent occupied points.

# Arguments
- `img::Union{Matrix{<:Real}, BitMatrix}`: A 2D array representing the image. Non-zero values are considered as occupied points.

# Returns
- A float representing the proportion of occupied points in the image.

"""
function occupied_points(img::Union{Matrix{<:Real}, BitMatrix})
    # Calculate the total number of pixels in the image
    total_pixels = size(img)[1]*size(img)[2]

    # Calculate the sum of all pixel values
    total_value = sum(img)

    # Calculate and return the proportion of occupied points
    return total_value / total_pixels
end


"""
    approx_radi_colo(img::Union{Matrix{<:Real}, BitMatrix})

Calculates the approximate diameter of a colony by summing up all the pixel values and taking the square root of the sum.

This function assumes that the pixel values represent the area of the colony. The diameter is then approximated using the formula for the diameter of a circle given its area.

# Arguments
- `img::Union{Matrix{<:Real}, BitMatrix}`: A 2D array representing the image. The pixel values are assumed to represent the area of the colony.

# Returns
- A float representing the approximate diameter of the colony.

"""
function approx_radi_colo(img::Union{Matrix{<:Real}, BitMatrix})
    return sum(img)^(1/2)
end
    


"""
    create_kernel(rad::Int;  geometry::String = "circle")

Creates a binary image kernel with a given radius. The kernel can be either a circle or a square.

# Arguments
- `rad::Int`: The radius of the kernel.
- `geometry::String`: The shape of the kernel. Can be either "circle" or "square". Defaults to "circle".

# Returns
- `kernel`: A 2D array representing the kernel.

"""
function create_kernel(rad::Int;  geometry::String = "circle")
    # Initialize the kernel with zeros
    kernel = zeros( 2*rad + 1, 2*rad + 1 )

    if geometry == "circle"
        # If the geometry is a circle, set the values inside the circle to 1
        for x in 1:2*rad+1, y in 1:2*rad+1
            dist = ( y - rad + 1 )^2 + ( x - rad + 1 )^2
            kernel[ y, x ] = ( rad^2 >= dist )
        end
    elseif geometry == "square"
        # If the geometry is a square, set all values to 1
        kernel .= 1
    end
    
    return kernel 
end


"""
    build_circle(center::Vector{Int}, img::Union{Matrix{<:Real}, BitMatrix}, points::Vector{Vector{Vector{Int}}}; threshold = 0.8::Float64)

Creates a binary image with the same size as the input image. The binary image is a circle with a given center. 
The circle is built by iterating over a set of points and setting the corresponding pixel in the binary image to 1 if the point is within the circle.
The occupation in the outermost circle band is calculated in each iteration and stored in the occupation vector. 
The function stops building the circle when the mean of the occupation vector is less than a given threshold.

# Arguments
- `center::Vector{Int}`: The center of the circle.
- `img::Union{Matrix{<:Real}, BitMatrix}`: The input image.
- `points::Vector{Vector{Vector{Int}}}`: A set of points used to build the circle.
- `threshold::Float64`: The threshold for the mean of the occupation vector. Defaults to 0.8.

# Returns
- `circle_kernel`: A binary image representing the circle.

"""
function build_circle(center::Vector{Int}, img::Union{Matrix{<:Real}, BitMatrix}, points::Vector{Vector{Vector{Int}}}; threshold = 0.8::Float64)
    # Initialize the binary image
    circle_kernel = zeros(Int, size(img))
    r = 1

    # Iterate until the radius is less than the minimum size of the image
    while r < minimum(size(img))
        # Initialize the occupation vector
        occupation_vec = zeros(Int,length(points[r]))

        # Iterate over the points
        for (i,point) in enumerate(points[r])
            # Calculate the coordinates of the point in the image
            point_c = point .+ center

            # Set the corresponding pixel in the binary image to 1
            circle_kernel[point_c...] = 1

            # If the corresponding pixel in the input image is 1, set the corresponding element in the occupation vector to 1
            if img[point_c...] == 1
                occupation_vec[i] = 1
            end
        end

        # Increase the radius
        r += 1

        # If the radius is greater than 10 and the mean of the occupation vector is less than the threshold, stop building the circle
        if r > 10 && mean(occupation_vec) < threshold
            break
        end
    end

    # Return the binary image
    return circle_kernel 
end




"""
    build_artifical_colony!(center::Vector{Int}, img::Union{Matrix{<:Real}, BitMatrix}, radius::Int, points::Vector{Vector{Vector{Int}}})

This function constructs an artificial spherical colony within a given image. The colony is represented as a circle with a specified center and radius. 
The function directly modifies the input image.

# Arguments
- `center::Vector{Int}`: A vector representing the center coordinates of the colony.
- `img`: The input image where the colony will be built.
- `radius::Int`: The radius of the colony.
- `points::Vector{Vector{Vector{Int}}}`: A nested vector containing the points used to construct the colony.

# Returns
- The image with the built colony.

"""
function build_artifical_colony!(center::Vector{Int}, img::Union{Matrix{<:Real}, BitMatrix}, radius::Int, points::Vector{Vector{Vector{Int}}})
    # Initialize the radius
    r = 0

    # Loop until the radius is less than or equal to the specified radius
    while r <= radius
        # If the radius is 0, set the pixel at the center to 1
        if r == 0
            img[center...] = 1
        else
            # If the radius is not 0, iterate over the points at the current radius
            for (i,point) in enumerate(points[r])
                # Calculate the actual coordinates of the point in the image
                point_c = point .+ center

                # Set the pixel at the calculated coordinates to 1
                img[point_c...] = 1
            end
        end

        # Increment the radius
        r += 1
    end

    # Return the modified image
    return img
end

"""
    expand_colony_circular!(img::Union{Matrix{<:Real}, BitMatrix}, points::Vector{Vector{Vector{Int}}}, center::Vector{Int}, pixels_to_add::Int)

This function expands a colony in a circular pattern within an image. The expansion starts from the center of the colony and proceeds outward. 
The function directly modifies the input image.

# Arguments
- `img`: The input image where the colony will be expanded.
- `points::Vector{Vector{Vector{Int}}}`: A nested vector containing the points used to expand the colony.
- `center::Vector{Int}`: A vector representing the center coordinates of the colony.
- `pixels_to_add::Int`: The number of pixels to add to the colony.

"""
function expand_colony_circular!(img::Union{Matrix{<:Real}, BitMatrix}, points::Vector{Vector{Vector{Int}}}, center::Vector{Int}, pixels_to_add::Int)
    # Initialize the pixel count
    pix_count = 0

    # Loop over the radius from 1 to the minimum dimension of the image
    for r in 1:(round(Int, minimum(size(img))))
        # For each radius, iterate over the points
        for (i, point) in enumerate(points[r])
            # Calculate the actual coordinates of the point in the image
            point_c = point .+ center

            # If the pixel at the calculated coordinates is not black yet, make it black
            if img[point_c...] == 0
                img[point_c...] = 1
                pix_count += 1
            end

            # If the pixel count exceeds the number of pixels to add, exit the loop
            if pix_count > pixels_to_add
                @goto escape
            end
        end
    end

    @label escape
end

"""
    expand_colony_radom!(img::Union{Matrix{<:Real}, BitMatrix}, pixels_to_add::Int)

This function expands a colony in a random pattern within an image. The expansion is performed by adding pixels to the border of the colony. 
The function directly modifies the input image.

# Arguments
- `img`: The input image where the colony will be expanded.
- `pixels_to_add::Int`: The number of pixels to add to the colony.

"""
function expand_colony_radom!(img::Union{Matrix{<:Real}, BitMatrix}, pixels_to_add::Int)
    # Initialize the pixel count
    pix_count = 0

    # Find the border points of the colony and shuffle them
    border_points = shuffle!(findall((distance_transform(feature_transform(img)).== 1)))
    shuffle_counter = 0

    # Loop until the desired number of pixels have been added
    while pix_count < pixels_to_add
        # For each point in the shuffled border points
        for point in border_points
            # If the pixel at the point is not black yet, make it black
            if img[point] == 0
                img[point] = 1
                pix_count += 1
            end

            # Increment the shuffle counter
            shuffle_counter += 1
            
            # If the shuffle counter reaches 10% of the length of the border points, break the loop
            if shuffle_counter >= length(border_points)*0.1
                break
            end            
        end

        # Find the new border points of the colony and shuffle them
        border_points = shuffle!(findall((distance_transform(feature_transform(img)).== 1)))
        shuffle_counter = 0 
    end
end
                

function expand_colony_radom_cov!(img,pixels_to_add)
    pix_count = 0
    laplac_kernel = [0 1 0; 1 -4 1; 0 1 0]
    cov_img = conv( img, laplac_kernel )
    border_points = shuffle!(findall(cov_img .> 0.1))
    shuffle_counter = 0
    while pix_count < pixels_to_add
        for point in border_points
            if img[point] == 0
                img[point] = 1
                #also modify cov image, as change in cov_image is just convolution of kernel with added point
                # which is only the kernel itself as a single added point is the identity element of convolution 
                cov_img[point[1]-1:point[1]+1,point[2]-1:point[2]+1] += laplac_kernel
                pix_count += 1
            end
            shuffle_counter += 1
            
            if shuffle_counter >= length(border_points)*0.1
                break
            end            
        end
        border_points = shuffle!(findall(cov_img .> 0.1))
        shuffle_counter = 0 
    end
end


""" Expands colony randomly but it is more likely to expand to points that only have 2 neigboors or less
"""
function expand_colony_fractal_radom_cov!(img,pixels_to_add; max_neigbours = 1, prob_neigbour = 0.3)
    pix_count = 0
    laplac_kernel = [0 1 0; 1 -4 1; 0 1 0]
    cov_img = conv( img, laplac_kernel )
    border_points = shuffle!(findall(cov_img .> 0.1))
    shuffle_counter = 0
    while pix_count < pixels_to_add
        for point in border_points
            
            if cov_img[point] > max_neigbours+ 0.1
                if rand() > prob_neigbour
                    shuffle_counter += 1
            
                    if shuffle_counter >= length(border_points)*0.1
                        break
                    end    
                    continue
                end
            end
            
            if img[point] == 0
                img[point] = 1
                #also modify cov image, as change in cov_image is just convolution of kernel with added point
                # which is only the kernel itself as a single added point is the identity element of convolution 
                cov_img[point[1]-1:point[1]+1,point[2]-1:point[2]+1] += laplac_kernel
                pix_count += 1
            end
            shuffle_counter += 1
            
            if shuffle_counter >= length(border_points)*0.1
                break
            end            
        end
        border_points = shuffle!(findall(cov_img .> 0.1))
        shuffle_counter = 0 
    end
end

function expand_colony_un_radom_cov!(img,pixels_to_add,dir)
    pix_count = 0
    laplac_kernel = [0 1 0; 1 -4 1; 0 1 0]
    cov_img = conv( img, laplac_kernel )
    border_points = shuffle!(findall(cov_img .> 0.1))
    shuffle_counter = 0
    center = size(img).÷2
    
    
    while pix_count < pixels_to_add
        for point in border_points
            # if dot product between preferred dir and chosen point is negative only proceed in 50% of the time
            if dot(dir,point.I .- center) .< 0 && rand() > 0.5
                continue
            else
                if img[point] == 0
                    img[point] = 1
                    #also modify cov image, as change in cov_image is just convolution of kernel 
                    #with added point which is only the kernel itself as a single
                    #added point is the identity element of convolution 
                    cov_img[point[1]-1:point[1]+1,point[2]-1:point[2]+1] += laplac_kernel
                    pix_count += 1
                end
                shuffle_counter += 1

                if shuffle_counter >= length(border_points)*0.1
                    break
                end     
            end
        end
        border_points = shuffle!(findall(cov_img .> 0.1))
        shuffle_counter = 0 
    end
end
        


function expand_colony_finger_radom_cov!(img,pixels_to_add,dir; dir_match_rate = 0.999, still_spawn_rate = 0.99, min_neigbour = 2 )
    pix_count = 0
    laplac_kernel = [0 1 0; 1 -4 1; 0 1 0]
    cov_img = conv( img, laplac_kernel )
    border_points = shuffle!(findall(cov_img .> 0.1))
    shuffle_counter = 0
    center = size(img).÷2
    
    
    while pix_count < pixels_to_add
        for point in border_points
            # two checks are performed to determine if a randomly chosen border points is filled
            # 1. It is checked if the dot produced of the points_dir (relative to center of colony) and of
            # the randomly generated vectors of dir is bigger than the dir_match_rate
            #2. If that's not the case there is a small percentage (still_spawn_rate) that the point is filled nevertheless
            # but only if the points neigbours at least 2 occupied points already
            norm_dir = normalize(collect(point.I .- center))
            if 1 ∉ (dot.(dir,[norm_dir for i in length(dir)]) .>dir_match_rate)
                if rand() < still_spawn_rate
                    
                    # check if more than x neighbour point are occupied
                    if cov_img[point] < (min_neigbour- 0.1)
                    
                        
                        shuffle_counter += 1

                    
                        continue
                    end
                else 
                     shuffle_counter += 1

                    
                        continue
                end
            end
                if img[point] == 0
                    img[point] = 1
                    #also modify cov image, as change in cov_image is just convolution of kernel 
                    #with added point which is only the kernel itself as a single
                    #added point is the identity element of convolution 
                    cov_img[point[1]-1:point[1]+1,point[2]-1:point[2]+1] += laplac_kernel
                    pix_count += 1
                end
                shuffle_counter += 1

                if shuffle_counter >= length(border_points)*0.1
                    break
                end     
            
        end
        border_points = shuffle!(findall(cov_img .> 0.1))
        shuffle_counter = 0 
    end
end
       



"""
    generate_dir_vec(number_fingers, rand_dist)

Generates vectors of x and y coordinates that span vectors ranging from 0 to 2π.

# Arguments
- `number_fingers`: The number of vectors to generate.
- `rand_dist`: The random distance to add to the vectors.

# Returns
- `dir`: A vector of vectors, each containing the y and x coordinates of a vector.

# Example
```julia
using CairoMakie
number_vecs = 20
y,x = generate_dir_vec(20,0.1)
yy = zeros(number_vecs)
arrows(yy,yy, y, x)
```
"""
function generate_dir_vec(number_fingers, rand_dist)
    # Generate the angles for the vectors vectus
    vectus = 2π/number_fingers.*[1:number_fingers...]+ rand(number_fingers) *rand_dist
    # Calculate the y and x coordinates of the vectors
    y = sin.(vectus)
    x = cos.(vectus)

    # Combine the y and x coordinates into a vector of vectors
    dir = [[y[i],x[i]] for i in 1:length(y)]

return dir
end
    

"""
    conv(img::Union{Matrix{<:Real}, BitMatrix}, kernel::Union{Matrix{<:Real}, BitMatrix})

Performs a convolution operation on an image using a given kernel. The input image and kernel are 2D arrays of Int or Float or Bool. 
The function returns a 2D Float64 array of the same size as the input image.

# Arguments
- `img::Union{Matrix{<:Real}, BitMatrix}`: The input image, a 2D array of Int, Float or Bool.
- `kernel::Union{Matrix{<:Real}, BitMatrix}`: The kernel used for the convolution, a smaller 2D array of Int,Float or Bool.

# Returns
- A 2D Float64 array representing the convolved image.

"""
function conv(img::Union{Matrix{<:Real}, BitMatrix}, kernel::Union{Matrix{<:Real}, BitMatrix})
    # Calculate the size of the convolution
    csize = size(img) .+ size(kernel) .- 1

    # Initialize padded versions of the image and kernel
    padimg = zeros(ComplexF32, csize)
    padkernel = zeros(ComplexF32, csize)

    # Copy the image and kernel into the padded versions
    padimg[1:size(img,1), 1:size(img,2)] .= img
    padkernel[1:size(kernel,1), 1:size(kernel,2)] .= kernel

    # Perform the Fast Fourier Transform on the padded image and kernel
    fft!(padimg)
    fft!(padkernel)

    # Multiply the transformed image and kernel
    padimg .*= padkernel

    # Perform the inverse Fast Fourier Transform on the result
    ifft!(padimg)

    # Extract the real part of the result
    output = real.(padimg)

    # Calculate the offset for the kernel
    off = div.(size(kernel), 2)

    # Extract the convolved image from the output
    h, w = size(img)
    return output[1+off[1]:h+off[1], 1+off[2]:w+off[2]]
end

b_w_n(img) = (round.(Int64,Float64.(Gray.(img))).-1).*-1
b_w_i(img )= round.(Int64,Float64.(Gray.(img)))


"""
    b_w(img)

Converts a grayscale colony image into a binary image/BitArray. If more than half of the image is black, it inverts the image. 
This ensures that in the output image, the pixels inside the colony are always set to 1 and the background pixels to 0, 
regardless of the inversion status of the input image.

# Arguments
- `img`: The input image.

# Returns
- A binary image where the colony pixels are set to 1 and the background pixels are set to 0.

"""
function b_w(img)
    # Convert the image to binary
    img = b_w_i(img)

    # If more than half of the image is black, invert the image
    if sum(img) >= prod(size(img))*0.5
        img = (img .-1) .*-1 
    end

    # Convert the image to a BitArray and return it
    return Bool.(img)
end


"""
    fill_holes(img, size_holes::Real)

Fills holes in a binary image. The size of the holes to be filled is determined by the `size_holes` parameter.

# Arguments
- `img: The input binary image.
- `size_holes::Real`: The relative size of the holes to be filled. This is a fraction of the total number of pixels in the image.

# Returns
- A binary image with the holes filled.

"""
function fill_holes(img, size_holes::Real)
    # Calculate the absolute size of the holes to be filled
    size_absolut = size(img)[1]*size(img)[2]*size_holes 

    # Invert the image, fill the holes, and then invert the image again
    return .!(imfill(.!(Bool.(img)),(0.1,size_absolut)))
end


"""
    rad2deg_discrete(ϕ::AbstractFloat; steps::Int =360)

Converts an angle in radians to a discrete angle in degrees. The output is the number of a circular sector on a unit circle divided into a given number of sectors. 
For example, if the angle is 0.0 and the number of steps is 360, the output is 1, which corresponds to the first circular sector on the unit circle. 
If the angle is (2pi - 0.01) and the number of steps is 360, the output is 360, which corresponds to the last circular sector on the unit circle.

# Arguments
- `ϕ::AbstractFloat`: The input angle in radians.
- `steps::Int`: The number of circular sectors into which the unit circle is divided. Default is 360.

# Returns
- The number of the circular sector to which the angle corresponds.

"""
function rad2deg_discrete(ϕ::AbstractFloat; steps::Int =360)
    # Calculate the size of each step in degrees
    stepsize = 360/steps 

    # Convert the angle to degrees and shift it by 180 degrees
    ϕ_a = rad2deg(ϕ) + 180

    # If the angle is 360 or more, reset it to 0
    if ϕ_a >= 360
        ϕ_a = 0.0
    end

    # Calculate the discrete angle
    ϕ_a_d = round(Int, (ϕ_a - 0.5) / stepsize) + 1

    # If the discrete angle is less than 1, set it to 1
    if ϕ_a_d < 1
        ϕ_a_d = 1 
    end

    # If the discrete angle is greater than or equal to the number of steps, set it to 1
    if ϕ_a_d >= steps
        ϕ_a_d = 1
    end

    # Return the discrete angle
    return ϕ_a_d
end



function angular_metric_shit(img, center; steps = 360)
    angle_vec = zeros(Int64, steps)
    for y in 1:size(img)[1], x in 1: size(img)[2]
        if img[y,x] != 0
            vec = [y,x] .- center
            ϕ = atan(vec...)
            step_index =rad2deg_discrete(ϕ, steps = steps)
            angle_vec[step_index] += 1
        end
    end
    #last and first circular segment are buggy due to discretions artefacts, mean of both of them
    mean_angle = round(Int,mean([angle_vec[1], angle_vec[end]]))
    angle_vec[1], angle_vec[end] = mean_angle,mean_angle
    return angle_vec
end

"""
    angular_metric(img::Union{Matrix{<:Real}, BitMatrix}, center::Vector{Int}; steps::Int = 360)

Calculates an angular metric for a given image. The metric is a vector where each element represents the number of pixels in a certain angular sector of the image. 
The sectors are determined by dividing a circle centered at a given point into a certain number of equal parts.

# Arguments
- `img::Union{Matrix{<:Real}, BitMatrix}`: The input image.
- `center::Vector{Int}`: The center of the circle.
- `steps::Int`: The number of sectors into which the circle is divided. Default is 360.

# Returns
- A vector where each element represents the number of pixels in a certain angular sector of the image.

"""
function angular_metric(img::Union{Matrix{<:Real}, BitMatrix}, center::Vector{Int}; steps::Int = 360)
    # Initialize the vector for the angular metric
    angle_vec = zeros(Int64, steps)

    # Calculate the size of each step in radians
    stepsize = 2π/steps 

    # For each pixel in the image
    for y in 1:size(img)[1], x in 1: size(img)[2]
        # If the pixel is not zero
        if img[y,x] != 0
            # Calculate the vector from the center to the pixel
            vec = [y,x] .- center

            # Calculate the angle of the vector
            ϕ = atan(vec...) + π

            # Determine the sector to which the angle belongs
            step_index = minimum([round(Int, ϕ / stepsize), steps - 1])

            # Increment the count for the sector
            angle_vec[step_index + 1] += 1
        end
    end

    # The first and last sectors may be inaccurate due to discretization artifacts, so take the mean of both
    mean_angle = round(Int, mean([angle_vec[1], angle_vec[end]]))
    angle_vec[1], angle_vec[end] = mean_angle, mean_angle

    # Return the angular metric
    return angle_vec
end


function pair_cor_metric(img, center_1; samples = 10000, steps = 360)
    angle_vec = zeros(Int64, steps)
    stepsize = π/steps
    y_size, x_size = size(img)
    i = 0
    while i < samples
    #for i in 1:samples
        v1 = [rand(1:y_size),rand(1:x_size)]
        v2 = [rand(1:y_size),rand(1:x_size)]
        if  img[v1...] != 0 && img[v2...] != 0
            v1 = v1 .- center_1
            v2 = v2 .- center_1
            x = dot(v1,v2)/(norm(v1)*norm(v2))
            angle = acos(x >1.0 ? 1.0 : x < -1.0 ? -1.0 : x)
            step_index = (minimum([round(Int64,replace_nan(angle/stepsize)), steps-1]))
            angle_vec[step_index+1] += 1
            i += 1
            #println(i)
        end
    end
    return angle_vec
end
    
function pair_cor_metric2(img, center; samples = 10000, steps = 360)
    angle_vec = zeros(Int64, steps+1)
    stepsize = π/steps
    indizies = Tuple.(findall(x-> x != 0 ,img))
    angles = Float64[]
    for i in 1:samples
        v1, v2  = sample(indizies, 2)
        v1_v = v1 .- center 
        v2_2 = v2 .- center
        x = dot(v1_v,v2_2)/(norm(v1_v)*norm(v2_2))
        angle = acos(x >1.0 ? 1.0 : x < -1.0 ? -1.0 : x)
        step_index = round(Int64,replace_nan(angle/stepsize))
        # sometimes step index can be bigger than 359 for unknown reasons, if so a random value between 0:359 is assigned
        #step_index_final = step_index > 359 ? rand(0:359) : step_index 
        angle_vec[step_index+1] += 1
           
       
    end
    return angle_vec
end

"""
    pair_cor_metric3(img::Union{Matrix{<:Real}, BitMatrix}, center::Vector{Int}; samples::Int = 10000, steps::Int = 360)

Calculates a pair correlation metric for a given image. The metric is a vector where each element represents the number of pairs of pixels that have a certain relative angle. 
The angles are determined by dividing a half-circle into a certain number of equal parts.

# Arguments
- `img::Union{Matrix{<:Real}, BitMatrix}`: The input image.
- `center::Vector{Int}`: The reference point for calculating the relative angles.
- `samples::Int`: The number of pairs of pixels to sample. Default is 10000.
- `steps::Int`: The number of sectors into which the half-circle is divided. Default is 360.

# Returns
- A vector where each element represents the number of pairs of pixels that have a certain relative angle.

"""
function pair_cor_metric3(img::Union{Matrix{<:Real}, BitMatrix}, center::Vector{Int}; samples::Int = 10000, steps::Int = 360,)
    # Initialize the vector for the pair correlation metric
    angle_vec = zeros(Int64, steps)

    # Calculate the size of each step in radians
    stepsize = π/steps

    # Find the indices of the non-zero pixels
    indizies = Tuple.(findall(x-> x != 0 ,img))

    # For each sample
    for i in 1:samples
        # Randomly select two non-zero pixels
        v1, v2  = sample(indizies, 2, replace = false)

        # If the reference point is one of the pixels, select two new pixels
        while Tuple(center) in [v1,v2]
            v1, v2  = sample(indizies, 2, replace = false)
        end

        # Calculate the vectors from the reference point to the pixels
        v1_v = v1 .- center 
        v2_v = v2 .- center

        # Calculate the cosine of the angle between the vectors
        x = dot(v1_v,v2_v)/(norm(v1_v)*norm(v2_v))

        # Calculate the angle between the vectors
        angle = acos(x >1.0 ? 1.0 : x < -1.0 ? -1.0 : x)

        # Determine the sector to which the angle belongs
        step_index = trunc(Int64,angle/stepsize)

        # If the step index is 360 due to rounding errors, set it to 359
        step_index_final = (step_index == 360 ? 359 : step_index)

        # Increment the count for the sector
        angle_vec[step_index_final+1] += 1
    end

    # Return the pair correlation metric
    return angle_vec
end


"""
    res_scaling(img_int_vec; factor = 3, plots = 1)

Scales the resolution of a plot based on a given factor.
The function counts the number of images in the given image_vec and and scales the resolution 
of plot containg all these images accordingly.

# Arguments
- `img_int_vec`: A vector of images.
- `factor`: The scaling factor. Defaults to 3.
- `plots`: The number of plots per images. Defaults to 1.

# Returns
- A tuple containing the scaled width and height of the image.

"""
function res_scaling(img_int_vec; factor = 3, plots = 1)
    # Initialize the counter
    c = 0

    # Iterate over the images
    for i in img_int_vec
        # Iterate over the slices in the image
        for j in 1:size(i,3)
            # Increment the counter
            c += 1
        end
    end

    # Calculate the scaled width and height of the image
    width = round(Int64, factor * 1000) * plots
    height = round(Int64, factor * 200 * (c ÷ 5 + ((c % 5 == 0) ? 0 : 1)))

    return width, height
end



replace_nan_1(x) = ismissing(x) || (x isa Number && isnan(x)) ? 1 : x

replace_nan(x) = ismissing(x) || (x isa Number && isnan(x)) ? 0 : x

function filter_fourier_alpha(vec::Vector{<:Real}; a = 5)
    vec2 = deepcopy(vec)
    for i in 1:length(vec2)
        if i > a
            vec2[i] = 0.0
        end
    end
    return vec2
end
    

function filter_fourier_beta(vec; b = 0.5)
    max = maximum(abs.(vec))
    vec2 = deepcopy(vec)
    for i in 1:length(vec2)
        if abs(vec2[i]) <= max*b
            vec2[i] = 0.0
    
        end 
    end
    return vec2
end

function find_freq(vec; ignore_latter_half = true )
    freq_vec = Int64[]
    len_lopp = length(vec)
    if ignore_latter_half
        len_lopp = len_lopp ÷2
    end
    
    ### change to 2 iteration, not sure if sensefull 
    
    for i in 2:len_lopp
        if vec[i] != 0.0
            push!(freq_vec,i)
        end
    end
    return freq_vec
end
        

function filter_fourier_beta2(vec; b = 0.5)
    max = maximum(abs.(vec))
    vec2 = deepcopy(vec)
    for i in 1:length(vec2)
        if abs(vec2[i]) <= max*b
                vec2[i] = 0.0
    
        end 
        
    end
    return vec2
end


"Finds the local maxima inside a 1D signal, in areas where the 
signal exceed an its mean value by a given factor"
function find_peaks(signal; threshold = 1.0)
    mean_s = mean(signal)*threshold
    nr_peaks = 0 
    inside_peak = 0
    position_peaks = typeof(signal)(undef,0)
    current_peak = typeof(signal)(undef,0)
    for (i,x) in enumerate(signal )
        if x > mean_s
            inside_peak = 1
            push!(current_peak, i)
        elseif x < mean_s && inside_peak == 1
            inside_peak = 0
            push!(position_peaks, findmax(signal[[current_peak...]])[2]+i-1-length(current_peak))
            nr_peaks += 1
            current_peak = typeof(signal)(undef,0)
        end
    end
    return position_peaks, nr_peaks
end


function expand_matrix(mat; annuli = 2 )
    if annuli ≈ 0
        return mat
    else
        mat = hcat(mat, zeros(eltype(mat),size(mat)[1]))
        mat = rotr90(mat)
    return expand_matrix(mat, annuli = annuli -0.25)
    end
end
"""
    initialize_colonies(para::parameters)

Initialize colonies for simulations.

# Arguments
- `para::parameters`: A `parameters` object containing simulation parameters.

# Returns
- `vec_of_sims`: A vector of vectors of `BitArray{3}` representing the initialized empty arrays which will later became colonies.

"""
function initialize_colonies(para::parameters)
    vec_of_sims = Vector{Vector{BitArray{3}}}(undef, 0)
    for i in 1:length(para.simulations)
        img_int_vec_1 = Vector{BitArray{3}}(undef, 0)
        for i in 1:para.colony_nr
            hans = BitArray(zeros(Bool, para.im_size..., para.stacks))
            push!(img_int_vec_1, hans)
        end
        push!(vec_of_sims, img_int_vec_1)
    end
    return vec_of_sims
end


"""
    build_artifical_colony!(center::Vector{Int}, img::AbstractArray}, radius::Int, points::Vector{Vector{Vector{Int}}})

This function constructs an artificial spherical colony within a given image. The colony is represented as a circle with a specified center and radius. 
The function directly modifies the input image.

# Arguments
- `center::Vector{Int}`: A vector representing the center coordinates of the colony.
- `img::AbstractArray`: The input image where the colony will be built.
- `radius::Int`: The radius of the colony.
- `points::Vector{Vector{Vector{Int}}}`: A nested vector containing the points used to construct the colony.

# Returns
- The image with the built colony.

"""
function build_artifical_colony!(center::Vector{Int}, img::AbstractArray, radius::Int, points::Vector{Vector{Vector{Int}}})
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
    expand_colony_circular!(img::AbstractArray, points::Vector{Vector{Vector{Int}}}, center::Vector{Int}, pixels_to_add::Int)

This function expands an artifical colony in a circular pattern within an image. The expansion starts from the center of the colony and proceeds outward. 
The function directly modifies the input image.

# Arguments
- `img::AbstractArray`: The input image where the colony will be expanded.
- `points::Vector{Vector{Vector{Int}}}`: A nested vector containing the points used to expand the colony.
- `center::Vector{Int}`: A vector representing the center coordinates of the colony.
- `pixels_to_add::Int`: The number of pixels to add to the colony.

"""
function expand_colony_circular!(img::AbstractArray, points::Vector{Vector{Vector{Int}}}, center::Vector{Int}, pixels_to_add::Int)
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
    expand_colony_radom_cov!(img::AbstractArray, pixels_to_add::Int)

Expand the colony in the image `img` by adding `pixels_to_add` pixels. The expansion is done randomly 
at the border of the colony. The border is determined by convolving the image with a Laplacian kernel 
and finding points where the convolution is greater than 0.1. The function modifies the input image in-place.

Compared to `expand_colony_radom!`, this function is faster for large images and many pixels to add, 
but slower for small images and fewer pixels to add. This is due to the fact that the computationally heavy convolution only needs to be 
calculated once for the whole image, whereas the distance transform in `expand_colony_radom` needs to be calculated for each iteration of the loop.

# Arguments
- `img::AbstractArray`: A 2D array representing the image of the colony. The colony is represented by 1s and the background by 0s.
- `pixels_to_add::Int`: The number of pixels to add to the colony.

# Example
```julia
img = zeros(100, 100)
img[50:55, 50:55] .= 1
expand_colony_radom_cov!(img, 100)
```
Compared to `expand_colony_radom!`, this function is faster for large images and many pixels to add, 
but slower for small images and fewer pixels to add. This is due to the fact that the computational heavy convolution only needs to be 
    calculated once for the whole image, whereas the distance transform in `expand_colony_radom` needs to be calculated for each iteration of the loop.
"""
function expand_colony_radom_cov!(img::AbstractArray,pixels_to_add::Int)
     # Initialize pixel count and Laplacian kernel 
    pix_count = 0 
    laplac_kernel = [0 1 0; 1 -4 1; 0 1 0]
    # Convolve the image with the Laplacian kernel
    cov_img = conv( img, laplac_kernel )

    # Find border points where the convolution is greater than 0.1 and shuffle them
    border_points = shuffle!(findall(cov_img .> 0.1))
    shuffle_counter = 0

    # Loop until the desired number of pixels have been added
    while pix_count < pixels_to_add
        point = rand(border_points)
        
        # If the point is in the background, add it to the colony
        if img[point] == 0
            img[point] = 1
            
            # Update the convolution image as well
            cov_img[point[1]-1:point[1]+1,point[2]-1:point[2]+1] += laplac_kernel
            pix_count += 1

            new_border_points = findall(cov_img[point[1]-1:point[1]+1,point[2]-1:point[2]+1] .> 0.1) .+( point - CartesianIndex(2,2))
            append!(border_points, new_border_points)
            # old border point stays in border point arrays, due to computational efficiency, but is never choosen anyways, as it does not pass the if statement
        end
        shuffle_counter += 1
        

        # recreate the border_points vector if  over 50% of the original border points have been checked
        if shuffle_counter >= length(border_points)*0.5
            border_points = shuffle!(findall(cov_img .> 0.1))
            shuffle_counter = 0 
        end            
        
    end
end




"""
    expand_colony_point!(img::AbstractArray, cov_img::AbstractArray, point::Vector{Int})

Expands a colony at a given point in an image. The function uses a Laplacian kernel to convolve the image and find border points. 
If the given point is in the border points and is in the background of the image, it is added to the colony and the convolution image is updated.

# Arguments
- `img::AbstractArray`: The input image.
- `cov_img::AbstractArray`: The convolution image.
- `point::Vector{Int}`: The point at which to expand the colony.

# Returns
- The updated convolution image.
"""
function expand_colony_point!(img::AbstractArray, cov_img::AbstractArray, point::CartesianIndex{2})
    # Define the Laplacian kernel
    laplac_kernel = [0 1 0; 1 -4 1; 0 1 0]

    # Find border points where the convolution is greater than 0.1
    border_points = (findall(cov_img .> 0.1))
    if point in border_points
        # If the point is in the background, add it to the colony
        if img[point] == 0
            img[point] = 1

            # Update the convolution image as well
            cov_img[point[1]-1:point[1]+1,point[2]-1:point[2]+1] += laplac_kernel
            # Note: pix_count is not defined in this function. It should be defined outside this function and passed as a reference if it's needed to be updated.
        end
    else println("Point is not in border points")
    end

    # Return the updated convolution image
    #return cov_img	
end



"""
    expand_colony_radom!(img::AbstractArray, pixels_to_add::Int)

This function expands a colony in a random pattern within an image. The expansion is performed by adding pixels to the border of the colony. 
The function directly modifies the input image. 

# Arguments
- `img::AbstractArray`: The input image where the colony will be expanded.
- `pixels_to_add::Int`: The number of pixels to add to the colony.

Compared to `expand_colony_radom_cov!`, this function is slower for large images and many pixels to add,
 but faster for small images and and fewer pixels to add .This is due to the fact that the computational heavy convolution only needs to be 
    calculated once for the whole image, whereas the distance transform in `expand_colony_radom` needs to be calculated for each iteration of the loop.
"""
function expand_colony_radom!(img::AbstractArray, pixels_to_add::Int)
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
    



"""
    expand_colony_un_radom_cov!(img::AbstractArray, pixels_to_add::Int, dir::Vector{Int}; still_spawn_rate::AbstractFloat = 0.5)

This function expands a colony in a given direction by adding pixels to the image.

# Arguments
- `img::AbstractArray`: A 2D array representing the image of the colony.
- `pixels_to_add::Int`: The number of pixels to add to the colony.
- `dir::Vector{Int}`: A vector representing the direction in which to expand the colony.
- `still_spawn_rate::AbstractFloat`: The probability of adding a pixel even if it is not in the desired direction. Default is 0.5.

# Example
```julia
img = zeros(100, 100)
img[50:55, 50:55] .= 1
expand_colony_un_radom_cov!(img, 100, [1, 0])
```
"""
function expand_colony_un_random_cov!(img::AbstractArray, pixels_to_add::Int, dir::Vector{Int}; still_spawn_rate::AbstractFloat = 0.5)
    pix_count = 0  # Initialize pixel count
    laplac_kernel = [0 1 0; 1 -4 1; 0 1 0]  # Define Laplacian kernel
    cov_img = conv( img, laplac_kernel )  # Calculate Laplacian of image
    border_points = shuffle!(findall(cov_img .> 0.1))  # Find border points and shuffle them
    shuffle_counter = 0  # Initialize shuffle counter
    center = size(img).÷2  # Calculate center of image
    
    while pix_count < pixels_to_add  # Loop until desired number of pixels have been added
        point = rand(border_points)  # Randomly select a border point
        
        # If the point is not in the desired direction and a random number is greater than still_spawn_rate, skip to the next iteration
        if dot(dir,point.I .- center) .< 0 && rand() > still_spawn_rate
            continue
        else
            # If the point is in the desired direction or a random number is less than still_spawn_rate, add a pixel at that point and update the border
            if img[point] == 0
                img[point] = 1
                cov_img[point[1]-1:point[1]+1,point[2]-1:point[2]+1] += laplac_kernel
                pix_count += 1
                new_border_points = findall(cov_img[point[1]-1:point[1]+1,point[2]-1:point[2]+1] .> 0.1) .+( point - CartesianIndex(2,2))
                append!(border_points, new_border_points)
            end
            shuffle_counter += 1
        end
        # Reshuffle the border points after every half of the border points have been considered
        if shuffle_counter >= length(border_points)*0.1
            border_points = shuffle!(findall(cov_img .> 0.1))
            shuffle_counter = 0 
        end
    end
end

"""
    generate_dir_vec(para::parameters)

Generates vectors of x and y coordinates that span vectors ranging from 0 to 2π. 
This function is used to generate the vectors used in `expand_colony_finger_radom_cov!`.

# Arguments
- `para::parameters`: A parameters object containing various parameters for the analysis. 
  - `number_finger`: The number of vectors to generate.
  - `finger_dist`: The random distance to add to the vectors.

# Returns
- `dir::Array{Array{Float64,1},1}`: A vector of vectors, each containing the y and x coordinates of a vector.

# Example
```julia
using CairoMakie
para = parameters(number_finger = 20, finger_dist = 0.1)
dir = generate_dir_vec(para)
yy = zeros(para.number_finger)
arrows(yy,yy, [d[1] for d in dir], [d[2] for d in dir])
```
"""
function generate_dir_vec(para::parameters)
    # Generate the angles for the vectors vectus
    vectus = 2π/para.number_finger.*[1:para.number_finger...].+ rand(para.number_finger) *para.finger_dist
    # Calculate the y and x coordinates of the vectors
    y = sin.(vectus)
    x = cos.(vectus)

    # Combine the y and x coordinates into a vector of vectors
    dir = [[y[i],x[i]] for i in 1:length(y)]

return dir
end
    


"""
    expand_colony_finger_radom_cov!(img::AbstractArray, pixels_to_add::Int, dir::Vector{Int}; dir_match_rate::AbstractFloat = 0.999, still_spawn_rate::AbstractFloat = 0.99, min_neigbour::Int = 2)

Expand the colony in the image `img` by adding `pixels_to_add` pixels. The expansion is done randomly 
at the border of the colony, but it is more likely to expand in the direction specified by `dir`. 
The border is determined by convolving the image with a Laplacian kernel and finding points where the convolution 
is greater than 0.1. The function modifies the input image in-place.

# Arguments
- `img::AbstractArray`: A 2D array representing the image of the colony. The colony is represented by 1s and the background by 0s.
- `pixels_to_add::Int`: The number of pixels to add to the colony.
- `dir::Vector{Int}`: A vector representing the preferred direction of expansion.
- `dir_match_rate::AbstractFloat`: A float representing the threshold for the dot product between the preferred direction and the point. Default is 0.999.
- `still_spawn_rate::AbstractFloat`: A float representing the probability of expanding in a direction opposite to `dir`. Default is 0.99.
- `min_neigbour::Int`: The minimum number of neighboring points that must be occupied for a point to be added to the colony. Default is 2.

# Example
```julia
img = zeros(100, 100)
img[50:55, 50:55] .= 1
expand_colony_finger_radom_cov!(img, 100, [1, 0])
```
"""	
function expand_colony_finger_radom_cov!(img::AbstractArray,pixels_to_add::Int,dir; dir_match_rate = 0.999, still_spawn_rate::AbstractFloat = 0.99, min_neigbour::Int = 1 )
    pix_count = 0
    laplac_kernel = [0 1 0; 1 -4 1; 0 1 0]
    cov_img = conv( img, laplac_kernel )
    border_points = shuffle!(findall(cov_img .> 0.1))
    shuffle_counter = 0
    center = size(img).÷2
    
    
    while pix_count < pixels_to_add
        point = rand(border_points)
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

            new_border_points = findall(cov_img[point[1]-1:point[1]+1,point[2]-1:point[2]+1] .> 0.1) .+( point - CartesianIndex(2,2))
            append!(border_points, new_border_points)
        end
        shuffle_counter += 1

        if shuffle_counter >= length(border_points)*0.1
            border_points = shuffle!(findall(cov_img .> 0.1))
            shuffle_counter = 0 
        end   
        
    end
    border_points = shuffle!(findall(cov_img .> 0.1))
    shuffle_counter = 0 

end
    





"""
    expand_colony_radom_cov!(img::AbstractArray, pixels_to_add::Int)

Expand the colony in the image `img` by adding `pixels_to_add` pixels. The expansion is done randomly 
at the border of the colony. The border is determined by convolving the image with a Laplacian kernel 
and finding points where the convolution is greater than 0.1. The function modifies the input image in-place.

Compared to `expand_colony_radom!`, this function is faster for large images and many pixels to add, 
but slower for small images and fewer pixels to add. This is due to the fact that the computationally heavy convolution only needs to be 
calculated once for the whole image, whereas the distance transform in `expand_colony_radom` needs to be calculated for each iteration of the loop.

# Arguments
- `img::AbstractArray`: A 2D array representing the image of the colony. The colony is represented by 1s and the background by 0s.
- `pixels_to_add::Int`: The number of pixels to add to the colony.

# Example
```julia
img = zeros(100, 100)
img[50:55, 50:55] .= 1
expand_colony_radom_cov!(img, 100)
```
Compared to `expand_colony_radom!`, this function is faster for large images and many pixels to add, 
but slower for small images and fewer pixels to add. This is due to the fact that the computational heavy convolution only needs to be 
    calculated once for the whole image, whereas the distance transform in `expand_colony_radom` needs to be calculated for each iteration of the loop.
"""
function expand_colony_radom_cov_show!(img::AbstractArray,pixels_to_add::Int)
     # Initialize pixel count and Laplacian kernel 
    pix_count = 0 
    laplac_kernel = [0 1 0; 1 -4 1; 0 1 0]
    # Convolve the image with the Laplacian kernel
    cov_img = conv( img, laplac_kernel )

    # Find border points where the convolution is greater than 0.1 and shuffle them
    border_points = shuffle!(findall(cov_img .> 0.1))
    shuffle_counter = 0

    # Loop until the desired number of pixels have been added
    while pix_count < pixels_to_add
        point = rand(border_points)
        
        # If the point is in the background, add it to the colony
        if img[point] == 0
            img[point] = 1
            
            # Update the convolution image as well
            cov_img[point[1]-1:point[1]+1,point[2]-1:point[2]+1] += laplac_kernel
            pix_count += 1

            new_border_points = findall(cov_img[point[1]-1:point[1]+1,point[2]-1:point[2]+1] .> 0.1) .+( point - CartesianIndex(2,2))
            append!(border_points, new_border_points)
            # old border point stays in border point arrays, due to computational efficiency, but is never choosen anyways, as it does not pass the if statement
        end
        shuffle_counter += 1
        

        # recreate the border_points vector if  over 50% of the original border points have been checked
        if shuffle_counter >= length(border_points)*0.5
            border_points = shuffle!(findall(cov_img .> 0.1))
            shuffle_counter = 0 
        end            
        
    end
    return cov_img
end


"""
    growth_colonies(vec_of_sims::Vector{Vector{BitArray{3}}}, para::parameters, Points::Vector{Vector{Vector{Int}}})

Simulates the growth of colonies based on the provided parameters and initial conditions.

# Arguments
- `vec_of_sims`: A vector of 3D bit arrays representing the initial state of each colony.
- `para`: An instance of the `parameters` struct containing the simulation parameters.
- `Points`: A vector of vectors of vectors of integers representing the points in the colony.

# Details
The function first creates an artificial colony using the `build_artifical_colony!` function. It then iterates over the simulations and colonies specified in `para`, expanding each colony at each time point based on the specified simulation type ("Random", "Finger_weak", or "Finger_strong").

# Example
```julia
growth_colonies(vec_of_sims, para, Points)
```
"""	
function growth_colonies(vec_of_sims::Vector{Vector{BitArray{3}}}, para::parameters,Points::Vector{Vector{Vector{Int}}})
    # Initialize the vector to store the pixels to add and the image of the colony
    pixel_to_add_vec = Int[];
    copy_int_img = BitArray(zeros(Bool, para.im_size...,));

    # Build the initial artificial colony
    build_artifical_colony!(para.Center, copy_int_img, para.radius_colony, Points);
    pixel_to_add_vec = para.pixel_to_add(copy_int_img);

    # Iterate over the simulations and colonies
    for (x, sim) in collect(enumerate(para.simulations))
        Threads.@threads for (i, colony) in collect(enumerate(vec_of_sims[x]))
            for (j, t) in enumerate(para.time_points)
                # Use a view to work on the array directly
                int_img  = @view colony[:,:,j];
                dir_vec = generate_dir_vec(para);

                if j == 1
                    # Build the initial artificial colony
                    build_artifical_colony!(para.Center, int_img, para.radius_colony, Points);
                else
                    # Copy the state of the colony from the previous time point
                    int_img[:] = colony[:,:,j-1];

                    # Expand the colony based on the simulation type
                    if sim == "Random"
                        expand_colony_radom_cov!(int_img, pixel_to_add_vec[j-1]);
                    elseif sim == "Finger_weak"
                        expand_colony_finger_radom_cov!(int_img, pixel_to_add_vec[j-1], dir_vec, still_spawn_rate= para.spawn_rate, dir_match_rate = para.dir_match_rate_B);
                    elseif sim == "Finger_strong"
                        expand_colony_finger_radom_cov!(int_img, pixel_to_add_vec[j-1], dir_vec, still_spawn_rate= 0.0, dir_match_rate = para.dir_match_rate_C);
                    end
                end
            end
        end   
    end
end
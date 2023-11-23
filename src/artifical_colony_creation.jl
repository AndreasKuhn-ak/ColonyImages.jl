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

This function expands an artifical colony in a circular pattern within an image. The expansion starts from the center of the colony and proceeds outward. 
The function directly modifies the input image.

# Arguments
- `img::Union{Matrix{<:Real}, BitMatrix}`: The input image where the colony will be expanded.
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
    expand_colony_radom_cov!(img::Union{Matrix{<:Real}, BitMatrix}, pixels_to_add::Int)

Expand the colony in the image `img` by adding `pixels_to_add` pixels. The expansion is done randomly 
at the border of the colony. The border is determined by convolving the image with a Laplacian kernel 
and finding points where the convolution is greater than 0.1. The function modifies the input image in-place.

Compared to `expand_colony_radom!`, this function is faster for large images and many pixels to add, 
but slower for small images and fewer pixels to add. This is due to the fact that the computationally heavy convolution only needs to be 
calculated once for the whole image, whereas the distance transform in `expand_colony_radom` needs to be calculated for each iteration of the loop.

# Arguments
- `img::Union{Matrix{<:Real}, BitMatrix}`: A 2D array representing the image of the colony. The colony is represented by 1s and the background by 0s.
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
function expand_colony_radom_cov!(img::Union{Matrix{<:Real}, BitMatrix},pixels_to_add::Int)
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
        for point in border_points
            # If the point is in the background, add it to the colony
            if img[point] == 0
                img[point] = 1
                
                # Update the convolution image as well
                cov_img[point[1]-1:point[1]+1,point[2]-1:point[2]+1] += laplac_kernel
                pix_count += 1
            end
            shuffle_counter += 1
            
            # Break the loop if we have checked 10% of the border points
            if shuffle_counter >= length(border_points)*0.1
                break
            end            
        end
        
        # Update the border points and reset the shuffle counter
        border_points = shuffle!(findall(cov_img .> 0.1))
        shuffle_counter = 0 
    end
end


"""
    expand_colony_radom!(img::Union{Matrix{<:Real}, BitMatrix}, pixels_to_add::Int)

This function expands a colony in a random pattern within an image. The expansion is performed by adding pixels to the border of the colony. 
The function directly modifies the input image. 

# Arguments
- `img::Union{Matrix{<:Real}, BitMatrix}`: The input image where the colony will be expanded.
- `pixels_to_add::Int`: The number of pixels to add to the colony.

Compared to `expand_colony_radom_cov!`, this function is slower for large images and many pixels to add,
 but faster for small images and and fewer pixels to add .This is due to the fact that the computational heavy convolution only needs to be 
    calculated once for the whole image, whereas the distance transform in `expand_colony_radom` needs to be calculated for each iteration of the loop.
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
    


"""
    expand_colony_fractal_radom_cov!(img::Union{Matrix{<:Real}, BitMatrix}, pixels_to_add::Int; max_neigbours::Int, = 1 prob_neigbour::AbstractFloat = 0.3)

Expand the colony in the image `img` by adding `pixels_to_add` pixels. The expansion is done randomly 
at the border of the colony, but it is more likely to expand to points that only have `max_neigbours` or less. 
The border is determined by convolving the image with a Laplacian kernel and finding points where the convolution 
is greater than 0.1. The function modifies the input image in-place.

# Arguments
- `img::Union{Matrix{<:Real}, BitMatrix}`: A 2D array of type Matrix with elements of any subtype of Real or a BitMatrix, representing the image of the colony. The colony is represented by 1s and the background by 0s.
- `pixels_to_add::Int`: An integer representing the number of pixels to add to the colony.
- `max_neigbours`: An integer representing the maximum number of neighbours a point can have to be considered for expansion. Default is 1.
- `prob_neigbour::AbstractFloat`: A float representing the probability of expanding to a point with more than `max_neigbours`. Default is 0.3.

# Example
```julia
img = zeros(100, 100)
img[50:55, 50:55] .= 1
expand_colony_fractal_radom_cov!(img, 100)
```
"""
function expand_colony_fractal_radom_cov!(img::Union{Matrix{<:Real}, BitMatrix},pixels_to_add::Int; max_neigbours::Int = 1, prob_neigbour::AbstractFloat = 0.3)
    pix_count = 0
    laplac_kernel = [0 1 0; 1 -4 1; 0 1 0]
    # Convolve the image with the Laplacian kernel
    cov_img = conv( img, laplac_kernel )

    # Find border points where the convolution is greater than 0.1 and shuffle them
    border_points = shuffle!(findall(cov_img .> 0.1))
    shuffle_counter = 0

    # Loop until the desired number of pixels have been added
    while pix_count < pixels_to_add
        for point in border_points
            # If the point has more than max_neigbours, skip it with a probability of prob_neigbour
            if cov_img[point] > max_neigbours+ 0.1
                if rand() > prob_neigbour
                    shuffle_counter += 1
            
                    if shuffle_counter >= length(border_points)*0.1
                        break
                    end    
                    continue
                end
            end
            
            # If the point is in the background, add it to the colony
            if img[point] == 0
                img[point] = 1
                
                # Update the convolution image as well
                cov_img[point[1]-1:point[1]+1,point[2]-1:point[2]+1] += laplac_kernel
                pix_count += 1
            end
            shuffle_counter += 1
            
            # Break the loop if we have checked 10% of the border points
            if shuffle_counter >= length(border_points)*0.1
                break
            end            
        end
        
        # Update the border points and reset the shuffle counter
        border_points = shuffle!(findall(cov_img .> 0.1))
        shuffle_counter = 0 
    end
end


"""
    expand_colony_un_radom_cov!(img::Union{Matrix{<:Real}, BitMatrix}, pixels_to_add::Int, dir::Vector{Int}; still_spawn_rate::AbstractFloat = 0.5)

Expand the colony in the image `img` by adding `pixels_to_add` pixels. The expansion is done randomly 
at the border of the colony, but it is more likely to expand in the direction specified by `dir`. 
The border is determined by convolving the image with a Laplacian kernel and finding points where the convolution 
is greater than 0.1. The function modifies the input image in-place.

# Arguments
- `img::Union{Matrix{<:Real}, BitMatrix}`: A 2D array representing the image of the colony. The colony is represented by 1s and the background by 0s.
- `pixels_to_add::Int`: The number of pixels to add to the colony.
- `dir::Vector{Int}`: A vector representing the preferred direction of expansion.
- `still_spawn_rate::AbstractFloat`: A float representing the probability of expanding in a direction opposite to `dir`. Default is 0.5.

# Example
```julia
img = zeros(100, 100)
img[50:55, 50:55] .= 1
expand_colony_un_radom_cov!(img, 100, [1, 0])
```
"""
function expand_colony_un_radom_cov!(img::Union{Matrix{<:Real}, BitMatrix},pixels_to_add::Int,dir::Vector{Int}; still_spawn_rate::AbstractFloat = 0.5,  )
    pix_count = 0
    laplac_kernel = [0 1 0; 1 -4 1; 0 1 0]
    cov_img = conv( img, laplac_kernel )
    border_points = shuffle!(findall(cov_img .> 0.1))
    shuffle_counter = 0
    center = size(img).÷2
    
    
    while pix_count < pixels_to_add
        for point in border_points
            # if dot product between preferred dir and chosen point is negative only proceed in 50% of the time
            if dot(dir,point.I .- center) .< 0 && rand() > still_spawn_rate
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


"""
    generate_dir_vec(number_fingers, rand_dist)

Generates vectors of x and y coordinates that span vectors ranging from 0 to 2π. 
This function is used to generate the vectors used in `expand_colony_finger_radom_cov!`.

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
    expand_colony_finger_radom_cov!(img::Union{Matrix{<:Real}, BitMatrix}, pixels_to_add::Int, dir::Vector{Int}; dir_match_rate::AbstractFloat = 0.999, still_spawn_rate::AbstractFloat = 0.99, min_neigbour::Int = 2)

Expand the colony in the image `img` by adding `pixels_to_add` pixels. The expansion is done randomly 
at the border of the colony, but it is more likely to expand in the direction specified by `dir`. 
The border is determined by convolving the image with a Laplacian kernel and finding points where the convolution 
is greater than 0.1. The function modifies the input image in-place.

# Arguments
- `img::Union{Matrix{<:Real}, BitMatrix}`: A 2D array representing the image of the colony. The colony is represented by 1s and the background by 0s.
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
function expand_colony_finger_radom_cov!(img::Union{Matrix{<:Real}, BitMatrix},pixels_to_add::Int,dir; dir_match_rate = 0.999, still_spawn_rate::AbstractFloat = 0.99, min_neigbour::Int = 2 )
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
    




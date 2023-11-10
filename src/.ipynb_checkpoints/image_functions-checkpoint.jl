function centroid(img)
    Y,X = size(img)
    centroid = Float64[0.0,0.0]
    counter = 0
    for x in 1:X
        for y in 1:Y
            if img[y,x] != 0
                counter += img[y,x]
                centroid += img[y,x].*[y,x]
            end
        end
    end
    centroid_norm = round.(Int64,centroid./counter)
    return centroid_norm
end


"""
Calcultes ascending integer pairs which fit into a circle with radius r.
The return type is a Vector{Vector{Vector{Int}}} where as the first index is corresponds to r and the second for a individual pair
"""
function lattice_points(r::Int)
    points = Vector{Vector{Real}}(undef,0)
    for x in 0:r , y in 0:r
        if !(x == 0 && y== 0)
            if sqrt(x^2 +y^2)< r
                push!(points,[x,y,sqrt(x^2 +y^2)],[-x,y,sqrt(x^2 +y^2)],[x,-y,sqrt(x^2 +y^2)],[-x,-y,sqrt(x^2 +y^2)])
            end
        end
    end
    sort!(points, by = x-> x[3])
    points_r = [round.(Int,x[1:2]) for x in points]
    
    points3 = Vector{Vector{Vector{Int}}}(undef,0)
    j = 1 
    val = 2
    
    for i in 1:length(points)
        if points[i][3] >= val
            append!(points3,[points_r[j:i-1]])
            j = i
            val += 1
        end
    end

    return points3
end

function occupied_points(img)
    return sum(img)/(size(img)[1]*size(img)[2])
end

function approx_radi_colo(img)
    return sum(img)^(1/2)
end
    
    
function create_kernel(rad;  geometry = "circle")
    kernel = zeros( 2*rad + 1, 2*rad + 1 )
    if geometry == "circle"
        for x in 1:2*rad+1, y in 1:2*rad+1
            dist = ( y - rad + 1 )^2 + ( x - rad + 1 )^2
            kernel[ y, x ] = ( rad^2 < dist )
        end
    elseif geometry == "square"
        kernel .= 1
    end
    
    return kernel 
end

function build_circle(center, img, points;threshold = 0.8)
    circle_kernel = zeros(Int, size(img))
    r = 1
    while r < minimum(size(img))
        occupation_vec = zeros(Int,length(points[r]))
        for (i,point) in enumerate(points[r])
            point_c = point .+ center
            circle_kernel[point_c...] = 1
            if img[point_c...] == 1
                
                occupation_vec[i] = 1
            end
        end
        r += 1
        if r > 10 && mean(occupation_vec) < threshold
            break
        end
    end
    return circle_kernel 
end

function build_artifical_colony!(center, img, radius, points)
    r = 0
    while r <= radius
        if r == 0
            img[center...] = 1
        else
            for (i,point) in enumerate(points[r])
                point_c = point .+center
                img[point_c...] = 1
            end
        end
        r += 1
    end
    return img
end

function expand_colony_circular!(img,points,center,pixels_to_add)
    pix_count = 0
    for r in 1:(round(Int,maximum(size(img))))
        for (i,point) in enumerate(points[r])
            point_c = point .+center
            #if pixel is not black yet, make it black
            if img[point_c...] == 0
                img[point_c...] = 1
                #println(pix_count)
                pix_count += 1
            end
            if pix_count > pixels_to_add
                @goto escape
            end
        end
    end
    @label escape
end

function expand_colony_radom!(img,pixels_to_add)
    pix_count = 0
    border_points = shuffle!(findall((distance_transform(feature_transform(img)).== 1)))
    shuffle_counter = 0
    while pix_count < pixels_to_add
        for point in border_points
            if img[point] == 0
                img[point] = 1
                pix_count += 1
            end
            shuffle_counter += 1
            
            if shuffle_counter >= length(border_points)*0.1
                break
            end            
        end
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
Generates vectors of x and y coordinates that span vectors ranging from 0 to 2π
    Execute the following code to understand how this function works.
```
    using CairoMakie
    number_vecs = 20
    y,x = generate_dir_vec(20,0.1)
    yy = zeros(number_vecs)
    arrows(yy,yy, y, x)
```
"""
function generate_dir_vec(number_fingers, rand_dist)
    vectus = 2π/number_fingers.*[1:number_fingers...]+ rand(number_fingers) *rand_dist
    y = sin.(vectus)
    x = cos.(vectus)
    dir = [[y[i],x[i]] for i in 1:length(y)]
    return dir
end
    

function conv( img, kernel )
    csize = size(img) .+ size(kernel) .- 1
    padimg = zeros( ComplexF32, csize )
    padkernel = zeros( ComplexF32, csize )
    
    padimg[ 1:size(img,1), 1:size(img,2) ] .= img
    padkernel[ 1:size(kernel,1), 1:size(kernel,2)  ] .= kernel
    
    fft!( padimg )
    fft!( padkernel) 
    padimg .*= padkernel
    ifft!( padimg )
    
    output = real.( padimg )
    
    off = div.( size(kernel), 2 )
    h, w = size(img)
    return output[ 1+off[1]:h+off[1], 1+off[2]:w+off[2] ]
end

b_w_n(img) = (round.(Int64,Float64.(Gray.(img))).-1).*-1
b_w_i(img )= round.(Int64,Float64.(Gray.(img)))

function b_w(img)
    img = b_w_i(img)
    if sum(img) >= prod(size(img))*0.5
        img = (img .-1) .*-1 
    end
    return img
end



function fill_holes(img,size_holes)
    size_absolut = size(img)[1]*size(img)[2]*size_holes 
    
    return .!(imfill(.!(Bool.(img)),(0.1,size_absolut)))
end


"""Helper functions for angular metric, takes an angle in radiants as input and a number of circular sectors into which a unit circle is divided. The output is the number of circular sector to which the angle corresponds. E.g angle = 0.0 ,steps = 360 ;output = 1 which is first circular sector on the unit circle, EG.  angle = (2pi -0.01), steps = 360, output = 360, last circular sector ranging from 2pi -0.01), steps = 360, output = 360, last circular sector ranging from [2pi/steps *(steps-1) - 2pi)
"""
function rad2deg_discrete(ϕ; steps=360)
    stepsize = 360/steps 
    ϕ_a = rad2deg(ϕ) +180
    if ϕ_a >= 360
        ϕ_a = 0.0
    end
    
    ϕ_a_d = round(Int,(ϕ_a-0.5)/stepsize)+1
    if ϕ_a_d < 1
        ϕ_a_d =1 
    end
    if ϕ_a_d >= steps
        ϕ_a_d = 1
    end
        
    
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

function angular_metric(img, center; steps = 360)
    angle_vec = zeros(Int64, steps)
    stepsize = 2π/steps 
    for y in 1:size(img)[1], x in 1: size(img)[2]
        if img[y,x] != 0
            vec = [y,x] .- center
            ϕ = atan(vec...)+π
            step_index = minimum([round(Int,ϕ/stepsize), steps-1])
            angle_vec[step_index+1] += 1
        end
    end
    #last and first circular segment are buggy due to discretions artefacts, mean of both of them
    mean_angle = round(Int,mean([angle_vec[1], angle_vec[end]]))
    angle_vec[1], angle_vec[end] = mean_angle,mean_angle
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

function pair_cor_metric3(img, center; samples = 10000, steps = 360,counter_while = 0 )
    angle_vec = zeros(Int64, steps)
    stepsize = π/steps
    indizies = Tuple.(findall(x-> x != 0 ,img))
    angles = Float64[]
    counter_while = 0
    for i in 1:samples
        v1, v2  = sample(indizies, 2, replace = false)
        # ensure that center of the colony is not part of the sampled pair of points 
        while Tuple(center) in [v1,v2]
            v1, v2  = sample(indizies, 2, replace = false)
            counter_while += 1
        end
        v1_v = v1 .- center 
        v2_v = v2 .- center
        x = dot(v1_v,v2_v)/(norm(v1_v)*norm(v2_v))
        angle = acos(x >1.0 ? 1.0 : x < -1.0 ? -1.0 : x)
        step_index = trunc(Int64,angle/stepsize)
        # sometimes step index can be bigger than 359 for rounding reasons, if so 359 is assigned 
        step_index_final = (step_index == 360 ? 359 : step_index)
        angle_vec[step_index_final+1] += 1
       
    end
    return angle_vec
end

function res_scaling(img_int_vec; factor = 3, plots = 1)
    c = 0
    for i in img_int_vec
        for j in 1:size(i,3)
            c += 1
        end
    end
    return round(Int64,factor *1000)*plots,round(Int64,factor*200*(c÷5+((c%5 == 0) ? 0 : 1)))
end



replace_nan_1(x) = ismissing(x) || (x isa Number && isnan(x)) ? 1 : x

replace_nan(x) = ismissing(x) || (x isa Number && isnan(x)) ? 0 : x

function filter_fourier_alpha(vec; a = 5)
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
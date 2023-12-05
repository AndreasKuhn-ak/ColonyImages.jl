"""
    parameters

A struct for storing various parameters used in the analysis of Trypanosoma colony images.

# Fields
- `time_points::Vector{Float64}`: The time points in hours at which the images were taken. Default is `[0,48,96,144]`.
- `threshold_conv::Float64`: The threshold for the for convertion the convoluted image into bit array,
    all pixels . Default is `0.8`.
- `threshold_c::Float64`: The threshold for the c value. Default is `0.8`.
- `kernel_ratio::Float64`: The ratio for the kernel operation. Default is `0.4`.
- `steps_angular::Int`: The number of steps for the angular metric calculation. Default is `360`.
- `samples_pair::Int`: efines how many random pair of lattice points are sampled for the pair correlation metric. A higher number here is always better for the precision.
    The only limiting factor here is computation time (and at some point probably memory). 
    Its default values is:  `2000000`.
- `timesteps::Int64`: The number of time steps. Default is `300`.
- `im_size::Vector{Int}`: The size of the image. Default is `[600,600]`.
- `stacks::Int`: The number of stacks. Default is the length of `time_points`.
- `radius_colony::Int`: The radius of the colony. Default is `round(Int,(im_size[1]*0.05))`.
- `Center::Vector{Int}`: The center of the colony. Default is `round.(Int,im_size./2)`.
- `growth_rate::Float64`: The growth rate of the colony. Default is `0.02971700864000873`.
- `colony_size::Function`: A function to calculate the size of the colony. Default is `t-> (1+growth_rate).^t`.
- `relative_size_filles_holes::Float64`: `relative_size_filles_holes` defines how big the holes 
    inside a colony can be relative to 
    the image size, that are filled automatically with the `fill_holes` function.  Default is `0.01`.
"""
@with_kw struct parameters
    time_points::Vector{Float64}        = Float64[0,48,96,144]
    threshold_conv::Float64             = 0.8
    threshold_c::Float64                = 0.8
    kernel_ratio::Float64               = 0.4
    steps_angular::Int                  = 360
    samples_pair::Int                   = 2000000
    timesteps::Int64                    = 300
    im_size::Vector{Int}                = [600,600]
    stacks::Int                         = length(time_points)
    radius_colony::Int                  = round(Int,(im_size[1]*0.05))
    Center::Vector{Int}                 = round.(Int,im_size./2)
    growth_rate::Float64                = 0.02971700864000873
    colony_size::Function               = t-> (1+growth_rate).^t 
    relative_size_filles_holes::Float64 = 0.01
end


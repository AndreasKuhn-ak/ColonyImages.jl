"""
    parameters

A struct that holds various parameters for colony image processing, creation and analysis.

# Fields
- `time_points::Vector{Float64}`: The time points for the analysis. Default is `[0,48,96,144]`.
- `threshold_conv::Float64`: The threshold for the convolution. Default is `0.8`.
- `threshold_c::Float64`: Another threshold parameter. Default is `0.8`.
- `kernel_ratio::Float64`: The ratio for the kernel. Default is `0.4`.
- `steps_angular::Int`: The number of angular steps. Default is `360`.
- `samples_pair::Int`: The number of sample pairs. Default is `2000000`.
- `timesteps::Int64`: The number of time steps. Default is `300`.
- `im_size::Vector{Int}`: The size of the image. Default is `[600,600]`.
- `stacks::Int`: The number of stacks, which is the length of `time_points`.
- `radius_colony::Int`: The radius of the colony. Default is `round(Int,(im_size[1]*0.05))`.
- `Center::Vector{Int}`: The center of the image. Default is `round.(Int,im_size./2)`.
- `growth_rate::Float64`: The growth rate of the colony. Default is `[0.02971700864000873,0.0,.0,.0]` The 3 additional terms are only relevant for an addtional sigmodal growth component.
- `colony_size::Function`: A function to calculate the size of the colony. Default is `t-> (1+growth_rate[1]).^t + growth_rate[2] ./ (1 .+ exp.(-growth_rate[3] .* (t .- growth_rate[4])))`.
- `relative_size_filles_holes::Float64`: The relative size of filled holes. Default is `0.01`.
- `laplac_kernel::Matrix{Int}`: The Laplacian kernel. Default is `[0 1 0; 1 -4 1; 0 1 0]`.
- `colony_nr::Int`: The number of colonies. Default is `4`.
- `colonies::Vector{String}`: The names of the colonies. Default is `["Colony x artifical" for x in 1:colony_nr]`.
- `plot_factor::AbstractFloat`: The factor for plotting. Default is `2.0`.
- `Points::Vector{Vector{Vector{Int}}}`: The lattice points. Default is `lattice_points(Int(maximum(im_size)÷2))`.
- `col_size_add::Vector{Float64}`: The additional size of the colony. Default is `colony_size.(time_points).-1`.
- `col_size_add_diff::Vector{Float64}`: The difference in the additional size of the colony. Default is `col_size_add[2:end]-col_size_add[1:end-1]`.
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
    growth_rate::Union{Float64, Vector{Float64}}        = [0.02971700864000873,0.0,.0,.0]
    colony_size::Function               = t-> (1+growth_rate[1]).^t + growth_rate[2] ./ (1 .+ exp.(-growth_rate[3] .* (t .- growth_rate[4])))
    relative_size_filles_holes::Float64 = 0.01
    laplac_kernel::Matrix{Int}          = [0 1 0; 1 -4 1; 0 1 0]
    colony_nr::Int                      = 4
    colonies::Vector{String}            = (["Colony $(x) artifical" for x in 1:colony_nr])
    plot_factor::AbstractFloat          = 2.0
    col_size_add::Vector{Float64}       = colony_size.(time_points).-1
    col_size_add_diff::Vector{Float64}  = col_size_add[2:end]-col_size_add[1:end-1]
    number_finger::Int                  = 15
    finger_dist::AbstractFloat          = 0.1
    pixel_to_add::Function              = colony ->round.(Int,sum(colony[:,:,1]).*(col_size_add_diff))
    spawn_rate::AbstractFloat           = 0.2
    dir_match_rate_B::AbstractFloat     = 0.9995
    dir_match_rate_C::AbstractFloat     = 0.9995
end


"""
    analysis_parameters

A struct for holding parameters related to the analysis.

# Fields
- `plot_theme::Attributes`: An Attributes object for setting the theme of the plots. Default is a Theme object with a fontsize of 25, a size of (1000,800), a markersize of 15 for Scatter plots, and a linewidth of 4 for Lines plots.

# Example
```julia
params = analysis_parameters()
```
"""
@with_kw struct analysis_parameters
    plot_theme::Attributes              = Theme(    
                                        fontsize = 30,
                                        size = (1000,800),
                                        Scatter = (
                                        markersize = 18,
                                        ),                
                                        Lines  = ( 
                                        linewidth =4,
                                        ),
                                        Errorbars = (whiskerwidth = 20, 
                                        color = :black)
                                        )
end

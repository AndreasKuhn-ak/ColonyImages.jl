module ColonyImages

using Images, FFTW 
using Statistics

export conv, 
        b_w,
        fill_holes,
        lattice_points,
        centroid,
        approx_radi_colo,
        create_kernel,
        build_circle,
        res_scaling

include("image_functions.jl")
include("artifical_colony_creation.jl")


end

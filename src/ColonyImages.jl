module ColonyImages

using Images, FFTW 
using Statistics, StatsBase

export conv, 
        b_w,
        fill_holes,
        lattice_points,
        centroid,
        approx_radi_colo,
        create_kernel,
        build_circle,
        res_scaling,
        angular_metric,
        pair_cor_metric3,
        build_artifical_colony!,
        expand_colony_circular!


include("image_functions.jl")
include("artifical_colony_creation.jl")


end

module ColonyImages

using Images, FFTW 
using Statistics, StatsBase, Random

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
        expand_colony_circular!,
        expand_colony_radom_cov!,
        filter_fourier_beta,
        find_freq,
        expand_colony_radom_cov_show!,
        expand_colony_point!



include("image_functions.jl")
include("artifical_colony_creation.jl")


end

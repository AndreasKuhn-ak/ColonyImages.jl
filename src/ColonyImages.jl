module ColonyImages

using Images, FFTW 

export conv, 
        b_w,
        fill_holes,
        lattice_points,
        centroid,
        approx_radi_colo,
        create_kernel

include("image_functions.jl")


end

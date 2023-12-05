"""
    plot_time_series_cov_centroid(img_vec; kernel_ratio = 0.4)

Generates a time series plot of the centroids of the images in `img_vec`. The centroids are calculated in three ways: 
1. The original centroid of the first timestep of the timeseries.
2. The centroid after applying a convolution operation.
3. The current centroid of the image (for t = 1, 1. == 3.).

# Arguments
- `img_vec`: A vector of 3D image stacks.
- `kernel_ratio::Float64`: The ratio for the kernel operation. Default is `0.4`.

# Returns
- A figure with the time series plot of the centroids.

"""
function plot_time_series_cov_centroid(img_vec; kernel_ratio = 0.4)
    fig_big = Figure(size = res_scaling(img_vec, factor =3))

    c = 0 
    for (i,img_stack) in enumerate(img_vec)
        int_img = img_stack[:,:,1]
        y1,x1 = centroid(int_img)
        kernel = create_kernel(round(Int64,approx_radi_colo(int_img)*kernel_ratio), geometry = "square")
        nneigh = sum(kernel)
        
        for z in 1:size(img_stack,3)
            int_img = img_stack[:,:,z]
            out = conv( int_img, kernel ) ./ nneigh
            y_c, x_c = centroid(out .> threshold_conv)
            y,x = centroid(int_img)
            ax = CairoMakie.Axis(fig_big[c÷5+1,c%5+1], title = colonies[i])
            #hidedecorations!(ax)
            heatmap!(ax,int_img,colormap = :algae)
            heatmap!(ax,out .> threshold_conv, colormap  =(:algae, 0.2))
            scatter!(ax,y1,x1, color = :blue, markersize = 15, label = "OG centroid")
            
            scatter!(ax,y,x, color = :red, markersize = 15, label = "Current centroid")
            scatter!(ax,y_c, x_c, color = :yellow, markersize = 15, label = "Convolut centroid")
            axislegend(ax)
            c += 1
                
        end

    end
    return fig_big
end


function plot_convolution_schematic3(colony,colony_cov,para; name = "cov_fig_kernel" )
    fontsize_theme = Theme(fontsize = 35)
    update_theme!(fontsize_theme)

    t = length(para.time_points)
    data = colony[:,:,t]
    data_cov = colony_cov[:,:,t]
    laplac_kernel = [0 1 0; 1 -4 1; 0 1 0]
    kernel  = zeros(Int,para.im_size...)
    kernel[(1:size(laplac_kernel,1)).+3 , (1:size(laplac_kernel,2)).+4] .= laplac_kernel

    cov_fig_k = Figure(size = (1800,600))
    col_ax_k = CairoMakie.Axis(cov_fig_k[1,1],aspect = DataAspect(),title = "Colony")
    col_heat_k = heatmap!(col_ax_k,data, colormap = (:blues,0.3))
    hidedecorations!(col_ax_k)
    cov_ax_k2 = CairoMakie.Axis(cov_fig_k[1,2],aspect = DataAspect(),title = "Kernel")
    cov_heat_k2 = heatmap!(cov_ax_k2,kernel, colormap = (:blues,0.3))
    hidedecorations!(cov_ax_k2)


    cov_ax = CairoMakie.Axis(cov_fig_k[1,3],aspect = DataAspect(),title = "Convoluted Colony")
    cov_heat = heatmap!(cov_ax,data_cov, colormap = (:blues,0.3))
    hidedecorations!(cov_ax)


    for i in 1:para.im_size[1]
        for j in 1:para.im_size[2]
            text!(cov_ax_k2,"$(round(Int,kernel[i,j]))", position = (i, j),
            color = "black", align = (:center, :center))
            text!(col_ax_k,"$(round(Int,data[i,j]))", position = (i, j),
            color = "black", align = (:center, :center))
            text!(cov_ax,"$(round(Int,data_cov[i,j]))", position = (i, j),
            color = "black", align = (:center, :center))

        end
    end

    save("plots/$(name).png", cov_fig_k)
    return cov_fig_k
end


function plot_convolution_schematic2(colony,colony_cov,para; name = "cov_fig1" )
    t = length(para.time_points)

    if typeof(colony) == BitArray{3}
        data = colony[:,:,t]
    else
        data = colony
    end
    if typeof(colony_cov) == Array{Float64, 3}
        data_cov = colony_cov[:,:,t]
    else
        data_cov = colony_cov
    end
    
    cov_fig = Figure(size = (1200,600))
    col_ax = CairoMakie.Axis(cov_fig[1,1],aspect = DataAspect(),title = "Colony")
    col_heat = heatmap!(col_ax,data, colormap = (:blues,0.3))
    hidedecorations!(col_ax)
    cov_ax = CairoMakie.Axis(cov_fig[1,2],aspect = DataAspect(),title = "Convoluted Colony")
    cov_heat = heatmap!(cov_ax,data_cov, colormap = (:blues,0.3))
    hidedecorations!(cov_ax)
    
    
    for i in 1:para.im_size[1]
        for j in 1:para.im_size[2]
            text!(col_ax,"$(round(Int,data[i,j]))", position = (i, j),
            color = "black", align = (:center, :center))
            text!(cov_ax,"$(round(Int,data_cov[i,j]))", position = (i, j),
            color = "black", align = (:center, :center))
        end
    end
    
    save("plots/$(name).png", cov_fig)
    return cov_fig
end
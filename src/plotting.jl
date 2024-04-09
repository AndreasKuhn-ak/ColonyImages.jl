"""
    res_scaling(img_int_vec; factor = 3, plots = 1)

Scales the resolution of a plot based on a given factor.
The function counts the number of images in the given image_vec and and scales the resolution 
of plot containg all these images accordingly.

# Arguments
- `img_int_vec`: A vector of images.
- `factor`: The scaling factor. Defaults to 3.
- `plots`: The number of plots per images. Defaults to 1.

# Returns
- A tuple containing the scaled width and height of the image.
"""
function res_scaling(img_int_vec::AbstractArray; factor::Real = 3.0, plots::Int = 1)
    # Initialize the counter
    c = 0

    # Iterate over the images
    for i in img_int_vec
        # Iterate over the slices in the image
        for j in 1:size(i,3)
            # Increment the counter
            c += 1
        end
    end

    # Calculate the scaled width and height of the image
    width = round(Int64, factor * 1000) * plots
    height = round(Int64, factor * 200 * (c ÷ 5 + ((c % 5 == 0) ? 0 : 1)))

    return width, height
end



"""
    plot_time_series_cov_centroid(img_vec::AbstractArray, para::parameters)

Generates a time series plot of the centroids of the images in `img_vec`. The centroids are calculated in three ways: 
1. The original centroid of the first timestep of the timeseries.
2. The centroid after applying a convolution operation.
3. The current centroid of the image.

# Arguments
- `img_vec::AbstractArray`: A vector of 3D image stacks.
- `para::parameters`: A parameters object containing various parameters for the analysis. 
  - `kernel_ratio::Float64`: The ratio for the kernel operation. Default is `0.4`.
  - `plot_factor`: The scaling factor for the plot size.
  - `threshold_conv`: The threshold for the convolution operation.
  - `colonies`: The names of the colonies.

# Returns
- `fig_big::Figure`: A Figure object with the time series plot of the centroids.

# Example
```julia
img_vec = [rand(100, 100, 10) for _ in 1:5]
para = parameters(kernel_ratio = 0.4, plot_factor = 1, threshold_conv = 0.5, colonies = ["Colony i" for i in 1:5])
fig = plot_time_series_cov_centroid(img_vec, para)
```
"""
function plot_time_series_cov_centroid(img_vec::AbstractArray, para::parameters)
    fig_big = Figure(size = res_scaling(img_vec, factor = para.plot_factor))

    c = 0 
    for (i,img_stack) in enumerate(img_vec)
        int_img = img_stack[:,:,1]
        y1,x1 = centroid(int_img)
        kernel = create_kernel(round(Int64,approx_radi_colo(int_img)*para.kernel_ratio), geometry = "square")
        nneigh = sum(kernel)
        
        for z in 1:size(img_stack,3)
            int_img = img_stack[:,:,z]
            out = conv( int_img, kernel ) ./ nneigh
            y_c, x_c = centroid(out .> para.threshold_conv)
            y,x = centroid(int_img)
            ax = CairoMakie.Axis(fig_big[c÷5+1,c%5+1], title = para.colonies[i])
            #hidedecorations!(ax)
            heatmap!(ax,int_img,colormap = :algae)
            heatmap!(ax,out .> para.threshold_conv, colormap  =(:algae, 0.2))
            scatter!(ax,y1,x1, color = :blue, markersize = 15, label = "OG centroid")
            
            scatter!(ax,y,x, color = :red, markersize = 15, label = "Current centroid")
            scatter!(ax,y_c, x_c, color = :yellow, markersize = 15, label = "Convolut centroid")
            axislegend(ax)
            c += 1
                
        end

    end
    return fig_big
end


"""
    plot_time_series_metrics(img_vec, para)

This function creates a series of plots for each image in a stack, showing the original image, the angular metric, and the pair correlation metric. 

# Arguments
- `img_vec::Array`: A 4D array where the first two dimensions are the image dimensions, the third dimension is the time point, and the fourth dimension is the image stack.
- `para::parameters`: A parameters object containing various parameters for the analysis.

# Returns
- `fig_big::Figure`: A Figure object containing the plots.

# Details
The function first creates a Figure object with a size determined by the `res_scaling` function. It then loops over each image stack in `img_vec`. For each stack, it calculates the centroid of the first image and creates a kernel based on the approximate radius of the colony in the image.

The function then loops over each time point in the image stack. For each time point, it calculates the centroid of the convoluted image and calculates the angular metric and pair correlation metric for the original image and the image minus the first image in the stack.

The function then fits a circle to the image and calculates the angular metric and pair correlation metric for the image minus the fitted circle.

The function then creates three Axis objects for the original image, the angular metric, and the pair correlation metric, and adds plots to each axis. If the time point is not the first, it also plots the angular metric and pair correlation metric for the image minus the fitted circle.

Finally, the function increments a counter and returns the Figure object after looping over all image stacks and time points.
"""
function plot_time_series_metrics(img_vec::AbstractArray, para::parameters)
    a = 1
    fig_big = Figure(size = res_scaling(img_vec, factor =para.plot_factor, plots = 3))
    c = 0 

    for (i,img_stack) in enumerate(img_vec)
        int_img = img_stack[:,:,1]
        y1,x1 = centroid(int_img)
        kernel = create_kernel(round(Int64,approx_radi_colo(int_img)*para.kernel_ratio), geometry = "square")
        nneigh = sum(kernel)
        
        for z in 1:size(img_stack,3)
            int_img = img_stack[:,:,z]
            out = conv( int_img, kernel ) ./ nneigh
            y_c, x_c = centroid(out .> para.threshold_conv)
            ang_mec_og = angular_metric(int_img .- img_stack[:,:,1],[y1,x1], steps = para.steps_angular)
            pair_mec_og = pair_cor_metric3(z == 1 ? int_img : int_img .- img_stack[:,:,1],[y1,x1], steps = para.steps_angular,samples = para.samples_pair)
            
            fitted_circle = build_circle([y_c, x_c], int_img, para.Points, threshold = para.threshold_c)
            ang_mec_conv = angular_metric(z == 1 ? int_img : int_img .- fitted_circle , [y1,x1], steps = para.steps_angular )
            
            pair_mec_conv = pair_cor_metric3(z == 1 ? int_img : int_img.- fitted_circle, [y_c, x_c], steps = para.steps_angular,samples = para.samples_pair )
            
            ax = CairoMakie.Axis(fig_big[c÷5+1,(c%5+1)*3-2], title = para.colonies[i])
            heatmap!(ax,int_img,colormap = :algae)
            #heatmap!(ax,out .> threshold_conv, colormap  =(:algae, 0.2))
            scatter!(ax,y1,x1, color = :blue, markersize = 10, label = "first centroid")
            scatter!(ax,y_c, x_c, color = :yellow, markersize = 10, label = "Convolut centroid")
            axislegend(ax)
            
            ax2 = CairoMakie.Axis(fig_big[c÷5+1,(c%5+1)*3-1], title = para.colonies[i])
            lines!(ax2,ang_mec_og, label = "OG angular metric" )
            
            ax3 = CairoMakie.Axis(fig_big[c÷5+1,(c%5+1)*3], title = para.colonies[i])
            lines!(ax3,pair_mec_og, label = "OG pair metric" )
            # convoluted angular metrix does not make sense for colony without expansion 
            if z != 1
                lines!(ax2,ang_mec_conv, label = "conv angular metric" )
                lines!(ax3,pair_mec_conv, label = "conv pair metric" )
            end
            axislegend(ax2)
            axislegend(ax3)
            #push!(df,[data_set,colonies[i],time_points[z],ang_mec_og,ang_mec_conv,pair_mec_og,pair_mec_conv,sum(img_stack[:,:,1])])
            c += 1
        end
    end

    return fig_big
end



"""
    plot_and_save_time_series_metrics!(img_vec::AbstractArray, para::parameters, df::DataFrame; name_data ="random_growth")

This function creates a series of plots for each image in a stack, showing the original image, the angular metric, and the pair correlation metric. It also updates a DataFrame with the metrics.

# Arguments
- `img_vec::AbstractArray`: A 4D array where the first two dimensions are the image dimensions, the third dimension is the time point, and the fourth dimension is the image stack.
- `para::parameters`: A parameters object containing various parameters for the analysis. 
  - `plot_factor`: The scaling factor for the plot size.
  - `kernel_ratio::Float64`: The ratio for the kernel operation.
  - `threshold_conv`: The threshold for the convolution operation.
  - `steps_angular`: The number of steps for the angular metric calculation.
  - `samples_pair`: The number of samples for the pair correlation metric calculation.
  - `Points`: The number of points for the circle fitting.
  - `threshold_c`: The threshold for the circle fitting.
  - `colonies`: The names of the colonies.
- `df::DataFrame`: A DataFrame to update with the metrics.
- `name_data::String`: A string to prepend to the date for the data set name in the DataFrame. Default is "random_growth".

# Returns
- `fig_big::Figure`: A Figure object containing the plots.

# Details
The function first creates a Figure object with a size determined by the `res_scaling` function. It then loops over each image stack in `img_vec`. For each stack, it calculates the centroid of the first image and creates a kernel based on the approximate radius of the colony in the image.

The function then loops over each time point in the image stack. For each time point, it calculates the centroid of the convoluted image and calculates the angular metric and pair correlation metric for the original image and the image minus the first image in the stack.

The function then fits a circle to the image and calculates the angular metric and pair correlation metric for the image minus the fitted circle.

The function then creates three Axis objects for the original image, the angular metric, and the pair correlation metric, and adds plots to each axis. If the time point is not the first, it also plots the angular metric and pair correlation metric for the image minus the fitted circle.

Finally, the function increments a counter, updates the DataFrame with the metrics, and returns the Figure object after looping over all image stacks and time points.
"""
function plot_and_save_time_series_metrics!(img_vec::AbstractArray, para::parameters,df::DataFrame; name_data ="random_growth" )
    a = 1
    fig_big = Figure(size = res_scaling(img_vec, factor =para.plot_factor, plots = 3))
    c = 0 
    data_set =  name_data * " $(Dates.format(now(), "yyyy_mm_dd"))"
    for (i,img_stack) in enumerate(img_vec)
        int_img = img_stack[:,:,1]
        y1,x1 = centroid(int_img)
        kernel = create_kernel(round(Int64,approx_radi_colo(int_img)*para.kernel_ratio), geometry = "square")
        nneigh = sum(kernel)
        
        for z in 1:size(img_stack,3)
            int_img = img_stack[:,:,z]
            out = conv( int_img, kernel ) ./ nneigh
            y_c, x_c = centroid(out .> para.threshold_conv)
            ang_mec_og = angular_metric(int_img .- img_stack[:,:,1],[y1,x1], steps = para.steps_angular)
            pair_mec_og = pair_cor_metric3(z == 1 ? int_img : int_img .- img_stack[:,:,1],[y1,x1], steps = para.steps_angular,samples = para.samples_pair)
            
            fitted_circle = build_circle([y_c, x_c], int_img, para.Points, threshold = para.threshold_c)
            ang_mec_conv = angular_metric(z == 1 ? int_img : int_img .- fitted_circle , [y1,x1], steps = para.steps_angular )
            
            pair_mec_conv = pair_cor_metric3(z == 1 ? int_img : int_img.- fitted_circle, [y_c, x_c], steps = para.steps_angular,samples = para.samples_pair )
            
            ax = CairoMakie.Axis(fig_big[c÷5+1,(c%5+1)*3-2], title = para.colonies[i])
            heatmap!(ax,int_img,colormap = :algae)
            #heatmap!(ax,out .> threshold_conv, colormap  =(:algae, 0.2))
            scatter!(ax,y1,x1, color = :blue, markersize = 10, label = "first centroid")
            scatter!(ax,y_c, x_c, color = :yellow, markersize = 10, label = "Convolut centroid")
            axislegend(ax)
            
            ax2 = CairoMakie.Axis(fig_big[c÷5+1,(c%5+1)*3-1], title = para.colonies[i])
            lines!(ax2,ang_mec_og, label = "OG angular metric" )
            
            ax3 = CairoMakie.Axis(fig_big[c÷5+1,(c%5+1)*3], title = para.colonies[i])
            lines!(ax3,pair_mec_og, label = "OG pair metric" )
            # convoluted angular metrix does not make sense for colony without expansion 
            if z != 1
                lines!(ax2,ang_mec_conv, label = "conv angular metric" )
                lines!(ax3,pair_mec_conv, label = "conv pair metric" )
            end
            axislegend(ax2)
            axislegend(ax3)
            
            push!(df,[data_set,para.colonies[i],para.time_points[z],ang_mec_og,ang_mec_conv,pair_mec_og,pair_mec_conv,sum(img_stack[:,:,1])])
            c += 1
        end
    end

    return fig_big, data_set
end





"""
    save_time_series_metrics!(img_vec::AbstractArray, para::parameters, df::DataFrame; name_data="random_growth")

This function calculates and saves time series metrics for a given set of images.

# Arguments
- `img_vec::AbstractArray`: An array of image stacks.
- `para::parameters`: An object containing various parameters for the analysis.
- `df::DataFrame`: The DataFrame to which the calculated metrics will be appended.
- `name_data::String`: (optional) A string representing the name of the data set. Default is "random_growth".

# Returns
- `data_set::String`: The name of the data set.
"""
function save_time_series_metrics!(img_vec::AbstractArray, para::parameters,df::DataFrame; name_data ="random_growth" )
    data_set =  name_data * " $(Dates.format(now(), "yyyy_mm_dd"))"
    for (i,img_stack) in enumerate(img_vec)
        int_img_og = img_stack[:,:,1]
        y1,x1 = centroid(int_img_og)
        kernel = create_kernel(round(Int64,approx_radi_colo(int_img_og)*para.kernel_ratio), geometry = "square")
        nneigh = sum(kernel)
        
        data_to_push = Vector{Any}(undef,size(img_stack,3) )
        Threads.@threads for z in 1:size(img_stack,3)
            int_img = img_stack[:,:,z]
            out = conv( int_img, kernel ) ./ nneigh
            y_c, x_c = centroid(out .> para.threshold_conv)
            ang_mec_og = angular_metric(int_img .- int_img_og,[y1,x1], steps = para.steps_angular)
            pair_mec_og = pair_cor_metric3(z == 1 ? int_img : int_img .- int_img_og,[y1,x1], steps = para.steps_angular,samples = para.samples_pair)
            
            fitted_circle = build_circle([y_c, x_c], int_img, para.Points, threshold = para.threshold_c)
            ang_mec_conv = angular_metric(z == 1 ? int_img : int_img .- fitted_circle , [y1,x1], steps = para.steps_angular )
            
            pair_mec_conv = pair_cor_metric3(z == 1 ? int_img : int_img.- fitted_circle, [y_c, x_c], steps = para.steps_angular,samples = para.samples_pair )
            
            data_to_push[z] = [data_set,para.colonies[i],para.time_points[z],ang_mec_og,ang_mec_conv,pair_mec_og,pair_mec_conv,sum(int_img_og)]
            
        end
        #push each time point of the colony to the dataframe
        for data in data_to_push
            push!(df,data)
        end
    end
    return data_set
end


"""
    plot_convolution_schematic3(colony, colony_cov, para; name = "cov_fig_kernel")

This function creates a schematic plot of a convolution operation on a colony image. 

# Arguments
- `colony::Array`: A 3D array representing the colony image stack.
- `colony_cov::Array`: A 3D array representing the convoluted colony image stack.
- `para::parameters`: A parameters object containing various parameters for the analysis.
- `name::String`: An optional name for the output plot. Default is "cov_fig_kernel".

# Returns
- `cov_fig_k::Figure`: A Figure object containing the plots.

# Details
The function first sets up a theme for the plot with a specific font size. It then extracts the data for the last time point from the colony and convoluted images. It creates a kernel and places it in a larger matrix of zeros.

The function then creates a Figure object and three Axis objects for the original colony, the kernel, and the convoluted colony. It adds heatmaps to each of these axes using the respective data.

The function then iterates over each point in the images and adds a text object to each point in the heatmaps, displaying the value at that point.

Finally, the function saves the figure as a PNG image and returns the Figure object.
"""
function plot_convolution_schematic3(colony::AbstractArray,colony_cov::AbstractArray,para::parameters; name = "cov_fig_kernel" )
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

"""
    plot_convolution_schematic2(colony::AbstractArray, colony_cov::AbstractArray, para::parameters; name = "cov_fig1")

This function creates a plot of a colony image and its convolution. 

# Arguments
- `colony::AbstractArray`: A 3D array representing the colony image stack or a 2D array representing a single colony image.
- `colony_cov::AbstractArray`: A 3D array representing the convoluted colony image stack or a 2D array representing a single convoluted colony image.
- `para::parameters`: A parameters object containing various parameters for the analysis.
- `name::String`: An optional name for the output plot. Default is "cov_fig1".

# Returns
- `cov_fig::Figure`: A Figure object containing the plots.

# Details
The function first determines the type of the `colony` and `colony_cov` inputs and extracts the data for the last time point if they are 3D arrays. 

It then creates a Figure object and two Axis objects for the original colony and the convoluted colony. It adds heatmaps to each of these axes using the respective data.

The function then iterates over each point in the images and adds a text object to each point in the heatmaps, displaying the value at that point.

Finally, the function saves the figure as a PNG image and returns the Figure object.
"""
function plot_convolution_schematic2(colony::AbstractArray,colony_cov::AbstractArray,para::parameters; name = "cov_fig1" )
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



"""
    plot_metric_schematic(colony::AbstractArray, para::parameters, colormap ; name = "metric_fig1")

This function creates a schematic plot of a colony image and how it is cut into angle section to be analysed later. 

# Arguments
- `colony::AbstractArray`: A  3D array representing the colony image stack.
- `para::parameters`: A parameters object containing various parameters for the analysis.
- `colormap`: A colormap to use for the heatmap.
- `name::String`: An optional name for the output plot. Default is "metric_fig1".

# Returns
- `fig::Figure`: A Figure object containing the plots.

# Details
The function first creates a Figure object and extracts the first and last images from the colony stack. It creates a kernel based on the approximate radius of the colony and the kernel ratio parameter.

The function then cuts a portion of the last image and rotates it 90 degrees. It creates two Axis objects on the Figure and hides their decorations.

The function then cuts a portion of the first image and rotates it 90 degrees. It adds a heatmap to the first Axis, showing the difference between the last and first images.

The function then creates a pie chart on the second Axis, with each slice having an equal value, representing the angle sections of the colony.

Finally, the function returns the Figure object.
"""
function plot_metric_schematic(colony::AbstractArray, para::parameters, colormap ; name = "metric_fig1")
    fig = Figure(size =(1500,round(Int,1500)), fontsize = 28)
    int_img = colony[:,:,1]
    kernel = create_kernel(round(Int64,approx_radi_colo(int_img)*para.kernel_ratio), geometry = "square")
    nneigh = sum(kernel)
    
    int_img = colony[:,:,end]
    
    cut_fac =450
    
    int_img = int_img[cut_fac:end-cut_fac,cut_fac:end-cut_fac ]
    int_img = rotr90(int_img)
    
    ax = CairoMakie.Axis(fig[1,1])
    ax2 = CairoMakie.Axis(fig[1,1])
    hidedecorations!(ax2)
    hidedecorations!(ax)
    
    int_img_1 = colony[:,:,1]
    int_img_1 = int_img_1[cut_fac:end-cut_fac,cut_fac:end-cut_fac ]
    
    int_img_1 = rotr90(int_img_1)
    
    heatmap!(ax,int_img.-int_img_1,colormap =colormap, alpha = 0.95)
    
    pie_data = [1 for i in 1:10:360]
    pie!(ax2,pie_data,radius = 10, color = (:red, 0.01))
    #axislegend(ax)
    return fig
end

"""
    plot_timeseries_heatmap(colony::AbstractArray, para::parameters; name = "Eden Growth Model", colormap  = reverse(co.Blues))

This function creates a heatmap plot of a colony image stack over time.

# Arguments
- `colony::AbstractArray`: A 3D array representing the colony image stack.
- `para::parameters`: A parameters object containing various parameters for the analysis.
- `name::String`: An optional name for the output plot. Default is "Eden Growth Model".
- `colormap`: An optional colormap to use for the heatmap. Default is the reverse of Blues.

# Returns
- `fig_eden::Figure`: A Figure object containing the heatmap plot.

# Details
The function first creates a Figure object and an Axis object with the title set to the provided name. It then initializes an intensity image with zeros.

The function then iterates over the z-dimension of the colony stack, adding each image to the intensity image.

It then creates a heatmap on the Axis, showing the intensity image with the maximum intensity subtracted and the sign reversed.

The function then hides the decorations of the Axis and adds a colorbar to the Figure, with the ticks set to the time points and the label set to "time [h]".

Finally, the function saves the Figure as a PDF in the plots folder and returns the Figure object.
"""
function plot_timeseries_heatmap(colony::AbstractArray, para::parameters; name = "Eden Growth Model", colormap  = reverse(co.Blues))

    fig_eden = Figure(size= Tuple(para.im_size))
    ax = CairoMakie.Axis(fig_eden[1,1],aspect = DataAspect(),title = name)
    int_img = zeros(Int,para.im_size...)
    for z in 1:size(colony,3)
        int_img += colony[:,:,z]
    end 
    heaty = heatmap!(ax,(int_img .-(maximum(int_img))) .*-1 , colormap = colormap)
    hidedecorations!(ax)
    Colorbar(fig_eden[1,2], heaty ,ticks = (0:(length(para.time_points)-1),string.(para.time_points)), label = "time [h]",tellheight=true)
    rowsize!(fig_eden.layout, 1, Aspect(1, 1))

    save("plots/$(name).pdf",fig_eden)
    return fig_eden
end


function mk_label(data_set::String)
    return split(data_set," ")[1]
end

"""
@h methodname

Outputs documentations in jupyternotenbooks in VScode as markdown without bugs.

Example of how to use the `@h` macro:
```julia
@h res_scaling(img_int_vec; factor = 3, plots = 1)
```

Outputs documentations in jupyternotenbooks in VScode as markdown without bugs.
"""
macro h(x)
    quote
        display("text/markdown", @doc $x)
    end    
end

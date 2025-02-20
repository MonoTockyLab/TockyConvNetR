#' Integrated Grad-CAM Feature Cell Identification
#'
#' Identifies significant feature cells using Grad-CAM data either with visualization (gating mode)
#' or as a lightweight operation using pre-converted images (base mode).
#'
#' @param mode Operation mode ("gating" for visualization, "base" for lightweight analysis)
#' @param x TockyPrepData object (required for "gating" mode)
#' @param results convert_to_image output (required for "base" mode)
#' @param feature_matrix Feature intensity matrix from Grad-CAM analysis
#' @param percentile Significance threshold percentile (0-1). If NULL, grad-CAM values will be returned, instead of feature cell designation.
#' @param n_resolution Binning resolution (for "gating" mode only)
#' @param ncol The number of the columns for the output plot
#' @param nrow The number of the rows for the output plot
#' @param filename Optional PDF output path (for "gating" mode only)
#' @param transpose Logical Whether to tranpose feature_matrix input. Note that TockyConvNetPy
#'  output Grad-CAM matrix for feature_matrix typically needs to be transposed. The default is TRUE.
#'
#' @return A binary numeric vector (1/0) for feature and other cells. If `percentile = NULL`, grad-CAM values are returned as a numeric vector.
#'
#' @examples
#' \dontrun{
#' # Gating mode with visualization
#' out <- inverseGradCAM(mode = "gating", x = prep_data,
#'                          feature_matrix = cam_matrix, filename = "output.pdf")
#'
#' # Base mode with pre-converted results
#' img_data <- convert_to_image(merged_data)
#' out <- inverseGradCAM(mode = "base", results = img_data,
#'                          feature_matrix = cam_matrix)
#' }
#' @export
#' @importFrom stats quantile
#' @importFrom grDevices dev.new pdf dev.off
#' @importFrom graphics par
#' @importFrom fields image.plot

inverseGradCAM <- function(x = NULL, results = NULL, feature_matrix,
    mode = c("gating", "base"),percentile = 0.9,n_resolution = 100,
    transpose = TRUE, filename = NULL, ncol = 2, nrow = 2) {

    mode <- match.arg(mode)
    
    if (mode == "gating") {

        if (!inherits(x, 'TockyPrepData')) {
            stop("x must be a TockyPrepData object in 'gating' mode")
        }
        
        expression_data <- x@Data
        allcells <- expression_data$ID <- seq_len(nrow(expression_data))
        sampledefinition <- x@sampledef$sampledef
        merged_data <- merge(expression_data, sampledefinition, by = "file")
        sample_size <- nrow(sampledefinition)
        results <- convert_to_image(merged_data, n_resolution)

    } else {
        # Base mode parameter checks
        if (is.null(results)) {
            stop("results must be provided in 'base' mode")
        }
        sample_size <- length(results$matrices_list)
    }
    
    feature_matrix <- as.matrix(feature_matrix)
    if(transpose){
        feature_matrix <- t(feature_matrix)
    }
    
    expression_data <- x@Data
    expression_data <- expression_data[!is.na(expression_data$Angle),]
    expression_data$ID <- 1:nrow(expression_data)
    results <- convert_to_image(expression_data, n_resolution)
    bin_edges_Intensity <- results$bin_edges_Intensity
    bin_edges_Angle <- results$bin_edges_Angle
    
    row_idx <- findInterval(expression_data[['Angle']],     bin_edges_Angle)
    col_idx <- findInterval(expression_data[['Intensity']], bin_edges_Intensity)
    
    expression_data[['row_index']] <- row_idx
    expression_data[['col_index']] <- col_idx
    
    valid_mask <- row_idx > 0 & row_idx <= nrow(feature_matrix) &
    col_idx > 0 & col_idx <= ncol(feature_matrix)
    
    expression_data[['gradCAM_value']][valid_mask] <- feature_matrix[cbind(row_idx[valid_mask], col_idx[valid_mask])]

    if(is.null(percentile)){
        return(expression_data[['gradCAM_value']])
    }else{
        threshold <- quantile(expression_data[['gradCAM_value']], probs = percentile, na.rm = TRUE)
        result_vector <- ifelse(expression_data[['gradCAM_value']] > threshold, 1, 0)
        return(as.numeric(result_vector))
    }
    
    
}


#' Inverse 2D Pixel Mapping
#'
#' This internal function maps matrix indices back to their corresponding
#' angle and intensity intervals based on the provided bin edges. It is
#' primarily used to interpret the indices of a 2D histogram or matrix
#' representation back to the original data space.
#'
#' @param row_idx An integer representing the row index in the matrix,
#' corresponding to a specific intensity bin.
#' @param col_idx An integer representing the column index in the matrix,
#' corresponding to a specific angle bin.
#' @param bin_edges_Angle A numeric vector of bin edges for the angle dimension.
#' The length of this vector should be one more than the number of angle bins.
#' @param bin_edges_Intensity A numeric vector of bin edges for the intensity
#' dimension. The length of this vector should be one more than the number of
#' intensity bins.
#'
#' @return A list containing two elements: angle_interval and
#' intensity_interval. Each element is a numeric vector of length 2,
#' representing the lower and upper bounds of the interval for the angle
#' and intensity, respectively.
#'
#' @examples
#' \dontrun{
#' pixel_2d_inverse(10, 20, bin_edges_Angle, bin_edges_Intensity)
#'}
#'
#' @keywords internal

pixel_2d_inverse <- function(row_idx, col_idx, bin_edges_Angle, bin_edges_Intensity) {
    if (row_idx < 1 || row_idx > length(bin_edges_Intensity) - 1 ||
        col_idx < 1 || col_idx > length(bin_edges_Angle) - 1) {
        stop("Row or Column index out of bounds!")
    }
    
    angle_interval <- c(bin_edges_Angle[col_idx], bin_edges_Angle[col_idx + 1])
    intensity_interval <- c(bin_edges_Intensity[row_idx], bin_edges_Intensity[row_idx + 1])
    
    return(list(angle_interval = angle_interval, intensity_interval = intensity_interval))
}



#' Convert Data to Image Matrices
#'
#' This function processes a dataset containing multiple subsets, each identified by a unique 'file' identifier, and converts each subset into a matrix based on binned 'Angle' and 'Intensity' values. It ensures consistent binning across all subsets by using global minimum and maximum values, facilitating direct comparison between the resulting matrices.
#'
#' @param data A data frame or list of data frames containing at least three columns: 'file' (unique identifier for each subset), 'Angle', and 'Intensity'. Each row represents an observation with a specific angle and intensity value.
#' @param n_resolution The number of bins

#'
#' @return A list containing:
#'   - matrices_list: A list where each element is a data frame corresponding to a subset identified by 'file', containing the original data.
#'   - bin_edges_Angle: Numeric vector of bin edges used for the 'Angle' dimension.
#'   - bin_edges_Intensity: Numeric vector of bin edges used for the 'Intensity' dimension.
#'   - counts_list: A list of matrices where each matrix represents the binned counts of 'Angle' and 'Intensity' for a subset.
#'
#' @examples
#' \dontrun{
#' # Assume 'data' is your dataset containing 'file', 'Angle', and 'Intensity' columns
#' result <- convert_to_image(data)
#'
#' # Access the binned counts matrix for the first file
#' first_counts_matrix <- result$counts_list[[1]]
#'
#' # Plot the matrix as an image
#' image(result$bin_edges_Angle, result$bin_edges_Intensity, first_counts_matrix)
#'}
#' @export

convert_to_image <- function(data, n_resolution = 100){
    if(is.null(data$ID)){
        stop('Include ID in your data \n')
    }
    long_matrix <- data
    files <- unique(long_matrix$file)
    matrices_list <- as.list(1:length(files))
    subset_data <- matrices_list

    for(i in 1:length(matrices_list)){
      subset_data[[i]] <- long_matrix[long_matrix$file == files[i],]
     # subset_data[[i]] <- mutate(subset_data[[i]], ID = rownames(long_matrix)[long_matrix$file == files[i]])

    }

    matrices_list <- subset_data
    names(matrices_list) <- files
    global_min_angle <- min(sapply(matrices_list, function(df) min(df$Angle, na.rm = TRUE)), na.rm = TRUE)
    global_max_angle <- max(sapply(matrices_list, function(df) max(df$Angle, na.rm = TRUE)), na.rm = TRUE)

    global_min_intensity <- min(sapply(matrices_list, function(df) min(df$Intensity, na.rm = TRUE)), na.rm = TRUE)
    global_max_intensity <- max(sapply(matrices_list, function(df) max(df$Intensity, na.rm = TRUE)), na.rm = TRUE)

    pixel_2d <- function(data, n_resolution = 100) {
        bin_edges_Angle <- seq(global_min_angle, global_max_angle, length.out = n_resolution)
        bin_edges_Intensity <- seq(global_min_intensity, global_max_intensity, length.out = n_resolution)
        bin_indices <- .bincode2D(data$Intensity, data$Angle, bin1 = bin_edges_Intensity, bin2 = bin_edges_Angle)
        counts_matrix <- matrix(0, nrow = n_resolution, ncol = n_resolution)
        for (i in 1:length(bin_indices$bin1)) {
          col_idx <- bin_indices$bin2[i]   # columns (i.e. x-axis in image.plot) =  Angle
          row_idx <- bin_indices$bin1[i]  # rows (i.e. y-axis in image.plot) =  Intensity
          counts_matrix[row_idx, col_idx] <- counts_matrix[row_idx, col_idx] + 1
        }
        return(list(counts = counts_matrix, bin_edges_Angle = bin_edges_Angle, bin_edges_Intensity = bin_edges_Intensity))
      }


    .bincode2D <- function(x, y, bin1, bin2) {
        bin1_indices <- findInterval(x, bin1)
        bin2_indices <- findInterval(y, bin2)
        return(list(bin1 = bin1_indices, bin2 = bin2_indices))
      }

    bin_matrices_list <- lapply(matrices_list, pixel_2d)

    counts_list <- lapply(bin_matrices_list, function(x) x$counts)
    data_matrix <- do.call(cbind, lapply(counts_list, as.vector))
    
    out <- list(matrices_list = matrices_list,
    bin_edges_Angle = bin_matrices_list[[1]]$bin_edges_Angle,
    bin_edges_Intensity = bin_matrices_list[[1]]$bin_edges_Intensity,
    counts_list = counts_list, n_resolution = n_resolution)
    
    return(out)
}



#' Plot GradCAM (or Heatmap) Values in Original Flow Cytometry Plot
#'
#' This function obtains GradCAM heatmap values for each cell in the original
#' flow cytometric space, mapping the `Angle` and `Intensity`
#' image dimensions (pixels) back to the corresponding cells.
#'
#' @param x A `TockyPrepData` object containing the original flow cytometry data.
#' @param feature_matrix A 100 x 100 matrix representing the GradCAM output.
#' @param xaxis The name of the dataframe column to be used as x-axis in the plot (default 'Red_log').
#' @param yaxis The name of the dataframe column to be used as y-axis in the plot (default 'Blue_log').
#' @param xlim Optional vector of length 2 defining the x-axis limits.
#' @param ylim Optional vector of length 2 defining the y-axis limits.
#' @param title Character string specifying the title of the plot.
#' @param color_bar Logical indicating whether to include a color bar in the plot (default TRUE).
#' @param n_resolution Binning resolution. The default is 100.
#' @param transpose Logical Whether to tranpose feature_matrix input. Note that TockyConvNetPy
#'  output Grad-CAM matrix for feature_matrix typically needs to be transposed. The default is TRUE.
#'
#' @return Invisibly returns the unmodified `TockyPrepData` object.
#' @export
#'
#' @examples
#' \dontrun{
#'   feature_matrix <- read.csv('heatmap.csv')
#'   par(mfrow= c(1,2))
#'   plotInverseGradCAM(x, feature_matrix, color_bar = TRUE)
#' }
#'
#' @importFrom dplyr filter mutate
#' @importFrom graphics rect plot segments text
#' @importFrom grDevices colorRampPalette
plotInverseGradCAM <- function(x, feature_matrix, xaxis = 'Red_log', yaxis = 'Blue_log',  xlim = NULL, ylim = NULL, title ='GradCAM', color_bar = TRUE, n_resolution = 100, transpose= TRUE){
    
    if(!inherits(x, 'TockyPrepData')){
        stop("Use a TockyPrepData object. \n")
        
    }
    
    if(transpose){
        feature_matrix <- as.matrix(t(feature_matrix))
        
    }
    expression_data <- x@Data
    expression_data <- expression_data[!is.na(expression_data$Angle),]
    expression_data$ID <- 1:nrow(expression_data)
    results <- convert_to_image(expression_data, n_resolution)
    bin_edges_Intensity <- results$bin_edges_Intensity
    bin_edges_Angle <- results$bin_edges_Angle
    
    row_idx <- findInterval(expression_data[['Angle']],     bin_edges_Angle)
    col_idx <- findInterval(expression_data[['Intensity']], bin_edges_Intensity)
    
    
    expression_data[['row_index']] <- row_idx
    expression_data[['col_index']] <- col_idx
    
    valid_mask <- row_idx > 0 & row_idx <= nrow(feature_matrix) &
    col_idx > 0 & col_idx <= ncol(feature_matrix)
    
    expression_data[['gradCAM_value']][valid_mask] <- feature_matrix[cbind(row_idx[valid_mask], col_idx[valid_mask])]
    
    coolwarm_colors <- c(
    rgb(59/255, 76/255, 192/255),   # Cool blue
    rgb(68/255, 142/255, 228/255),  # Lighter blue
    rgb(149/255, 165/255, 217/255), # Even lighter blue, approaching grey
    rgb(223/255, 229/255, 238/255), # Very light blue/grey - near white
    rgb(255/255, 255/255, 255/255), # White (central value)
    rgb(254/255, 224/255, 210/255), # Very light red/pink - near white
    rgb(252/255, 174/255, 145/255), # Light red
    rgb(251/255, 106/255, 74/255),  # Redder
    rgb(222/255, 45/255, 38/255),   # Deep red
    rgb(165/255, 0/255, 38/255)     # Dark red
    )
    
    heatmap_palette <- colorRampPalette(coolwarm_colors)(100)
    rng <- range(expression_data[['gradCAM_value']], na.rm = TRUE)
    norm_vals <- (expression_data[['gradCAM_value']] - rng[1]) / diff(rng)
    col_indices <- floor(norm_vals * 99) + 1  # Map to 1:100
    expression_data$cell_color <- heatmap_palette[col_indices]
    
    if(is.null(xlim)){
        xlim <- c(min(expression_data[[xaxis]]), max(expression_data[[xaxis]]))
    }
    
    if(is.null(ylim)){
        ylim <- c(min(expression_data[[yaxis]]), max(expression_data[[yaxis]]))
    }
    
    if(color_bar){
        layout_matrix <- matrix(c(1, 1, 1, 2), nrow=1, byrow=TRUE)
        layout(layout_matrix)
    }
    
    plot(expression_data[[xaxis]], expression_data[[yaxis]], col = expression_data$cell_color, pch = 16, xlim= xlim, ylim =ylim, main = title, xlab = xaxis, ylab = yaxis, cex.lab = 1.6, cex.main = 2)
    
    if(color_bar){
        plot(1, type = "n", xlab = "", ylab = "", xaxt = 'n', yaxt = 'n', bty = 'n', xlim = c(0, 1), ylim = c(0, 1))
        
        xleft <- 0.1
        xright <- 0.2
        ybottom <- 0.1
        ytop <- 0.9
        num_colors <- length(heatmap_palette)
        
        for (i in 1:num_colors) {
            rect(xleft, ybottom + (i - 1) / num_colors * (ytop - ybottom),
            xright, ybottom + i / num_colors * (ytop - ybottom),
            col = heatmap_palette[i], border = NA)
        }
        
        rect(xleft, ybottom, xright, ytop, col = NA, border = "black")
        
        num_ticks <- 5
        tick_positions <- seq(ybottom, ytop, length.out = num_ticks)
        label_values <- seq(rng[1], rng[2], length.out = num_ticks)
        
        for (i in 1:num_ticks) {
            segments(xright, tick_positions[i], xright + 0.01, tick_positions[i], col = "black")
            text(xright + 0.02, tick_positions[i], labels = format(label_values[i], digits = 2), pos = 4)
        }
        
    }
    return(invisible(x))
}



#' Generate a boxplot of MFI (median fluorescence intensity) for Grad-CAM feature cells,
#' other Timer+ cells, and Timer negative cells following inverseGradCAM.
#'
#' @param x A `TockyPrepData` object containing the original flow cytometry data.
#' @param feature_vector A vector output from inverseGradCAM.
#' @param group A character vector specifying the group(s) to use for analysis. If NULL, all groups are used.
#' @param select A logical indicating whether to allow interactive selection of variables for analysis.
#' @param variables A character vector specifying the variables to analyze. Only used if `select` is TRUE.
#'
#' @return A list containing the following elements:
#'   * `plot`: A ggplot object of the boxplot.
#'   * `summary_data`: A data frame containing the summarized data used to create the plot.
#'
#' @examples
#' \dontrun{
#'   feature_matrix <- read.csv('heatmap.csv')
#' plotGradCAMFeatureMFI(x, feature_matrix)
#' }
#' @export
#' @import ggplot2
#' @importFrom tidyr pivot_longer all_of
#' @importFrom dplyr group_by summarize %>%
#' @importFrom utils select.list
#' @importFrom rlang sym

plotGradCAMFeatureMFI  <-  function(x, feature_vector, group = NULL, select = FALSE, variables = NULL, title = 'GradCAM Feature Cells'){

    clusters <- ifelse(feature_vector==1, 'Feature', 'Others')
    data <- x@Data
    data <- data[!is.na(data$Angle),]
    data$Cluster <-  clusters
    
    nadata <- x@Data[is.na(x@Data$Angle),]
    
    
    data <- data[!is.na(data$Angle),]
    data <- as.data.frame(data)
    nadata$Cluster <- rep('Timer_Neg', nrow(nadata))
    data <- rbind(data, nadata)
    data <- merge(data, x@sampledef$sampledef, by = 'file')
    
    if(!is.null(group)){
        data <- data[data$group ==group,]
    }
    
    choices <- colnames(data)
    choices <- choices[grepl(pattern = 'logdata', choices)]
    
    if(select){
        choices <- sub(choices, pattern = '.logdata', replacement = '')
        var <- select.list(choices, graphics = TRUE, title = "Data to be analysed", multiple =TRUE)
        var_name <- paste0(var, '.logdata')
    }else{
        if(is.null(variables)){
            var_name <- choices
        }else{
            var_name <- variables
        }
        
    }
    
    data <- data[,c('file','Angle', 'Intensity', 'group', 'Cluster', var_name)]
    
    long_data <- pivot_longer(
    data,
    cols = all_of(var_name),
    names_to = "variable",
    values_to = "Expression"
    )
    
    if(!is.null(group)){
        long_data <- long_data[long_data$group %in% group, ]
    }
    
    long_data$Cluster <- as.factor(long_data$Cluster)
    
    df_long_all_summary <- long_data %>%
    group_by(file, !!sym('variable'), !!sym('Cluster')) %>%
    summarize(Mean = mean(!!sym('Expression'), na.rm = TRUE), .groups = 'drop')
    
    
    df_long_all_summary$variable <- sub(df_long_all_summary$variable, pattern = '.logdata', replacement ='')
    
    p <-   ggplot(df_long_all_summary, aes(x = !!sym('Cluster'), y = !!sym('Mean'),fill = !!sym('Cluster')))+
    geom_violin(trim = FALSE, alpha = 0.7) +
    geom_boxplot(width = 0.1, outlier.shape = NA, fill = "white", colour = "black", alpha = 0.5) +
    geom_jitter(width = 0.1, color = "black", size = 0.75, alpha = 0.6, shape = 21, show.legend = FALSE) +
    facet_wrap(~ variable, scales = "free_y", ncol = length(var_name)) +
    theme_bw() +
    labs(title = title,
    x = "Cluster",
    y = "Mean Fluorescence Intensity (MFI)")+
    theme(axis.text.x = element_text(angle = 45, hjust = 1))
    
    plot(p)
    out <- list(plot = p, summary_data = df_long_all_summary)
    return(invisible(out))
    
}



#' Generate Plots for Analysing Feature Cell Abundance
#'
#' This function processes clustering results, plots each cluster, and overlays each cluster's convex hull.
#' It is adaptable to any number of cell_cluster_id.
#'
#' @param res_tockyrf A list object output from `TockyKmeansRF`, which has been further processed
#'        using `ClusteringFeatureCells`.
#' @param p_adjust_method A method for p-value adjustment in multiple testing using Mann Whitney.
#' clusteringFeatureCells cen be used.
#' @param min_cells Numeric. The minimum nunmber of cells within a cluster to be analysed. The default is 10.
#' @param scatter_plot Logical. If TRUE, scatter plot for Angle and Intensity is generated.
#' @param ncol Number of columns in output figure panel.
#' @param Timer_positive Logical. Whether to remove Timer negative cells.
#' @param ylim Optional. A numeric vector of the length 2 for specifying ylim.
#' @importFrom graphics plot polygon points text
#' @importFrom grDevices rgb col2rgb rainbow
#' @importFrom grDevices chull
#' @importFrom stats wilcox.test p.adjust sd setNames
#' @importFrom ggplot2 ggplot geom_boxplot geom_violin labs geom_jitter theme_bw aes ylim theme element_text
#' @importFrom dplyr group_by summarise mutate ungroup select distinct left_join n
#' @importFrom magrittr %>%
#' @importFrom gridExtra grid.arrange
#' @importClassesFrom TockyPrep TockyPrepData
#'
#' @examples
#' \dontrun{
#'   data <- data.frame(Angle = runif(100), Intensity = runif(100))
#'   cell_cluster_id <- dbscan(data, eps = 0.1, minPts = 5)$cluster
#'   plotGradCAMFeatureCells(data, cell_cluster_id)
#' }
#' @export
plotGradCAMFeatureCells <- function(x, feature_cells, p_adjust_method = "BH", ncol = 3, min_cells = 10, scatter_plot = FALSE, title = 'GradCAM Feature Cells', Timer_positive = TRUE, ylim = NULL) {
    
    feature_cells <- as.numeric(feature_cells)
    clusters <- ifelse(feature_cells == 1, "Feature", 'Others')
    
    data <- x@Data
    if(Timer_positive){
        lg <- !is.na(data$Angle)
        data <- data[lg,]
        if(length(clusters) == nrow(x@Data)){
            clusters <- clusters[lg]
        }
        
    }
    sampledef <- x@sampledef$sampledef
    data$Cluster <- clusters
    
    data <- merge(data, sampledef, merge = 'file')
    
    data_percent <- data %>%
    group_by(file, Cluster, group) %>%
    summarise(Count = n(), .groups = 'drop') %>%
    group_by(file) %>%
    mutate(Total = sum(Count)) %>%
    ungroup() %>%
    mutate(Percentage = (Count / Total) * 100) %>%
    select(file, Cluster, group, Percentage)
    
    data_percent <- data_percent[data_percent$Cluster == "Feature",]
    
    summary_data <- data_percent %>%
    select(file, Cluster, group) %>%
    distinct() %>%
    left_join(data_percent, by = c("file", "Cluster", "group")) %>%
    group_by(group, Cluster) %>%
    summarise(
    Ave_Percentage = mean(Percentage, na.rm = TRUE),
    SD = sd(Percentage, na.rm = TRUE)
    )
    
    mw_test <- wilcox.test(Percentage ~ group, data = data_percent)
    
    
    
    p <- ggplot(data_percent, aes(x = group, y = Percentage, fill = group)) +
    geom_violin(trim = FALSE, alpha = 0.7, adjust = 1.5, width = 1.1)+
    geom_boxplot(width = 0.1, outlier.shape = NA, fill = "white", colour = "black", alpha = 0.5) +
    geom_jitter(width = 0.1, color = "black", size = 1.5, alpha = 0.6, shape = 21, show.legend = FALSE) +
    labs(title = title,
    subtitle = paste0("p-value: ", round(mw_test$p.value, 8)),
    y = "Percentage",
    x = "Group",
    fill = "Group") +
    theme_bw() +
    theme(axis.text.x = element_text(angle = 45, hjust = 1),
    legend.position = "right")  # E
    
    if(!is.null(ylim)){
        p <- p + ylim
    }
    
    plot(p)
    return(invisible(p))
}






    

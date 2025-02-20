#' Sample Splitting for CNN model Using Image Conversion
#' The function produces 100x100 image matrices, where the position (i, j) in the matrix corresponds to the count of points that fall into the i-th 'Angle' bin and j-th 'Intensity' bin.
#' @param x A TockyPrepData object
#' @param train_size The ratio of training/validation dataset: test_dataset (e.g. 0.8 for training/validation:test = 80:20).
#' @param valid_size The ratio of training dataset: valid dataset
#' @param test Whether to generate test data
#' @param n_resolution The number of bins to be used. The default is 100, i.e. the resolution of output images is 100 x 100.
#' @return An updated TockyPrepData object containing train, validation, and test datasets
#' @export
#' @examples
#' \dontrun{
#' x  <-  PrepTockyCNN(x)
#'}
#' @importFrom reticulate import array_reshape
#' @importFrom dplyr group_by ungroup filter n_distinct pull rowwise mutate summarise %>% .data
#' @importFrom utils read.table write.table
#' @importFrom methods show
#' @importClassesFrom TockyPrep TockyPrepData

PrepTockyCNN <- function(x, train_size = 0.8, test = TRUE, valid_size = 0.2, n_resolution = 100){
    if(!inherits(x, 'TockyPrepData')){
        stop("Use a TockyPrepData object. \n")
        
    }
    
    data <- x@Data
    sampledef <- x@sampledef$sampledef
    merged_data <- merge(data, sampledef, by = "file")
    
    if(test){
        tmp <- train_size
    }else{
        tmp <- 1
    }
    
    split_data <- function(data) {
        files_per_group <- data %>%
            dplyr::group_by(.data$group) %>%
            dplyr::summarise(num_files = n_distinct(.data$file), .groups = 'drop') %>%
            dplyr::mutate(num_train_valid_files = round(tmp * .data$num_files))
        
        train_valid_files <- files_per_group %>%
            dplyr::rowwise() %>%
            dplyr::mutate(sampled_files = list(sample(unique(data$file[data$group == .data$group]),
            size = .data$num_train_valid_files, replace = FALSE))) %>%
            dplyr::ungroup() %>%
            dplyr::pull(.data$sampled_files) %>%
            unlist()
        
        train_valid <- data %>%
            dplyr::filter(.data$file %in% train_valid_files)
        test_data <- data %>%
            dplyr::filter(!(.data$file %in% train_valid_files))
        
        list(train_valid = train_valid, test = test_data)
    }
    
    data_splits <- split_data(merged_data)
    data_train_valid <- data_splits$train_valid
    test_data <- data_splits$test
    
    
    cat("Number of cells in training/validation set:", nrow(data_train_valid), "\n")
    cat("Number of cells in testing set:", nrow(test_data), "\n")
    cat("Number of unique samples in training/validation set:", n_distinct(data_train_valid$file), "\n")
    cat("Number of unique samples in testing set:", n_distinct(test_data$file), "\n")
    
    print_distribution <- function(data, dataset_name) {
        cat("\nDistribution in", dataset_name, "\n")
        cat("-------------------------------\n")
        cat("Group:\n")
        print(table(data$group) / nrow(data) * n_resolution)
    }
    print_distribution(data_train_valid, "Training/Validation")
    
    if(test){
        print_distribution(test_data, "Test")
    }
    
    set.seed(123)
    unique_files_df <- unique(data_train_valid[, c("file", "group")])
    
    sample_train_files <- function(group_name) {
        sample(unique_files_df$file[unique_files_df$group == group_name],
        size = round((1-valid_size) * sum(unique_files_df$group == group_name)))
    }
    
    train_files <- unlist(lapply(unique(unique_files_df$group), sample_train_files))
    valid_files <- setdiff(unique_files_df$file, train_files)
    
    train_data <- data_train_valid[data_train_valid$file %in% train_files, ]
    valid_data <- data_train_valid[!data_train_valid$file %in% train_files, ]
    
    if(test){
        show(table(valid_data$group[duplicated(test_data$file) == FALSE]))
        cat("\n")
        cat("Test data: \n")
        show(table(test_data$group))
        cat("\n")
    }else{
        show(table(valid_data$group))
        
        cat("\n")
    }
    
    cat("Dimensions of training data:", dim(train_data), "\n")
    cat("Dimensions of validation data:", dim(valid_data), "\n")
    
    generate_label_matrix <- function(data) {
        unique_files <- unique(data$file)
        unique_groups <- unique(data$group)
        label_matrix <- matrix(0, nrow=length(unique_files), ncol=length(unique_groups))
        colnames(label_matrix) <- unique_groups
        rownames(label_matrix) <- unique_files
        
        for(i in 1:nrow(label_matrix)) {
            file <- rownames(label_matrix)[i]
            group_for_file <- data$group[data$file == file]
            label_matrix[i, group_for_file] <- 1
        }
        
        return(label_matrix)
    }
    
    train_labels <- generate_label_matrix(train_data)
    valid_labels <- generate_label_matrix(valid_data)
    cat("Sample numbers of training data: \n")
    show(colSums(train_labels))
    cat("Sample numbers of validation data: \n")
    show(colSums(valid_labels))
    test_labels <- NULL
    if (test) {
        test_labels <- generate_label_matrix(test_data)
        cat("Sample numbers of test data: \n")
        show(colSums(test_labels))
    }
    
    cn <- colnames(train_labels)
    train_labels <- train_labels[,cn]
    valid_labels <- valid_labels[,cn]
    
    if (test) {
        test_labels <- test_labels[,cn]
        
        out <- list(train_data = train_data,
        valid_data = valid_data,
        test_data = test_data,
        train_labels = train_labels,
        valid_labels = valid_labels,
        test_labels = test_labels)
    }else{
        out <- list(train_data = train_data,
        valid_data = valid_data,
        train_labels = train_labels,
        valid_labels = valid_labels)
    }
    
    np <- import('numpy')
    np$save('train_labels.npy', train_labels)
    np$save('valid_labels.npy', valid_labels)
    
    write.csv(train_labels, file = 'train_labels.csv')
    write.csv(valid_labels, file = 'valid_labels.csv')
    if (test) {
        np$save('test_labels.npy', test_labels)
        write.csv(test_labels, file='test_labels.csv')
    }
    
    x@Tocky[['PrepTockyCNN']] <- out
    
    return(invisible(x))
}



#' Normalize a matrix
#'
#' This function normalizes a matrix by multiplying each element by 10,000
#' and then dividing by the sum of the matrix.
#'
#' @param matrix A numeric matrix to be normalized.
#' @param n_resolution The number of bins to be used. The default is 100, i.e. the resolution of output images is 100 x 100.
#' @return A normalized matrix.
#' @examples
#' \dontrun{
#' matrix <- matrix(1:4, 2, 2)
#' normalize_matrix(matrix)
#' }
#' @keywords internal

normalize_matrix <- function(matrix, n_resolution = 100) {
    n_resolution * n_resolution * matrix / sum(matrix)
}

#' Normalize an array
#'
#' This function normalizes a 3-dimensional array by processing only the first layer,
#' multiplying each element by 10,000 and then dividing by the sum of the first layer.
#' It then replicates this normalized first layer if additional layers exist.
#'
#' @param array A 3-dimensional array to be normalized.
#' @param n_resolution The number of bins to be used. The default is 100, i.e. the resolution of output images is 100 x 100.
#' @return A normalized 3-dimensional array.
#' @examples
#' \dontrun{
#' array <- array(1:24, c(2, 2, 3))
#' normalize_array(array)
#' }
#' @keywords internal
#' @importFrom abind abind


normalize_array <- function(array, n_resolution = 100) {
    out <- n_resolution * n_resolution * array[, , 1] / sum(array[, , 1])
    out <- array(out, dim = c(dim(array)[1], dim(array)[2], 1))
    if (dim(array)[3] > 1) {
        out <- abind(out, array[, , 2:dim(array)[3]], along = 3)
    }
    
    return(out)
}

#' Reshape a list of matrices into an array
#'
#' This function reshapes a list of matrices into a 4-dimensional array, assuming
#' each matrix represents a single layer of an image or a similar data structure.
#'
#' @param data_list A list of numeric matrices to be reshaped.
#' @param n_resolution The number of bins to be used. The default is 100, i.e. the resolution of output images is 100 x 100.
#' @return A 4-dimensional array.
#' @examples
#' \dontrun{
#' data_list <- list(matrix(1:10000, 100, 100), matrix(101:20000, 100, 100))
#' reshape_data(data_list)
#' }
#' @keywords internal

reshape_data <- function(data_list, n_resolution = 100) {
  array_reshape(unlist(data_list), dim = c(length(data_list), n_resolution, n_resolution, 1))
}

#' Reshape a list of color matrices into an array
#'
#' This function reshapes a list of color matrices (each representing a different image)
#' into a 4-dimensional array, preserving the color channels.
#'
#' @param data_list A list of numeric matrices with color dimensions to be reshaped.
#' @return A 4-dimensional array.
#' @examples
#' \dontrun{
#' data_list <- list(array(1:30000, c(100, 100, 3)), array(301:60000, c(100, 100, 3)))
#' reshape_array_colour(data_list)
#' }
#' @export
#' @keywords internal

reshape_array_colour <- function(data_list) {
    num_images <- length(data_list)
    image_dims <- dim(data_list[[1]])  # Assuming all images have the same dimension
    
    result_array <- array(dim = c(num_images, image_dims[1], image_dims[2], image_dims[3]))
    
    for (i in seq_along(data_list)) {
         result_array[i, , , ] <- data_list[[i]]
    }
    
    return(result_array)
}



#' Image Conversion for Direct Export of Samples
#' @param x A TockyPrepData object
#' @param n_resolution The number of bins to be used. The default is 100, i.e. the resolution of output images is 100 x 100.
#' @param selected_markers A character vector with the length two for defining  marker data to be converted into 2D images.
#' @param output A character to specify the output directory.
#' @return An updated TockyPrepData object. In addition, the function exports image data as numpy array files, which are to be used for TockyMLPy analysis and TockyCNN model construction.
#' @export
#' @examples
#' \dontrun{
#' x <- ImageConversion(x, n_resolution = 100)
#'}
#' @importFrom reticulate array_reshape import
#' @importFrom utils write.csv


ImageConversion <- function(x, n_resolution = 100, selected_markers = NULL, output = NULL){
    
    if(!inherits(x, 'TockyPrepData')){
        stop("Use a TockyPrepData object. \n")
        
    }
    transformed_data <- x@Data
    
    if(is.null(selected_markers)){
        selected_markers <- c("Angle", "Intensity")
        transformed_data <- transformed_data[!is.na(transformed_data$Angle), ]
        
    }else{
        if(length(selected_markers)!=2){
            stop("Choose exactly two markers. \n")
            
        }
        
        if(!all(selected_markers %in% colnames(x@Data))){
            stop("Some marker data are not found in TockyPrepData. \n")
            
        }
        
    }
    
    sampledef <- x@sampledef$sampledef
    unique_files <- unique(transformed_data$file)
    sampledef <- sampledef[sampledef$file %in% unique_files, ]
    df <- sampledef
    unique_groups <- unique(sampledef$group)
    label_matrix <- matrix(0, nrow = length(unique_files), ncol = length(unique_groups),
    dimnames = list(unique_files, unique_groups))
    
    for (i in 1:nrow(df)) {
        file_row <- which(rownames(label_matrix) == df$file[i])
        group_col <- which(colnames(label_matrix) == df$group[i])
        label_matrix[file_row, group_col] <- 1
    }
    sample_labels <- label_matrix
    
    
    
    sample_out <- image_conversion_raw(transformed_data, selected_markers, n_resolution = n_resolution)
    
    normalize_matrix <- function(matrix) {
        n_resolution*n_resolution*matrix / sum(matrix)
    }
    
    sample_out$normalized_counts <- lapply(sample_out$counts, normalize_matrix)
    
    
    sample_images <- reshape_data(sample_out$normalized_counts, n_resolution = n_resolution)
    out <- list(sample_images = sample_images,
    sample_labels = sample_labels, selected_markers = sample_out$selected_markers)
    
    np <- import('numpy')
    
    if(!is.null(output) && length(output)==1){
        if (!dir.exists(output)) {
            dir.create(output)
        }
        
        npy_path <- file.path(output, 'sample_images.npy')
        np$save(npy_path, sample_images)
        npy_path <- file.path(output, 'sample_labels.npy')
        np$save(npy_path, sample_labels)
        csv_path <- file.path(output, 'sample_images.csv')
        write.csv(sample_images, file = csv_path)
        csv_path <- file.path(output, 'sample_labels.csv')
        write.csv(sample_labels, file = csv_path)
        
    }else{
        np$save('sample_images.npy', sample_images)
        np$save('sample_labels.npy', sample_labels)
        write.csv(sample_images, file = 'sample_images.csv')
        write.csv(sample_labels, file = 'sample_labels.csv')
        
    }
    
    x@Tocky[['TockyCNNimages']] <- out
    return(invisible(x))
    
}

#' Export TockyCNN Image Data from TockyPrepData Object
#'
#' Export CNN images data from a TockyPrepData object and prints diagnostic information about the image data structure.
#'
#' @param x A \code{TockyPrepData} object containing processed CNN images data. This object must have already undergone ImageConversion.
#' @param numpy Logical. Whether to export your data in the numpy format. If FALSE, csv files are exported instead.
#' @return The CNN images data stored in the object's \code{TockyCNNimages} slot. The returned value maintains the original structure from the slot, typically including:
#' \itemize{
#'   \item Image arrays in the first element
#'   \item Sample definitions in the second element
#'   \item Variables used in the third element
#' }
#'
#' @details This function performs the following actions:
#' \itemize{
#'   \item Checks if the \code{TockyCNNimages} slot contains data
#'   \item Prints diagnostic information including:
#'   \itemize{
#'     \item Dimensions of the image arrays
#'     \item Variables used in the analysis
#'     \item Sample definitions
#'   }
#' }
#'
#' @examples
#' \dontrun{
#' # After successful ImageConversion:
#' images_data <- getCNNimages(tocky_prep_object)
#' dim(images_data[[1]])  # Access image array dimensions
#' }
#'
#' @export

writeImages <- function(x, numpy = TRUE){
    
    if(length(x@Tocky[['TockyCNNimages']])==0){
        stop('Use ImageConversion \n.')
    }
    
    out <- x@Tocky[['TockyCNNimages']]

    
    if(numpy){
        np <- import('numpy')
        np$save('sample_images.npy',  out[['sample_images']])
        np$save('sample_labels.npy', out[['sample_labels']])
    }else{
        
        write.csv(out[['sample_images']], file = 'sample_images.csv')
        write.csv(out[['sample_labels']], file = 'sample_labels.csv')
    }
    
}





#' Get CNN Images from TockyPrepData Object
#'
#' Retrieves CNN images data from a TockyPrepData object and prints diagnostic information about the image data structure.
#'
#' @param x A \code{TockyPrepData} object containing processed CNN images data. This object must have already undergone ImageConversion.
#' @return The CNN images data stored in the object's \code{TockyCNNimages} slot. The returned value maintains the original structure from the slot, typically including:
#' \itemize{
#'   \item Image arrays in the first element
#'   \item Sample definitions in the second element
#'   \item Variables used in the third element
#' }
#'
#' @details This function performs the following actions:
#' \itemize{
#'   \item Checks if the \code{TockyCNNimages} slot contains data
#'   \item Prints diagnostic information including:
#'   \itemize{
#'     \item Dimensions of the image arrays
#'     \item Variables used in the analysis
#'     \item Sample definitions
#'   }
#' }
#'
#' @examples
#' \dontrun{
#' # After successful ImageConversion:
#' images_data <- getCNNimages(tocky_prep_object)
#' dim(images_data[[1]])  # Access image array dimensions
#' }
#'
#' @export

getImages <- function(x){
    
    if(length(x@Tocky[['TockyCNNimages']])==0){
        stop('Use ImageConversion \n.')
    }
    
    tmp <- x@Tocky[['TockyCNNimages']]
    
    cat(paste("Dimension of converted images:", dim(tmp), '\n'))
    
    cat(paste("Variables used:", tmp[[3]], '\n'))
    
    cat(paste("Sample definition:", tmp[[2]], '\n'))
    
    return(tmp)
    
    
}


#' Convert data into image
#' @param data A data frame
#' @param n_resolution The number of bins
#' @return A list object containing bin matrices list and counts matrices list
#' @export
#' @examples
#' \dontrun{
#' out <- image_conversion(data)
#'}
#' @keywords internal

image_conversion <- function(data, n_resolution = 100){
    long_matrix <- data[,c("Angle", "Intensity", "file")]
    files <- unique(long_matrix$file)
    matrices_list <- as.list(1:length(files))
    subset_data <- matrices_list

    for(i in 1:length(matrices_list)){
      subset_data[[i]] <- long_matrix[long_matrix$file == files[i],c("Angle", "Intensity")]
      subset_data[[i]] <- mutate(subset_data[[i]], ID = rownames(long_matrix)[long_matrix$file == files[i]])

    }

    matrices_list <- subset_data
    names(matrices_list) <- files
    global_min_angle <- min(sapply(matrices_list, function(df) min(df$Angle, na.rm = TRUE)), na.rm = TRUE)
    global_max_angle <- max(sapply(matrices_list, function(df) max(df$Angle, na.rm = TRUE)), na.rm = TRUE)

    global_min_intensity <- min(sapply(matrices_list, function(df) min(df$Intensity, na.rm = TRUE)), na.rm = TRUE)
    global_max_intensity <- max(sapply(matrices_list, function(df) max(df$Intensity, na.rm = TRUE)), na.rm = TRUE)

    bin_edges_Angle <- seq(from = global_min_angle, to = global_max_angle, length.out = n_resolution)
    bin_edges_Intensity <- seq(from = global_min_intensity, to = global_max_intensity, length.out = n_resolution)

    bin_matrices_list <- lapply(matrices_list, function(df) pixel_2d(df, bin_edges_Angle = bin_edges_Angle, bin_edges_Intensity = bin_edges_Intensity, n_resolution = n_resolution))

    counts_list <- lapply(bin_matrices_list, function(x) x$counts)
    
    out <- list(bin_matrices_list = bin_matrices_list,
    counts = counts_list)
    
    return(out)
}




#' Convert data into pixel-style data
#' @param data A data frame
#' @param bin_edges_Angle A numeric vector defining bin adges for Angle
#' @param bin_edges_Intensity A numeric vector defining bin adges for Angle
#' @param n_resolution The number of bins
#' @return A list object containing bin matrices list and counts matrices list
#' @keywords internal

pixel_2d <- function(data, bin_edges_Angle, bin_edges_Intensity, n_resolution = 100) {
    
    bin_indices <- .bincode2D(data$Intensity, data$Angle, bin1 = bin_edges_Intensity, bin2 = bin_edges_Angle)
    counts_matrix <- matrix(0, nrow = n_resolution, ncol = n_resolution)
    
    for (i in 1:length(bin_indices$bin1)) {
        col_idx <- bin_indices$bin1[i]#Intensity = cols
        row_idx <- bin_indices$bin2[i]#Angle = rows
        counts_matrix[row_idx, col_idx] <- counts_matrix[row_idx, col_idx] + 1
    }
    
    return(list(counts = counts_matrix))
    
}

#' Retain bin code data
#' @param x A numeric vector for which binning is performed for the number of bin1
#' @param bin1 The number of bin for x
#' @param y A numeric vector for which binning is performed for the number of bin2
#' @param bin2 The number of bin for y

#' @return A list object containing bin matrices list and counts matrices list
#' @keywords internal

.bincode2D <- function(x, y, bin1, bin2) {
    
    bin1_indices <- findInterval(x, bin1)
    bin2_indices <- findInterval(y, bin2)
    
    return(list(bin1 = bin1_indices, bin2 = bin2_indices))
}


#' Convert data into image
#' @param data A data frame
#' @param selected_markers A character vector with the length two for defining  marker data to be converted into 2D images.
#' @param n_resolution The number of bins to be used. The default is 100, i.e. the resolution of output images is 100 x 100.
#' @return A list object containing bin matrices list and counts matrices list
#' @examples
#' \dontrun{
#' out <- image_conversion_raw(data)
#'}
#' @keywords internal
image_conversion_raw <- function(data, selected_markers, n_resolution = 100){
    long_matrix <- data[,c(selected_markers, "file")]
    files <- unique(long_matrix$file)
    matrices_list <- as.list(1:length(files))
    subset_data <- matrices_list

    for(i in 1:length(matrices_list)){
      subset_data[[i]] <- long_matrix[long_matrix$file == files[i],c(selected_markers)]
      subset_data[[i]] <- mutate(subset_data[[i]], ID = rownames(long_matrix)[long_matrix$file == files[i]])

    }

    matrices_list <- subset_data
    names(matrices_list) <- files
    global_min_angle <- min(sapply(matrices_list, function(df) min(df[,selected_markers[1]], na.rm = TRUE)), na.rm = TRUE)
    global_max_angle <- max(sapply(matrices_list, function(df) max(df[,selected_markers[1]], na.rm = TRUE)), na.rm = TRUE)

    global_min_intensity <- min(sapply(matrices_list, function(df) min(df[,selected_markers[2]], na.rm = TRUE)), na.rm = TRUE)
    global_max_intensity <- max(sapply(matrices_list, function(df) max(df[,selected_markers[2]], na.rm = TRUE)), na.rm = TRUE)

    bin_edges_Angle <- seq(from = global_min_angle, to = global_max_angle, length.out = n_resolution)
    bin_edges_Intensity <- seq(from = global_min_intensity, to = global_max_intensity, length.out = n_resolution)

    bin_matrices_list <- lapply(matrices_list, function(df) pixel_2d_raw(df, bin_edges_Angle = bin_edges_Angle, bin_edges_Intensity = bin_edges_Intensity, selected_markers = selected_markers, n_resolution = n_resolution))

    counts_list <- lapply(bin_matrices_list, function(x) x$counts)
    
    out <- list(bin_matrices_list = bin_matrices_list,
    counts = counts_list, selected_markers= selected_markers, n_resolution = n_resolution)
    
    return(out)
}




#' Convert data into pixel-style data
#' @param data A data frame
#' @param bin_edges_Angle A numeric vector defining bin adges for Angle
#' @param bin_edges_Intensity A numeric vector defining bin adges for Angle
#' @param selected_markers A character vector with the length two for defining  marker data to be converted into 2D images.
#' @param n_resolution The number of bins
#' @return A list object containing bin matrices list and counts matrices list
#' @keywords internal
pixel_2d_raw <- function(data, bin_edges_Angle, bin_edges_Intensity, n_resolution = 100, selected_markers) {
    
    bin_indices <- .bincode2D(data[,selected_markers[2]], data[,selected_markers[1]], bin1 = bin_edges_Intensity, bin2 = bin_edges_Angle)
    counts_matrix <- matrix(0, nrow = n_resolution, ncol = n_resolution)
    
    for (i in 1:length(bin_indices$bin1)) {
        col_idx <- bin_indices$bin1[i]#Intensity = cols
        row_idx <- bin_indices$bin2[i]#Angle = rows
        counts_matrix[row_idx, col_idx] <- counts_matrix[row_idx, col_idx] + 1
    }
    
    return(list(counts = counts_matrix))
    
}



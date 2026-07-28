
library(tiff)
library(ggplot2)

set.seed(619)

norm1 <- readTIFF("ITIA1141.tif")
norm2 <- readTIFF("ITIA1142.tif")
norm3 <- readTIFF("ITIA1143.tif")
norm4 <- readTIFF("ITIA1144.tif")
norm5 <- readTIFF("ITIA1145.tif")

anom1 <- readTIFF("ITIA1101.tif")
anom2 <- readTIFF("ITIA1102.tif")
anom3 <- readTIFF("ITIA1103.tif")
anom4 <- readTIFF("ITIA1104.tif")
anom5 <- readTIFF("ITIA1105.tif")
anom6 <- readTIFF("ITIA1106.tif")
anom7 <- readTIFF("ITIA1107.tif")
anom8 <- readTIFF("ITIA1108.tif")
anom9 <- readTIFF("ITIA1109.tif")
anom10 <- readTIFF("ITIA1110.tif")
anom11 <- readTIFF("ITIA1111.tif")
anom12 <- readTIFF("ITIA1112.tif")
anom13 <- readTIFF("ITIA1113.tif")
anom14 <- readTIFF("ITIA1114.tif")
anom15 <- readTIFF("ITIA1115.tif")
anom16 <- readTIFF("ITIA1116.tif")
anom17 <- readTIFF("ITIA1117.tif")
anom18 <- readTIFF("ITIA1118.tif")
anom19 <- readTIFF("ITIA1119.tif")
anom20 <- readTIFF("ITIA1120.tif")
anom21 <- readTIFF("ITIA1121.tif")
anom22 <- readTIFF("ITIA1122.tif")
anom23 <- readTIFF("ITIA1123.tif")
anom24 <- readTIFF("ITIA1124.tif")
anom25 <- readTIFF("ITIA1125.tif")
anom26 <- readTIFF("ITIA1126.tif")
anom27 <- readTIFF("ITIA1127.tif")
anom28 <- readTIFF("ITIA1128.tif")
anom29 <- readTIFF("ITIA1129.tif")
anom30 <- readTIFF("ITIA1130.tif")
anom31 <- readTIFF("ITIA1131.tif")
anom32 <- readTIFF("ITIA1132.tif")
anom33 <- readTIFF("ITIA1133.tif")
anom34 <- readTIFF("ITIA1134.tif")
anom35 <- readTIFF("ITIA1135.tif")
anom36 <- readTIFF("ITIA1136.tif")
anom37 <- readTIFF("ITIA1137.tif")
anom38 <- readTIFF("ITIA1138.tif")
anom39 <- readTIFF("ITIA1139.tif")
anom40 <- readTIFF("ITIA1140.tif")

scan <- function(img, window_size) {
  rows <- nrow(img)
  cols <- ncol(img)
  max_so_far <- 0
  for (i in 1:(rows-window_size+1)) {
    for (j in 1:(cols-window_size+1)) {
      window_sum <- sum(img[i:(i+window_size-1), j:(j+window_size-1)])
      if (window_sum > max_so_far) {
        max_so_far <- window_sum
      }
    }
  }
  scan_stat <- max_so_far
  return(scan_stat)
}

get_null <- function(norm_list, window_size, N, block_rows, block_cols) {
  scan_stats <- c()
  stat_index <- 1
  
  for (i in 1:length(norm_list)) {
    img <- norm_list[[i]]
    n_rows <- nrow(img)
    n_cols <- ncol(img)
    
    # divide images into blocks
    n_row_blocks <- n_rows / block_rows
    n_col_blocks <- n_cols / block_cols
    total_blocks <- n_row_blocks * n_col_blocks
    
    # fills blocks matrix with top-left coordinates of each block
    blocks <- matrix(0, nrow = total_blocks, ncol = 2)
    idx <- 1
    for (r in 0:(n_row_blocks - 1)) {
      for (c in 0:(n_col_blocks - 1)) {
        blocks[idx, 1] <- r * block_rows + 1  
        blocks[idx, 2] <- c * block_cols + 1
        idx <- idx + 1
      }
    }
    
    # shuffle the block positions
    for (j in 1:N) {
      permuted_indices <- sample(1:total_blocks)
      permuted_img <- matrix(0, nrow = n_rows, ncol = n_cols)
      
    # move each original block to new location using the permuted indices  
      for (k in 1:total_blocks) {
        from <- blocks[k, ]
        to <- blocks[permuted_indices[k], ]
        
        block <- img[from[1]:(from[1] + block_rows - 1),
                     from[2]:(from[2] + block_cols - 1)]
        
        permuted_img[to[1]:(to[1] + block_rows - 1),
                     to[2]:(to[2] + block_cols - 1)] <- block
      }
      
      # calculate scan stat
      scan_stats[stat_index] <- scan(permuted_img, window_size)
      stat_index <- stat_index + 1
    }
  }
  
  return(scan_stats)
}

find_clusters <- function(anom_img, window_size, null_dist, alpha) {
  rows <- nrow(anom_img)
  cols <- ncol(anom_img)
  critical_value <- quantile(null_dist, 1 - alpha, na.rm = TRUE)
  clusters <- data.frame(row = integer(), col = integer(), scan_value = numeric(), p_val = numeric())
  
  for (i in 1:(rows - window_size + 1)) {
    for (j in 1:(cols - window_size + 1)) {
      window_sum <- sum(anom_img[i:(i + window_size - 1), j:(j + window_size - 1)])
      if (window_sum > critical_value) {
        p_val <- mean(null_dist >= window_sum, na.rm = TRUE)
        clusters <- rbind(clusters, data.frame(
          row = i,
          col = j,
          scan_value = window_sum,
          p_val = p_val
        ))
      }
    }
  }
  return(clusters)
}

norm_list <- list(norm1, norm2, norm3, norm4, norm5)

anom_list <- list(
  anom1, anom2, anom3, anom4, anom5, anom6, anom7, anom8, anom9, anom10,
  anom11, anom12, anom13, anom14, anom15, anom16, anom17, anom18, anom19, anom20,
  anom21, anom22, anom23, anom24, anom25, anom26, anom27, anom28, anom29, anom30,
  anom31, anom32, anom33, anom34, anom35, anom36, anom37, anom38, anom39, anom40
)

# make processed images

norm_values <- c()
for (i in 1:length(norm_list)) {
  img <- norm_list[[i]]
  norm_values <- c(norm_values, as.vector(img))
}

threshold <- mean(norm_values)

processed_norm_list <- list()
for (j in 1:length(norm_list)) {
  img <- norm_list[[j]]
  processed_img <- ifelse(img >= threshold, 1, 0)
  processed_norm_list[[j]] <- processed_img
}

processed_anom_list <- list()
for (k in 1:length(anom_list)) {
  img <- anom_list[[k]]
  processed_img <- ifelse(img >= threshold, 1, 0)
  processed_anom_list[[k]] <- processed_img
}

# calculate scan statistic for each anomalous image for each window size

N <- 200
block_rows <- 140
block_cols <- 128
anom_imgs <- anom_list

scan_stats <- matrix(NA, nrow = 17, ncol = length(anom_imgs))

for (i in 1:length(anom_imgs)) {
  window_size <- 10
  image <- anom_imgs[[i]]
  for (j in c(1:17)) {
    scan_stats[j, i] <- scan(image, window_size)
    window_size <- window_size + 5
  }
}

##  change name here
scan_stats_sq_org <- scan_stats

# find p values for the scan statistics for each window size for each image
 
window_size <- 10
norm_imgs <- norm_list

null_dists <- vector("list", 17)
p_vals <- matrix(NA, nrow = 17, ncol = length(anom_imgs))

for (k in c(1:17)) {
  null_dist <- get_null(norm_imgs, window_size, N, block_rows, block_cols)
  null_dists[[k]] <- null_dist
  for (l in 1:length(anom_imgs)) {
    obs_stat <- scan_stats[k, l]
    p_vals[k, l] <- mean(null_dist >= obs_stat)
  }
  window_size <- window_size + 5
}

## change name here
null_dists_sq_org_140128 <- null_dists
p_vals_sq_org_140128  <- p_vals

# find minimum p value and its corresponding window size for each image and each box size

window_sizes <- seq(20, 80, by = 10)

min_p_vals <- c()
best_window_sizes_list <- vector("list", length(anom_imgs))

for (m in 1:length(anom_imgs)) {
  col_vals <- p_vals_sq_org_140128[c(3,5,7,9,11,13,15), m]
  
  min_val <- min(col_vals)
  min_indices <- which(col_vals == min_val)
  
  best_window_sizes_list[[m]] <- window_sizes[min_indices]
  min_p_vals[m] <- min_val
}

best_window_sizes_sq_org_140128_new <- best_window_sizes_list

# find clusters at best window size for an image

img_index <- 40
img <- anom_list[[img_index]]
alpha <- significance_level_sq_org_1416_new
null_dists <- null_dists_sq_org_1416[c(3,5,7,9,11,13,15)]
window_sizes <- seq(20, 80, by = 10)

window_size_vec <- best_window_sizes_sq_org_1416_new[[img_index]]

pval_mask <- matrix(NA, nrow = nrow(img), ncol = ncol(img))

for (w in window_size_vec) {
  null_index <- which(window_sizes == w)
  null_dist <- null_dists[[null_index]]
  
  clusters <- find_clusters(img, w, null_dist, alpha)
  
  for (j in 1:nrow(clusters)) {
    m <- clusters$row[j]
    n <- clusters$col[j]
    p <- clusters$p_val[j]
    
    rows <- m:(m + w - 1)
    cols <- n:(n + w - 1)
    
    current_vals <- pval_mask[rows, cols]
    current_vals[is.na(current_vals)] <- Inf
    
    pval_mask[rows, cols] <- pmin(current_vals, p)
  }
}

img_df <- expand.grid(x = 1:ncol(img), y = 1:nrow(img))
img_df$value <- as.vector(t(img[nrow(img):1, ])) 
img_df$pval <- as.vector(t(pval_mask[nrow(pval_mask):1, ]))

ggplot(img_df, aes(x = x, y = y)) +
  geom_tile(aes(fill = value)) +
  scale_fill_gradient(low = "black", high = "white") +
  geom_tile(data = subset(img_df, !is.na(pval)),
            aes(alpha = 1 - pval),
            fill = "red") +
  scale_alpha_continuous(range = c(0.1, 1)) +
  scale_y_reverse() +
  coord_fixed() +
  theme_minimal() +
  guides(fill = "none", alpha = guide_legend(title = "1 - p-value"))


# power and distribution of minimum p values

p_vals_null <- c()
null_dists <- null_dists_sq_org_140128[c(3,5,7,9,11,13,15)]

for (k in c(1:1000)) {
  obs_p_vals <- c()
  for (l in c(1:7)) {
    obs_vals <- null_dists[[l]][k]
    obs_p_vals[l] <- mean(null_dists[[l]] >= obs_vals)
  }
  p_vals_null[k] <- min(obs_p_vals) 
}

significance_level_sq_org_140128_new <- quantile(p_vals_null, 0.05)

# calculate power of each permutation box size over all window sizes

p_vals <- p_vals_sq_org_140128
mean(p_vals <= significance_level)

# calculate power of each permutation box size for best window size

p_vals <- p_vals_sq_org_140128[c(3,5,7,9,11,13,15),]
window_sizes <- seq(20, 80, by = 10)

best_pvals <- c()

for (m in 1:ncol(p_vals)) {
  best_sizes <- best_window_sizes_sq_org_140128_new[[m]]
  best_indices <- which(window_sizes %in% best_sizes)
  tied_pvals <- p_vals[best_indices, m]
  best_pvals[m] <- min(tied_pvals)
}

mean(best_pvals <= significance_level_sq_org_140128_new)

## circle scan statistic

circle_scan <- function(img, radius) {
  rows <- nrow(img)
  cols <- ncol(img)
  max_so_far <- 0
  offset <- floor(radius)
  
  mask <- outer(-offset:offset, -offset:offset, function(x, y) sqrt(x^2 + y^2) <= radius)
  
  for (i in (offset + 1):(rows - offset)) {
    for (j in (offset + 1):(cols - offset)) {
      region <- img[(i - offset):(i + offset), (j - offset):(j + offset)]
      window_sum <- sum(region * mask)
      if (window_sum > max_so_far) {
        max_so_far <- window_sum
      }
    }
  }
  
  return(max_so_far)
}

get_circle_null <- function(norm_list, radius, N, block_rows, block_cols) {
  scan_stats <- c()
  stat_index <- 1
  
  for (i in 1:length(norm_list)) {
    img <- norm_list[[i]]
    n_rows <- nrow(img)
    n_cols <- ncol(img)
    n_row_blocks <- n_rows / block_rows
    n_col_blocks <- n_cols / block_cols
    total_blocks <- n_row_blocks * n_col_blocks
    
    blocks <- matrix(0, nrow = total_blocks, ncol = 2)
    idx <- 1
    for (r in 0:(n_row_blocks - 1)) {
      for (c in 0:(n_col_blocks - 1)) {
        blocks[idx, ] <- c(r * block_rows + 1, c * block_cols + 1)
        idx <- idx + 1
      }
    }
    
    for (j in 1:N) {
      permuted_indices <- sample(1:total_blocks)
      permuted_img <- matrix(0, nrow = n_rows, ncol = n_cols)
      
      for (k in 1:total_blocks) {
        from <- blocks[k, ]
        to <- blocks[permuted_indices[k], ]
        
        block <- img[from[1]:(from[1] + block_rows - 1),
                     from[2]:(from[2] + block_cols - 1)]
        
        permuted_img[to[1]:(to[1] + block_rows - 1),
                     to[2]:(to[2] + block_cols - 1)] <- block
      }
      
      scan_stats[stat_index] <- circle_scan(permuted_img, radius)
      stat_index <- stat_index + 1
    }
  }
  
  return(scan_stats)
}

find_circle_clusters <- function(anom_img, radius, null_dist, alpha) {
  rows <- nrow(anom_img)
  cols <- ncol(anom_img)
  critical_value <- quantile(null_dist, 1 - alpha, na.rm = TRUE)
  clusters <- data.frame(row = integer(), col = integer(), scan_value = numeric(), p_val = numeric())
  
  offset <- floor(radius)
  mask <- outer(-offset:offset, -offset:offset, function(x, y) sqrt(x^2 + y^2) <= radius)
  
  for (i in (offset + 1):(rows - offset)) {
    for (j in (offset + 1):(cols - offset)) {
      region <- anom_img[(i - offset):(i + offset), (j - offset):(j + offset)]
      window_sum <- sum(region * mask)
      if (window_sum > critical_value) {
        p_val <- mean(null_dist >= window_sum, na.rm = TRUE)
        clusters <- rbind(clusters, data.frame(row = i, col = j, scan_value = window_sum, p_val = p_val))
      }
    }
  }
  
  return(clusters)
}

# calculate circle scan stat at radius 10-30 seq 5

N <- 200
block_rows <- 14
block_cols <- 16
anom_imgs <- anom_list

scan_stats <- matrix(NA, nrow = 1, ncol = length(anom_imgs))

for (i in 1:length(anom_imgs)) {
  radius <- 28
  image <- anom_imgs[[i]]
  for (j in c(1:1)) {
    scan_stats[j, i] <- circle_scan(image, radius)
    radius <- radius + 5
  }
}

## CHANGE NAME HERE TO SAVE RESULT
scan_stats_crcl_one <- scan_stats

# find p values for the scan statistics for each window size for each image

radius <- 28
norm_imgs <- norm_list

null_dists <- vector("list", length = 1)
p_vals <- matrix(NA, nrow = 1, ncol = length(anom_imgs))

for (k in c(1:1)) {
  null_dist <- get_circle_null(norm_imgs, radius, N, block_rows, block_cols)
  null_dists[[k]] <- null_dist
  for (l in 1:length(anom_imgs)) {
    obs_stat <- scan_stats[k, l]
    p_vals[k, l] <- mean(null_dist >= obs_stat)
  }
  radius <- radius + 5
}

## CHANGE NAME HERE TO SAVE RESULT
null_dists_crcl_one_1416 <- null_dists
p_vals_crcl_one_1416 <- p_vals

# find minimum p value and its corresponding window size for each image

radii <- c(28)

min_p_vals <- c()
best_radii_list <- vector("list", length(anom_imgs))

for (m in 1:length(anom_imgs)) {
  col_vals <- p_vals_crcl_one_1416[, m]
  
  min_val <- min(col_vals)
  min_indices <- which(col_vals == min_val)
  
  best_radii_list[[m]] <- radii[min_indices]
  min_p_vals[m] <- min_val
}

## CHANGE NAME HERE TO SAVE RESULT
best_radii_crcl_one_1416 <- best_radii_list

# find clusters at best radii for an image

i <- 2
img <- anom_list[[i]]
alpha <- 0.05
null_dists <- null_dists_crcl_one_1416

# __x__ box permutation (change best_window_size name and null dists)
radii_vec <- best_radii_crcl_one_1416[[i]]
size_idx_limit <- length(radii_vec)

# 1 - size of window_size_vec (choose which best window size to look at)
size_idx <- 1

radius <- radii_vec[size_idx]
null_index <- which(radii == radius)

null_dist <- null_dists[[null_index]]
clusters <- find_circle_clusters(img, radius, null_dist, alpha)

# heatmap

img_df <- expand.grid(x = 1:ncol(img), y = 1:nrow(img))
img_df$value <- as.vector(t(img[nrow(img):1, ])) 

pval_mask <- matrix(NA, nrow = nrow(img), ncol = ncol(img))

offset <- floor(radius)
for (i in 1:nrow(clusters)) {
  m <- clusters$row[i]
  n <- clusters$col[i]
  p <- clusters$p_val[i]
  
  for (dx in -offset:offset) {
    for (dy in -offset:offset) {
      if (sqrt(dx^2 + dy^2) <= radius) {
        x <- m + dx
        y <- n + dy
        if (x >= 1 && x <= nrow(img) && y >= 1 && y <= ncol(img)) {
          if (is.na(pval_mask[x, y])) {
            pval_mask[x, y] <- p
          } 
          else {
            pval_mask[x, y] <- min(pval_mask[x, y], p, na.rm = TRUE)
          }
        }
      }
    }
  }
}


img_df$pval <- as.vector(t(pval_mask[nrow(pval_mask):1, ]))

ggplot(img_df, aes(x = x, y = y)) +
  geom_tile(aes(fill = value)) +
  scale_fill_gradient(low = "black", high = "white") +
  geom_tile(data = subset(img_df, !is.na(pval)),
            aes(fill = NULL, alpha = 1 - pval),
            fill = "red") +
  scale_alpha_continuous(range = c(0.1, 1)) +
  scale_y_reverse() +
  coord_fixed() +
  theme_minimal() +
  guides(fill = "none", alpha = guide_legend(title = "1 - p-value"))

# power calc

sum(p_vals_crcl_one_1416 <= 0.05) / length(anom_imgs)

## weight scan statistic

# calculate scan statistic for each anomalous image for each window size

# scan stat with weight functions

weighted_scan <- function(img, window_size) {
  rows <- nrow(img)
  cols <- ncol(img)
  epsilon = 1e-6
  max_so_far <- 0
  for (i in 1:(rows-window_size+1)) {
    for (j in 1:(cols-window_size+1)) {
      window <- img[i:(i + window_size - 1), j:(j + window_size - 1)]
      window_sum <- sum(window)
      window_var <- sd(as.vector(window))
      weighted_sum <- window_sum / (window_var + epsilon)
      if (weighted_sum > max_so_far) {
        max_so_far <- weighted_sum
      }
    }
  }
  scan_stat <- max_so_far
  return(scan_stat)
}

get_weighted_null <- function(norm_list, window_size, N, block_rows, block_cols) {
  scan_stats <- c()
  stat_index <- 1
  
  for (i in 1:length(norm_list)) {
    img <- norm_list[[i]]
    n_rows <- nrow(img)
    n_cols <- ncol(img)
    
    # divide images into blocks
    n_row_blocks <- n_rows / block_rows
    n_col_blocks <- n_cols / block_cols
    total_blocks <- n_row_blocks * n_col_blocks
    
    # fills blocks matrix with top-left coordinates of each block
    blocks <- matrix(0, nrow = total_blocks, ncol = 2)
    idx <- 1
    for (r in 0:(n_row_blocks - 1)) {
      for (c in 0:(n_col_blocks - 1)) {
        blocks[idx, 1] <- r * block_rows + 1  
        blocks[idx, 2] <- c * block_cols + 1
        idx <- idx + 1
      }
    }
    
    # shuffle the block positions
    for (j in 1:N) {
      permuted_indices <- sample(1:total_blocks)
      permuted_img <- matrix(0, nrow = n_rows, ncol = n_cols)
      
      # move each original block to new location using the permuted indices  
      for (k in 1:total_blocks) {
        from <- blocks[k, ]
        to <- blocks[permuted_indices[k], ]
        
        block <- img[from[1]:(from[1] + block_rows - 1),
                     from[2]:(from[2] + block_cols - 1)]
        
        permuted_img[to[1]:(to[1] + block_rows - 1),
                     to[2]:(to[2] + block_cols - 1)] <- block
      }
      
      # calculate scan stat
      scan_stats[stat_index] <- weighted_scan(permuted_img, window_size)
      stat_index <- stat_index + 1
    }
  }
  
  return(scan_stats)
}

find_weighted_clusters <- function(anom_img, window_size, null_dist, alpha) {
  epsilon <- 1e-6
  rows <- nrow(anom_img)
  cols <- ncol(anom_img)
  critical_value <- quantile(null_dist, 1 - alpha)
  
  clusters <- data.frame(row = integer(), col = integer(), scan_value = numeric(), p_val = numeric())
  
  for (i in 1:(rows - window_size + 1)) {
    for (j in 1:(cols - window_size + 1)) {
      window <- anom_img[i:(i + window_size - 1), j:(j + window_size - 1)]
      window_vec <- as.vector(window)
      
      window_sum <- sum(window_vec)
      window_var <- sd(window_vec)
      
      weighted_sum <- window_sum / (window_var + epsilon)
      
      if (weighted_sum > critical_value) {
        p_val <- mean(null_dist >= weighted_sum)
        clusters <- rbind(clusters, data.frame(
          row = i,
          col = j,
          scan_value = weighted_sum,
          p_val = p_val
        ))
      }
    }
  }
  
  return(clusters)
}

# list of images not detected by multiple window

anom_list_undetected <- list(
  anom4, anom15, anom16, anom20, anom23, anom28, anom31, anom32, anom35
)

# calculations

N <- 100
block_rows <- 14
block_cols <- 16
anom_imgs <- norm_list # here

scan_stats <- matrix(NA, nrow = 7, ncol = length(anom_imgs))

for (i in 1:length(anom_imgs)) {
  window_size <- 20
  image <- anom_imgs[[i]]
  for (j in c(1:7)) {
    scan_stats[j, i] <- weighted_scan(image, window_size)
    window_size <- window_size + 10
  }
}

##  change name here
scan_stats_norm <- scan_stats

# find p values for the scan statistics for each window size for each image

window_size <- 20
norm_imgs <- norm_list

null_dists <- vector("list", 7)
p_vals <- matrix(NA, nrow = 7, ncol = length(anom_imgs))

for (k in c(1:7)) {
  null_dist <- get_weighted_null(norm_imgs, window_size, N, block_rows, block_cols)
  null_dists[[k]] <- null_dist
  for (l in 1:length(anom_imgs)) {
    obs_stat <- scan_stats[k, l]
    p_vals[k, l] <- mean(null_dist >= obs_stat)
  }
  window_size <- window_size + 10
}

## change name here
null_dists_sq_weight_1416_sd <- null_dists
p_vals_sq_weight_1416_sd  <- p_vals

# find minimum p value and its corresponding window size for each image and each box size

window_sizes <- seq(20, 80, by = 10)

min_p_vals <- c()
best_window_sizes_list <- vector("list", length(anom_imgs))

for (m in 1:length(anom_imgs)) {
  col_vals <- p_vals_sq_weight_1416_sd[, m]
  
  min_val <- min(col_vals)
  min_indices <- which(col_vals == min_val)
  
  best_window_sizes_list[[m]] <- window_sizes[min_indices]
  min_p_vals[m] <- min_val
}

best_window_sizes_sq_weight_1416_sd <- best_window_sizes_list

# find clusters at best window size for an image

img_index <- 7
img <- anom_list[[img_index]]
alpha <- significance_level_sq_weight_1416_sd/100
null_dists <- null_dists_sq_weight_1416_sd
window_sizes <- seq(20, 80, by = 10)

window_size_vec <- best_window_sizes_sq_weight_1416_sd[[img_index]]

pval_mask <- matrix(NA, nrow = nrow(img), ncol = ncol(img))

for (w in window_size_vec) {
  null_index <- which(window_sizes == w)
  null_dist <- null_dists[[null_index]]
  
  clusters <- find_weighted_clusters(img, w, null_dist, alpha)
  
  for (j in 1:nrow(clusters)) {
    m <- clusters$row[j]
    n <- clusters$col[j]
    p <- clusters$p_val[j]
    
    rows <- m:(m + w - 1)
    cols <- n:(n + w - 1)
    
    current_vals <- pval_mask[rows, cols]
    current_vals[is.na(current_vals)] <- Inf
    
    pval_mask[rows, cols] <- pmin(current_vals, p)
  }
}

img_df <- expand.grid(x = 1:ncol(img), y = 1:nrow(img))
img_df$value <- as.vector(t(img[nrow(img):1, ])) 
img_df$pval <- as.vector(t(pval_mask[nrow(pval_mask):1, ]))

ggplot(img_df, aes(x = x, y = y)) +
  geom_tile(aes(fill = value)) +
  scale_fill_gradient(low = "black", high = "white") +
  geom_tile(data = subset(img_df, !is.na(pval)),
            aes(alpha = 1 - pval),
            fill = "red") +
  scale_alpha_continuous(range = c(0.1, 1)) +
  scale_y_reverse() +
  coord_fixed() +
  theme_minimal() +
  guides(fill = "none", alpha = guide_legend(title = "1 - p-value"))


# power and distribution of minimum p values

p_vals_null <- c()
null_dists <- null_dists_sq_weight_1416_sd

for (k in c(1:500)) {
  obs_p_vals <- c()
  for (l in c(1:7)) {
    obs_vals <- null_dists[[l]][k]
    obs_p_vals[l] <- mean(null_dists[[l]] >= obs_vals)
  }
  p_vals_null[k] <- min(obs_p_vals) 
}

significance_level_sq_weight_1416_sd <- quantile(p_vals_null, 0.05)

# calculate power of each permutation box size over all window sizes

p_vals <- p_vals_sq_weight_1416_sd
mean(p_vals <= significance_level_sq_weight_1416_sd)

# calculate power of each permutation box size for best window size

window_sizes <- seq(20, 80, by = 10)

best_pvals <- c()

for (m in 1:ncol(p_vals)) {
  best_sizes <- best_window_sizes_sq_weight_1416_sd[[m]]
  best_indices <- which(window_sizes %in% best_sizes)
  tied_pvals <- p_vals[best_indices, m]
  best_pvals[m] <- min(tied_pvals)
}

mean(best_pvals <= significance_level_sq_weight_1416_sd)
# 1

# false positives

false_positives <- matrix(FALSE, nrow = 7, ncol = ncol(scan_stats_norm))

for (j in 1:7) {
  crit_val <- quantile(null_dists_sq_weight_1416[[j]], (1-significance_level_sq_weight_1416_sd))
  false_positives[j, ] <- scan_stats_norm[j, ] > crit_val
}

sum(false_positives)

false_positive_rate <- sum(false_positives) / length(scan_stats_norm)

# single window scan stat

N <- 200
block_rows <- 14
block_cols <- 16
anom_imgs <- anom_list

scan_stats <- matrix(NA, nrow = 1, ncol = length(anom_imgs))

for (i in 1:length(anom_imgs)) {
  window_size <- 50
  image <- anom_imgs[[i]]
  for (j in c(1:1)) {
    scan_stats[j, i] <- scan(image, window_size)
  }
}

##  change name here
scan_stats_sq_one <- scan_stats

# find p values for the scan statistics for each window size for each image

window_size <- 50
norm_imgs <- norm_list

null_dists <- vector("list", 1)
p_vals <- matrix(NA, nrow = 1, ncol = length(anom_imgs))

for (k in c(1:1)) {
  null_dist <- get_null(norm_imgs, window_size, N, block_rows, block_cols)
  null_dists[[k]] <- null_dist
  for (l in 1:length(anom_imgs)) {
    obs_stat <- scan_stats[k, l]
    p_vals[k, l] <- mean(null_dist >= obs_stat)
  }
}

## change name here
null_dists_sq_one_1416 <- null_dists
p_vals_sq_one_1416  <- p_vals

# find minimum p value and its corresponding window size for each image and each box size

window_sizes <- c(50)

min_p_vals <- c()
best_window_sizes_list <- vector("list", length(anom_imgs))

for (m in 1:length(anom_imgs)) {
  col_vals <- p_vals_sq_one_1416[, m]
  
  min_val <- min(col_vals)
  min_indices <- which(col_vals == min_val)
  
  best_window_sizes_list[[m]] <- window_sizes[min_indices]
  min_p_vals[m] <- min_val
}

best_window_sizes_sq_one_1416 <- best_window_sizes_list

# find clusters at best window size for an image

img_index <- 18
img <- anom_list[[img_index]]
alpha <- 0.05
null_dists <- null_dists_sq_one_1416
window_sizes <- c(50)

window_size_vec <- best_window_sizes_sq_one_1416[[img_index]]

pval_mask <- matrix(NA, nrow = nrow(img), ncol = ncol(img))

for (w in window_size_vec) {
  null_index <- which(window_sizes == w)
  null_dist <- null_dists[[null_index]]
  
  clusters <- find_clusters(img, w, null_dist, alpha)
  
  for (j in 1:nrow(clusters)) {
    m <- clusters$row[j]
    n <- clusters$col[j]
    p <- clusters$p_val[j]
    
    rows <- m:(m + w - 1)
    cols <- n:(n + w - 1)
    
    current_vals <- pval_mask[rows, cols]
    current_vals[is.na(current_vals)] <- Inf
    
    pval_mask[rows, cols] <- pmin(current_vals, p)
  }
}

img_df <- expand.grid(x = 1:ncol(img), y = 1:nrow(img))
img_df$value <- as.vector(t(img[nrow(img):1, ])) 
img_df$pval <- as.vector(t(pval_mask[nrow(pval_mask):1, ]))

ggplot(img_df, aes(x = x, y = y)) +
  geom_tile(aes(fill = value)) +
  scale_fill_gradient(low = "black", high = "white") +
  geom_tile(data = subset(img_df, !is.na(pval)),
            aes(alpha = 1 - pval),
            fill = "red") +
  scale_alpha_continuous(range = c(0.1, 1)) +
  scale_y_reverse() +
  coord_fixed() +
  theme_minimal() +
  guides(fill = "none", alpha = guide_legend(title = "1 - p-value"))

# power calc

sum(p_vals_sq_one_1416 <= 0.05) / length(anom_imgs)

# KDE scan statistic

library(EBImage)

create_gaussian_kernel <- function(size, sigma) {
  coords <- seq(-(size - 1) / 2, (size - 1) / 2, length.out = size)
  kernel <- outer(coords, coords, function(x, y) exp(- (x^2 + y^2) / (2 * sigma^2)))
  kernel / sum(kernel)
}

# ---- Simple 2D Convolution ----
convolve2d <- function(img, kernel) {
  kr <- floor(nrow(kernel) / 2)
  kc <- floor(ncol(kernel) / 2)
  pad_img <- matrix(0, nrow = nrow(img) + 2 * kr, ncol = ncol(img) + 2 * kc)
  pad_img[(kr + 1):(nrow(pad_img) - kr), (kc + 1):(ncol(pad_img) - kc)] <- img
  result <- matrix(0, nrow = nrow(img), ncol = ncol(img))
  for (i in 1:nrow(img)) {
    for (j in 1:ncol(img)) {
      region <- pad_img[i:(i + 2 * kr), j:(j + 2 * kc)]
      result[i, j] <- sum(region * kernel)
    }
  }
  result
}

# ---- Build Pixelwise Permutation Null ----
get_null_kde_pixelwise <- function(norm_list, kernel, N, block_rows, block_cols) {
  img_dim <- dim(norm_list[[1]])
  n_rows <- img_dim[1]
  n_cols <- img_dim[2]
  total_perms <- N * length(norm_list)
  stack <- array(NA, dim = c(n_rows, n_cols, total_perms))
  idx <- 1
  for (img in norm_list) {
    n_row_blocks <- n_rows / block_rows
    n_col_blocks <- n_cols / block_cols
    total_blocks <- n_row_blocks * n_col_blocks
    blocks <- matrix(0, nrow = total_blocks, ncol = 2)
    bidx <- 1
    for (r in 0:(n_row_blocks - 1)) {
      for (c in 0:(n_col_blocks - 1)) {
        blocks[bidx, 1] <- r * block_rows + 1
        blocks[bidx, 2] <- c * block_cols + 1
        bidx <- bidx + 1
      }
    }
    for (perm in 1:N) {
      permuted_indices <- sample(1:total_blocks)
      permuted_img <- matrix(0, nrow = n_rows, ncol = n_cols)
      for (k in 1:total_blocks) {
        from <- blocks[k, ]
        to   <- blocks[permuted_indices[k], ]
        block <- img[
          from[1]:(from[1] + block_rows - 1),
          from[2]:(from[2] + block_cols - 1)
        ]
        permuted_img[
          to[1]:(to[1] + block_rows - 1),
          to[2]:(to[2] + block_cols - 1)
        ] <- block
      }
      stack[,,idx] <- convolve2d(permuted_img, kernel)
      idx <- idx + 1
    }
  }
  stack
}

# ---- Run KDE on an Image ----
kde_scan <- function(img, kernel) {
  kde_img <- convolve2d(img, kernel)
  kde_img
}

# ---- Pixelwise P-value & Cluster Masking ----
find_clusters_pixelwise <- function(anom_img, kernel, null_stack, alpha = 0.05) {
  kde_img <- kde_scan(anom_img, kernel)
  dims <- dim(kde_img)
  pval_img <- matrix(NA, nrow = dims[1], ncol = dims[2])
  for (i in 1:dims[1]) {
    for (j in 1:dims[2]) {
      null_dist <- null_stack[i, j, ]
      pval_img[i, j] <- mean(null_dist >= kde_img[i, j])
    }
  }
  mask <- (pval_img < alpha)
  cluster_labels <- NULL
  if ("EBImage" %in% rownames(installed.packages())) {
    cluster_labels <- EBImage::bwlabel(mask)
  }
  list(
    kde_img = kde_img,
    pval_img = pval_img,
    sig_mask = mask,
    clusters = cluster_labels
  )
}

kernel_size <- 91
sigma <- 15
kernel <- create_gaussian_kernel(kernel_size, sigma)
N <- 50            # Lower for tests, increase for bigger jobs!
block_rows <- 14
block_cols <- 16
alpha <- 0.05

# ---- BUILD NULL ----
cat("Building pixelwise permutation null... (this may take time)\n")
null_stack <- get_null_kde_pixelwise(norm_list, kernel, N, block_rows, block_cols)
cat("Null distribution done.\n")

# ---- ANALYZE ONE ANOMALOUS IMAGE ----
img_idx <- 1 # Choose which
cat(sprintf("Analyzing anomalous image #%d...\n", img_idx))
res <- find_clusters_pixelwise(norm_list[[img_idx]], kernel, null_stack, alpha)

# ---- VISUALIZATION ----
img <- norm_list[[img_idx]]
img_df <- expand.grid(x = 1:ncol(img), y = 1:nrow(img))
img_df$value <- as.vector(t(img[nrow(img):1, ]))
img_df$pval  <- as.vector(t(res$pval_img[nrow(res$pval_img):1, ]))
img_df$mask  <- as.vector(t(res$sig_mask[nrow(res$sig_mask):1, ]))

ggplot(img_df, aes(x = x, y = y)) +
  geom_tile(aes(fill = value)) +
  scale_fill_gradient(low = "black", high = "white") +
  geom_tile(data = subset(img_df, mask), fill = "red", alpha = 0.3) +
  scale_y_reverse() +
  coord_fixed() +
  theme_minimal()

## weighted new plots

find_weighted_clusters_zscore <- function(anom_img, window_size, null_dist, top_k = 5) {
  epsilon <- 1e-6
  rows <- nrow(anom_img)
  cols <- ncol(anom_img)
  
  null_mean <- mean(null_dist)
  null_sd <- sd(null_dist)
  
  clusters <- data.frame(row = integer(), col = integer(), scan_value = numeric(), z_score = numeric())
  
  for (i in 1:(rows - window_size + 1)) {
    for (j in 1:(cols - window_size + 1)) {
      window <- anom_img[i:(i + window_size - 1), j:(j + window_size - 1)]
      window_vec <- as.vector(window)
      
      window_sum <- sum(window_vec)
      window_var <- sd(window_vec)
      weighted_sum <- window_sum / (window_var + epsilon)
      
      z_score <- (weighted_sum - null_mean) / null_sd
      
      clusters <- rbind(clusters, data.frame(
        row = i,
        col = j,
        scan_value = weighted_sum,
        z_score = z_score
      ))
    }
  }
  
  # Return top K most extreme clusters
  clusters <- clusters[order(-clusters$z_score), ]
  return(head(clusters, top_k))
}

img_index <- 7
img <- anom_list[[img_index]]
zscore_mask <- matrix(NA, nrow = nrow(img), ncol = ncol(img))

for (w in window_size_vec) {
  null_index <- which(window_sizes == w)
  null_dist <- null_dists[[null_index]]
  
  clusters <- find_weighted_clusters_zscore(img, w, null_dist, top_k = 600)
  
  for (j in 1:nrow(clusters)) {
    m <- clusters$row[j]
    n <- clusters$col[j]
    z <- clusters$z_score[j]
    
    rows <- m:(m + w - 1)
    cols <- n:(n + w - 1)
    
    current_vals <- zscore_mask[rows, cols]
    current_vals[is.na(current_vals)] <- -Inf
    
    zscore_mask[rows, cols] <- pmax(current_vals, z)
  }
}

img_df <- expand.grid(x = 1:ncol(img), y = 1:nrow(img))
img_df$value <- as.vector(t(img[nrow(img):1, ])) 
img_df$zscore <- as.vector(t(zscore_mask[nrow(zscore_mask):1, ]))

ggplot(img_df, aes(x = x, y = y)) +
  geom_tile(aes(fill = value)) +
  scale_fill_gradient(low = "black", high = "white") +
  geom_tile(data = subset(img_df, !is.na(zscore)),
            aes(alpha = zscore),
            fill = "red") +
  scale_alpha_continuous(range = c(0.1, 1), name = "Z-score") +
  scale_y_reverse() +
  coord_fixed() +
  theme_minimal() +
  guides(fill = "none")

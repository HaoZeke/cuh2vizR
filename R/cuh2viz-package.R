## usethis namespace: start
#' @useDynLib cuh2vizR, .registration = TRUE
## usethis namespace: end
#' @export cuh2pot_single_con
#' @export cuh2pot_df
#' @export cuh2_pdat_con
#' @export cuh2_pdat_df
#' @export cuh2_scan_grid
#' @export get_energy
NULL
cuh2vizR_get_energy_df <- function(df,
                                   hcu_dists = seq(-0.05, 5, length.out = 60),
                                   hh_dists = seq(0.4, 3, length.out = 60)) {
  ## Initialize vectors to store the results
  energies <- vector("double", length(hcu_dists) * length(hh_dists))
  hcu_dists_vec <- vector("double", length(hcu_dists) * length(hh_dists))
  hh_dists_vec <- vector("double", length(hcu_dists) * length(hh_dists))
  ## Loop over the distances and compute the energies
  counter <- 1
  for (hcu_dist in hcu_dists) {
    for (hh_dist in hh_dists) {
      energies[counter] <- get_energy(dfCon$atom_data,
        hcu_dist = hcu_dist,
        hh_dist = hh_dist
      )
      hcu_dists_vec[counter] <- hcu_dist
      hh_dists_vec[counter] <- hh_dist
      counter <- counter + 1
    }
  }
  ## Create a data frame with the results
  df <- data.frame(
    hcu_dist = hcu_dists_vec,
    hh_dist = hh_dists_vec,
    energy = energies
  )
}

cuh2vizR_get_energy_contours <- function(df, clip_max = 5, clip_min = 0) {
  library("ggplot2")
  library("dplyr")
  df %>%
    mutate(energy = replace(energy, energy > 5, 5)) %>%
    mutate(energy = replace(energy, energy < 0, 0)) %>%
    ggplot(aes(x = hh_dist, y = hcu_dist, z = energy)) +
    ## Can also be geom_tile()
    geom_raster(interpolate = T, aes(fill = energy)) +
    geom_contour(color = "white") +
    scale_fill_gradientn(colors = hcl.colors(10, palette = "Blue-Red")) +
    labs(x = "H-H distance", y = "Cu-H2 distance")
}

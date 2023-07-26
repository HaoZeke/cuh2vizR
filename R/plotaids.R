#' @export cuh2vizR_get_energy_contours
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

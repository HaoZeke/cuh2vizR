#' Generate and plot energy contours
#'
#' This function takes a data frame and plots the energy contours
#' using ggplot2. The energy levels can be clipped with the `clip_max`
#' and `clip_min` parameters.
#'
#' @param df A data frame containing the columns 'hh_dist', 'hcu_dist' and 'energy'.
#' @param clip_max The maximum energy level. Default is 5.
#' @param clip_min The minimum energy level. Default is 0.
#' @return A ggplot object.
#' @export
#' @examples
#' \dontrun{
#' cuh2vizR_get_energy_contours(df, clip_max = 5, clip_min = 0)
#' }
cuh2vizR_get_energy_contours <- function(df, clip_max = 5, clip_min = 0) {
  library("ggplot2")
  library("dplyr")
  df %>%
    mutate(energy = replace(energy, energy > clip_max, clip_max)) %>%
    mutate(energy = replace(energy, energy < clip_min, clip_min)) %>%
    ggplot(aes(x = hh_dist, y = hcu_dist, z = energy)) +
    ## Can also be geom_tile()
    geom_raster(interpolate = T, aes(fill = energy)) +
    geom_contour(color = "white") +
    scale_fill_gradientn(colors = hcl.colors(10, palette = "Blue-Red")) +
    labs(x = "H-H distance", y = "Cu-H2 distance")
}

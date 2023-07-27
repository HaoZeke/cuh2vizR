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
cuh2vizR_get_energy_contours <- function(a_df, clip_max = 5, clip_min = 0) {
  library("ggplot2")
  library("dplyr")
  a_df %>%
    mutate(energy = replace(energy, energy > clip_max, clip_max)) %>%
    mutate(energy = replace(energy, energy < clip_min, clip_min)) %>%
    ggplot(aes(x = hh_dist, y = hcu_dist, z = energy)) +
    ## Can also be geom_tile()
    ggplot2::geom_raster(interpolate = T, aes(fill = energy)) +
    ggplot2::geom_contour(color = "white") +
    ggplot2::scale_fill_gradientn(colors = hcl.colors(10, palette = "Blue-Red")) +
    ggplot2::labs(x = "H-H distance", y = "Cu-H2 distance", title = "CuH2 Potential Energy (True)")
}


#' Generate and plot energy contours with animation
#'
#' This function takes a data frame and plots the energy contours
#' using ggplot2 and gganimate. The energy levels can be clipped with the `clip_max`
#' and `clip_min` parameters. The animation is saved as a GIF file.
#'
#' @param dfx A data frame containing the columns 'hh_dist', 'hcu_dist' and 'energy'.
#' @param df_list A list of data frames for each iteration.
#' @param clip_max The maximum energy level. Default is 5.
#' @param clip_min The minimum energy level. Default is 0.
#' @param filename The name of the output file for the animation.
#' @param anim_opts A list of other options for the animate function. Default is list().
#' @return NULL
#' @export
#' @examples
#' \dontrun{
#' cuh2vizR_generate_animation(
#' dfx,
#' df_list,
#' clip_max = 5,
#' clip_min = 0,
#' filename = "neb_cuh2.gif",
#' anim_opts = list(duration = 10, fps = 60, width = 800, height = 800)
#' )
#' }
cuh2vizR_generate_animation <- function(dfx, df_list, clip_max = 5, clip_min = 0, filename = "neb_cuh2.gif", anim_opts = list()) {
  # Adjust the energy values
  dfx <- dfx %>%
    mutate(energy = replace(energy, energy > clip_max, clip_max)) %>%
    mutate(energy = replace(energy, energy < clip_min, clip_min))

  # Define the ggplot
  p <- ggplot(dfx, aes(x = hh_dist, y = hcu_dist, z = energy)) +
    ggplot2::geom_raster(interpolate = T, aes(fill = energy)) +
    ggplot2::geom_contour(color = "white") +
    ggplot2::scale_fill_gradientn(colors = hcl.colors(10, palette = "Blue-Red")) +
    ggplot2::geom_point(data = df_list %>% bind_rows(), aes(x = hh_dist, y = hcu_dist), color = "black") +
    ggplot2::geom_line(data = df_list %>% bind_rows(), aes(x = hh_dist, y = hcu_dist), color = "black") +
    gganimate::transition_manual(iteration) +
    ggplot2labs(title = "Iteration: {current_frame}", x = "H-H distance", y = "Cu-H2 distance")

  # Generate the animation with the provided options
  gganimate::animate(p, renderer = gifski_renderer(file = filename), !!(anim_opts))

  return(NULL)
}

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

#' Generate energy dataframe
#'
#' This function generates a data frame that represents the energy
#' levels based on the distances specified by 'hcu_dists' and 'hh_dists'.
#' The energy levels are computed by the function 'get_energy'.
#'
#' @param df A data frame containing the atom data required by 'get_energy'.
#' @param hcu_dists A vector of H-Cu distances for which the energies are to be computed.
#' @param hh_dists A vector of H-H distances for which the energies are to be computed.
#' @return A data frame with columns 'hcu_dist', 'hh_dist' and 'energy'.
#' @export
#' @examples
#' \dontrun{
#' cuh2vizR_get_energy_df(df, hcu_dists = seq(-0.05, 5, length.out = 60), hh_dists = seq(0.4, 3, length.out = 60))
#' }
cuh2vizR_get_energy_df <- function(a_df,
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
      energies[counter] <- get_energy(a_dfCon$atom_data,
                                      hcu_dist = hcu_dist,
                                      hh_dist = hh_dist
      )
      hcu_dists_vec[counter] <- hcu_dist
      hh_dists_vec[counter] <- hh_dist
      counter <- counter + 1
    }
  }

  ## Create a data frame with the results
  a_df <- data.frame(
    hcu_dist = hcu_dists_vec,
    hh_dist = hh_dists_vec,
    energy = energies
  )
}

#' CuH2 visualization CLI
#'
#' This function provides a command-line interface for the CuH2 visualization
#' tools. It requires a directory path and a ground state configuration file.
#' Additionally, it accepts optional arguments for file pattern, range of file
#' indices to include, and animation parameters.
#'
#' The "ground state configuration" refers to a configuration file that describes
#' the initial state of the system, typically a reactant or product configuration.
#' This configuration is used as a base for generating the energy surface that
#' guides the animation.
#'
#' @param dir The directory path where the files are stored. [required]
#' @param ground The path to the file containing the base system configuration (ground state). [required]
#' @param pattern The file pattern to match in the directory. Default is "neb_path_%03d.con".
#' @param range The range of file indices to include. Default is c(0, 10).
#' @param filename The filename for the output animation. Default is "neb_cuh2.gif".
#' @param width The width of the output animation in pixels. Default is 800.
#' @param height The height of the output animation in pixels. Default is 800.
#' @param fps The frames per second for the output animation. Default is 60.
#' @param res The resolution for the output animation. Default is 150.
#' @param duration The duration of the output animation in seconds. Default is 10.
#' @return NULL
#' @export
#' @examples
#' cuh2vizR_cli(dir = "/path/to/files", ground = "/path/to/ground_file.con")
cuh2vizR_cli <- function() {
  option_list <- list(
    optparse::make_option(c("-d", "--dir"),
      type = "character", default = NULL,
      help = "Directory path [required]", metavar = "character"
    ),
    optparse::make_option(c("-g", "--ground"),
      type = "character", default = NULL,
      help = "Base system configuration [required]", metavar = "character"
    ),
    optparse::make_option(c("-p", "--pattern"),
      type = "character", default = "neb_path_%03d.con",
      help = "File pattern to match in the directory", metavar = "character"
    ),
    optparse::make_option(c("-r", "--range"),
      type = "integer", default = c(0, 10),
      help = "Range of file indices to include", metavar = "integer"
    ),
    optparse::make_option(c("-f", "--filename"),
      type = "character", default = "neb_cuh2.gif",
      help = "Filename for the output animation", metavar = "character"
    ),
    optparse::make_option("--width",
      type = "integer", default = 800,
      help = "Width of the output animation", metavar = "integer"
    ),
    optparse::make_option("--height",
      type = "integer", default = 800,
      help = "Height of the output animation", metavar = "integer"
    ),
    optparse::make_option("--fps",
      type = "integer", default = 60,
      help = "Frames per second for the output animation", metavar = "integer"
    ),
    optparse::make_option("--res",
      type = "integer", default = 150,
      help = "Resolution of the output animation", metavar = "integer"
    ),
    optparse::make_option("--duration",
      type = "integer", default = 10,
      help = "Duration of the output animation", metavar = "integer"
    )
  )

  opt_parser <- optparse::OptionParser(option_list = option_list)
  opt <- optparse::parse_args(opt_parser)

  # Check that the directory was provided
  if (is.null(opt$dir)) {
    stop("You must provide a directory path.")
  }

  # Check that an initial configuration was provided This is used to generate
  # the energy surface, it may be reactant, product, or any configuration
  if (is.null(opt$ground)) {
    stop("You must provide an initial system.")
  }

  # Generate list of file paths
  file_paths <- sprintf(
    file.path(opt$dir, opt$pattern),
    opt$range[1]:opt$range[2]
  )

  # Use lapply to create a list of data frames
  df_list <- lapply(file_paths, cuh2vizR::cuh2_pdat_con)

  # Add an 'iteration' column to each data frame
  for (i in seq_along(df_list)) {
    df_list[[i]] <- df_list[[i]] %>%
      dplyr::mutate(iteration = i)
  }


  dfCon <- readConR::readCon(opt$ground)
  cuh2_scan_grid(dfCon$atom_data, hcu_dists = seq(-0.05, 5.2, length.out = 60), seq(0.4, 3.3, length.out = 60)) -> dfx

  # Generate the animation
  cuh2vizR::cuh2vizR_generate_animation(
    dfx,
    df_list,
    clip_max = 5,
    clip_min = 0,
    filename = opt$filename,
    duration = opt$duration,
    fps = opt$fps,
    width = opt$width,
    height = opt$height,
    res = opt$res,
  )
}

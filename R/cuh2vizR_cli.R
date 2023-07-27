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
#' @param rangeidx A comma-separated string specifying the range of file indices to include. Example: "0,10". Default is "0,10".
#' @param filename The filename for the output animation, without the file type. Default is "neb_cuh2".
#' @param clip_max The maximum energy of the contour plot. Default is 5.
#' @param clip_min The minimum energy of the contour plot. Default is 0.
#' @param width The width of the output animation in pixels. Default is 800.
#' @param height The height of the output animation in pixels. Default is 800.
#' @param fps The frames per second for the output animation. Default is 60.
#' @param res The resolution for the output animation. Default is 150.
#' @param duration The duration of the output animation in seconds. Default is 10.
#' @param hh_range A comma-separated string specifying the range for H-H distances. Example: "0.4,3.2". Default is "0.4,3.2".
#' @param hcu_range A comma-separated string specifying the range for H-Cu distances. Example: "-0.05,5.1". Default is "-0.05,5.1".
#' @return Prints the success message and saves the generated animation to the specified file.
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
      help = "File pattern to match in the directory (e.g. neb_path_%03d.con)", metavar = "character"
    ),
    optparse::make_option(c("-r", "--rangeidx"),
      type = "character", default = "0,10",
      help = "Range of file indices to include (comma-separated)", metavar = "character"
    ),
    optparse::make_option(c("-f", "--filename"),
      type = "character", default = "neb_cuh2",
      help = "Filename for the output animation", metavar = "character"
    ),
    optparse::make_option("--clip_max",
      type = "numeric", default = 5,
      help = "Maximum energy of the contour plot", metavar = "numeric"
    ),
    optparse::make_option("--clip_min",
      type = "numeric", default = 0,
      help = "Minimum energy of the contour plot", metavar = "numeric"
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
    ),
    optparse::make_option("--hh_range",
      type = "character", default = "0.4,3.2",
      help = "Range for H-H distances (comma separated)", metavar = "character"
    ),
    optparse::make_option("--hcu_range",
      type = "character", default = "-0.05,5.1",
      help = "Range for H-Cu distances (comma separated)", metavar = "character"
    ),
    optparse::make_option("--num_points",
      type = "integer", default = 60,
      help = "Number of points for the scan grid", metavar = "integer"
    )
  )

  opt_parser <- optparse::OptionParser(option_list = option_list)
  opt <- optparse::parse_args(opt_parser)

  # Check that the directory was provided
  if (is.null(opt$dir)) {
    stop("You must provide a directory path.")
  } else {
    cli::cli_alert_success("Directory path: {opt$dir}")
  }

  # Check that an initial configuration was provided
  if (is.null(opt$ground)) {
    stop("You must provide an initial system.")
  } else {
    cli::cli_alert_success("Ground state configuration: {opt$ground}")
  }

  # Convert character ranges to numeric vectors
  rangeidx <- as.numeric(strsplit(opt$rangeidx, ",")[[1]])
  hh_range <- as.numeric(strsplit(opt$hh_range, ",")[[1]])
  hcu_range <- as.numeric(strsplit(opt$hcu_range, ",")[[1]])


  # Generate list of file paths
  file_paths <- sprintf(
    file.path(opt$dir, opt$pattern),
    rangeidx[1]:rangeidx[2]
  )

  cli::cli_alert_info("Processing Files")

  # Use lapply to create a list of data frames
  pb <- progress::progress_bar$new(
    format = "[:bar] :percent  Processing :filename",
    total = length(file_paths)
  )

  # Add an 'iteration' column to each data frame
  df_list <- lapply(seq_along(file_paths), function(i) {
    pb$tick(tokens = list(filename = file_paths[i]))
    df <- cuh2vizR::cuh2_pdat_con(file_paths[i])
    df %>%
      dplyr::mutate(iteration = i)
  })


  cli::cli_alert_info("Generating Animation")

  dfCon <- readConR::readCon(opt$ground)
  cuh2_scan_grid(dfCon$atom_data,
    hcu_dists = seq(hcu_range[1],
      hcu_range[2],
      length.out = opt$num_points
    ),
    hh_dists = seq(hh_range[1],
      hh_range[2],
      length.out = opt$num_points
    )
  ) -> dfx

  # Generate the animation
  cuh2vizR::cuh2vizR_generate_animation(
    dfx,
    df_list,
    startrange = rangeidx[1],
    clip_max = opt$clip_max,
    clip_min = opt$clip_min,
    filename = opt$filename,
    duration = opt$duration,
    fps = opt$fps,
    width = opt$width,
    height = opt$height,
    res = opt$res
  ) -> final_filename

  cli::cli_alert_success("Animation generated and saved as: {final_filename}")
}

#include "helpers.hpp"

void peturb_positions(rgpot::AtomMatrix &positions, Eigen::VectorXi &atmNumVec,
                      double hcu_dist, double hh_dist) {
  std::vector<int> hIndices, cuIndices;
  for (int i = 0; i < atmNumVec.size(); ++i) {
    if (atmNumVec[i] == 1) { // Hydrogen atom
      hIndices.push_back(i);
    } else if (atmNumVec[i] == 29) { // Copper atom
      cuIndices.push_back(i);
    } else {
      throw std::runtime_error("Unexpected atomic number");
    }
  }

  if (hIndices.size() != 2) {
    throw std::runtime_error("Expected exactly two hydrogen atoms");
  }

  std::cout << "H Coords before: " << positions.row(hIndices[0]) << " "
            << positions.row(hIndices[1]) << std::endl;

  // Compute the midpoint of the hydrogens
  Eigen::VectorXd hMidpoint =
      (positions.row(hIndices[0]) + positions.row(hIndices[1])) / 2;
  std::cout << "Midpoint of hydrogens: " << hMidpoint.transpose() << std::endl;

  // Compute the HH direction
  Eigen::VectorXd hh_direction;
  if (positions(hIndices[0], 0) < positions(hIndices[1], 0)) {
    hh_direction =
        (positions.row(hIndices[1]) - positions.row(hIndices[0])).normalized();
  } else {
    hh_direction =
        (positions.row(hIndices[0]) - positions.row(hIndices[1])).normalized();
  }
  std::cout << "HH direction: " << hh_direction.transpose() << std::endl;

  std::cout << "HH distance: " << hh_dist << std::endl;
  // Set the new position of the hydrogens
  positions.row(hIndices[0]) = hMidpoint - (0.5 * hh_dist) * hh_direction;
  positions.row(hIndices[1]) = hMidpoint + (0.5 * hh_dist) * hh_direction;

  std::cout << "H Coords after changing HH distance: "
            << positions.row(hIndices[0]) << " " << positions.row(hIndices[1])
            << std::endl;

  // Find the z-coordinate of the topmost Cu layer
  double maxCuZ = std::numeric_limits<double>::lowest();
  for (int cuIndex : cuIndices) {
    maxCuZ = std::max(maxCuZ, positions(cuIndex, 2));
  }

  // Compute the new z-coordinate for the H atoms
  double new_z = maxCuZ + hcu_dist;

  // Update the z-coordinates of the H atoms
  for (int hIndex : hIndices) {
    positions(hIndex, 2) = new_z;
  }

  std::cout << "maxCuZ: " << maxCuZ << std::endl;
  std::cout << "hcu_dist: " << hcu_dist << std::endl;
  std::cout << "new_z: " << new_z << std::endl;

  std::cout << "H Coords after moving towards the slab: "
            << positions.row(hIndices[0]) << " " << positions.row(hIndices[1])
            << std::endl;
}

[[cpp11::register]] double get_energy(const cpp11::data_frame &ref_df,
                                      double hh_dist, double hcu_dist) {
  // Get data from reference dataframe
  std::vector<double> xStd = cpp11::as_cpp<std::vector<double>>(ref_df["x"]);
  std::vector<double> yStd = cpp11::as_cpp<std::vector<double>>(ref_df["y"]);
  std::vector<double> zStd = cpp11::as_cpp<std::vector<double>>(ref_df["z"]);
  std::vector<int> atmNumStd =
      cpp11::as_cpp<std::vector<int>>(ref_df["atmNum"]);

  // Convert std::vectors to Eigen vectors
  Eigen::VectorXd xVec = Eigen::Map<Eigen::VectorXd>(xStd.data(), xStd.size());
  Eigen::VectorXd yVec = Eigen::Map<Eigen::VectorXd>(yStd.data(), yStd.size());
  Eigen::VectorXd zVec = Eigen::Map<Eigen::VectorXd>(zStd.data(), zStd.size());
  Eigen::VectorXi atmNumVec =
      Eigen::Map<Eigen::VectorXi>(atmNumStd.data(), atmNumStd.size());

  // Stack them into a matrix
  rgpot::AtomMatrix positions(xVec.size(), 3);
  positions << xVec, yVec, zVec;

  // Create a potential
  auto cuh2pot = rgpot::CuH2Pot();

  // Peturb the positions
  peturb_positions(positions, atmNumVec, hcu_dist, hh_dist);

  auto [energy, forces] =
      cuh2pot(positions, atmNumVec, cuh2vizR::constants::DEFAULT_BOX);
      energy -= cuh2vizR::constants::CUH2_GLOBAL_MIN;

  return energy;
}

[[cpp11::register]] cpp11::writable::data_frame
cuh2_scan_grid(const cpp11::data_frame &ref_df,
               const std::vector<double> &hcu_dists,
               const std::vector<double> &hh_dists) {
  std::vector<double> energyVector;
  std::vector<double> hcu_distVector;
  std::vector<double> hh_distVector;

  // Get data from reference dataframe
  std::vector<double> xStd = cpp11::as_cpp<std::vector<double>>(ref_df["x"]);
  std::vector<double> yStd = cpp11::as_cpp<std::vector<double>>(ref_df["y"]);
  std::vector<double> zStd = cpp11::as_cpp<std::vector<double>>(ref_df["z"]);
  std::vector<int> atmNumStd =
      cpp11::as_cpp<std::vector<int>>(ref_df["atmNum"]);

  // Convert std::vectors to Eigen vectors
  Eigen::VectorXd xVec = Eigen::Map<Eigen::VectorXd>(xStd.data(), xStd.size());
  Eigen::VectorXd yVec = Eigen::Map<Eigen::VectorXd>(yStd.data(), yStd.size());
  Eigen::VectorXd zVec = Eigen::Map<Eigen::VectorXd>(zStd.data(), zStd.size());
  Eigen::VectorXi atmNumVec =
      Eigen::Map<Eigen::VectorXi>(atmNumStd.data(), atmNumStd.size());

  // Create a potential
  auto cuh2pot = rgpot::CuH2Pot();
  for (double hcu_dist : hcu_dists) {
    for (double hh_dist : hh_dists) {
      // Stack them into a matrix
      rgpot::AtomMatrix positions(xVec.size(), 3);
      positions << xVec, yVec, zVec;

      // Here, adjust the positions of the atoms based on hcu_dist and hh_dist
      peturb_positions(positions, atmNumVec, hcu_dist, hh_dist);

      // Compute the energy for this configuration
      auto [energy, forces] =
          cuh2pot(positions, atmNumVec, cuh2vizR::constants::DEFAULT_BOX);
      energy -= cuh2vizR::constants::CUH2_GLOBAL_MIN;

      // Store the results
      energyVector.push_back(energy);
      hcu_distVector.push_back(hcu_dist);
      hh_distVector.push_back(hh_dist);
    }
  }

  // Convert std::vectors to cpp11::writable::doubles and create the dataframe
  using namespace cpp11::literals;
  cpp11::writable::data_frame result(
      {"energy"_nm =
           cpp11::writable::doubles(energyVector.begin(), energyVector.end()),
       "hcu_dist"_nm = cpp11::writable::doubles(hcu_distVector.begin(),
                                                hcu_distVector.end()),
       "hh_dist"_nm = cpp11::writable::doubles(hh_distVector.begin(),
                                               hh_distVector.end())});

  return result;
}

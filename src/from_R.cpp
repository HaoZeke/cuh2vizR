#include "helpers.hpp"

[[cpp11::register]]
cpp11::writable::list cuh2pot_df(const cpp11::data_frame& df) {
  // Convert dataframe columns to std::vectors then to Eigen vectors
  std::vector<double> xStd = cpp11::as_cpp<std::vector<double>>(df["x"]);
  std::vector<double> yStd = cpp11::as_cpp<std::vector<double>>(df["y"]);
  std::vector<double> zStd = cpp11::as_cpp<std::vector<double>>(df["z"]);
  std::vector<int> atmNumStd = cpp11::as_cpp<std::vector<int>>(df["atmNum"]);

  Eigen::VectorXd xVec = Eigen::Map<Eigen::VectorXd>(xStd.data(), xStd.size());
  Eigen::VectorXd yVec = Eigen::Map<Eigen::VectorXd>(yStd.data(), yStd.size());
  Eigen::VectorXd zVec = Eigen::Map<Eigen::VectorXd>(zStd.data(), zStd.size());
  Eigen::VectorXi atmNumVec = Eigen::Map<Eigen::VectorXi>(atmNumStd.data(), atmNumStd.size());

  const size_t framesize = xVec.size();

  // Stack them into a matrix
  rgpot::AtomMatrix positions(framesize, 3);
  positions << xVec, yVec, zVec;

  // Compute the energy
  auto cuh2pot = rgpot::CuH2Pot();
  auto [energy, forces] = cuh2pot(positions, atmNumVec, cuh2vizR::helpers::DEFAULT_BOX);

  // Calculate distances
  auto [hDistance, minCuDistance] = cuh2vizR::helpers::calculateDistances(positions, atmNumVec);

  // Return a named List with the energy, hDistance and minCuDistance
  cpp11::writable::list result;
  result.push_back("energy"_nm = cpp11::writable::doubles({energy}));
  result.push_back("hDistance"_nm = cpp11::writable::doubles({hDistance}));
  result.push_back("minCuDistance"_nm = cpp11::writable::doubles({minCuDistance}));

  return result;
}
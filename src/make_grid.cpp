#include "helpers.hpp"


[[cpp11::register]] cpp11::writable::data_frame
cuh2_scan_grid(const cpp11::data_frame& ref_df, const std::vector<double>& hcu_dists, const std::vector<double>& hh_dists) {
  std::vector<double> energyVector;
  std::vector<double> hcu_distVector;
  std::vector<double> hh_distVector;

  // Get data from reference dataframe
  std::vector<double> xStd = cpp11::as_cpp<std::vector<double>>(ref_df["x"]);
  std::vector<double> yStd = cpp11::as_cpp<std::vector<double>>(ref_df["y"]);
  std::vector<double> zStd = cpp11::as_cpp<std::vector<double>>(ref_df["z"]);
  std::vector<int> atmNumStd = cpp11::as_cpp<std::vector<int>>(ref_df["atmNum"]);

  // Convert std::vectors to Eigen vectors
  Eigen::VectorXd xVec = Eigen::Map<Eigen::VectorXd>(xStd.data(), xStd.size());
  Eigen::VectorXd yVec = Eigen::Map<Eigen::VectorXd>(yStd.data(), yStd.size());
  Eigen::VectorXd zVec = Eigen::Map<Eigen::VectorXd>(zStd.data(), zStd.size());
  Eigen::VectorXi atmNumVec = Eigen::Map<Eigen::VectorXi>(atmNumStd.data(), atmNumStd.size());

  // Stack them into a matrix
  rgpot::AtomMatrix positions(xVec.size(), 3);
  positions << xVec, yVec, zVec;

  // Create a potential
  auto cuh2pot = rgpot::CuH2Pot();
  for (double hcu_dist : hcu_dists) {
    for (double hh_dist : hh_dists) {
      // Here, adjust the positions of the atoms based on hcu_dist and hh_dist
      peturb_positions(positions, atmNumVec, hcu_dist, hh_dist);

      // Compute the energy for this configuration
      auto [energy, forces] = cuh2pot(positions, atmNumVec, cuh2vizR::helpers::DEFAULT_BOX);

      // Store the results
      energyVector.push_back(energy);
      hcu_distVector.push_back(hcu_dist);
      hh_distVector.push_back(hh_dist);
    }
  }

  // Convert std::vectors to cpp11::writable::doubles and create the dataframe
  using namespace cpp11::literals;
  cpp11::writable::data_frame result({
    "energy"_nm = cpp11::writable::doubles(energyVector.begin(), energyVector.end()),
    "hcu_dist"_nm = cpp11::writable::doubles(hcu_distVector.begin(), hcu_distVector.end()),
    "hh_dist"_nm = cpp11::writable::doubles(hh_distVector.begin(), hh_distVector.end())
  });

  return result;
}

#include "../subprojects/potlib/CppCore/src/CuH2/CuH2Pot.hpp"
#include "../subprojects/potlib/CppCore/src/base_types.hpp"
#include "../subprojects/potlib/CppCore/src/pot_types.hpp"
#include "../subprojects/readCon/CppCore/include/BaseTypes.hpp"
#include "../subprojects/readCon/CppCore/include/ReadCon.hpp"
#include "../subprojects/readCon/CppCore/include/helpers/StringHelpers.hpp"
#include "cpp11.hpp"
using namespace cpp11::literals;
using namespace cpp11;
using namespace yodecon::types;

// Define box (assuming it's constant for now)
const Eigen::Matrix3d BOX{{15.345599999999999, 0, 0},
                          {0, 21.702000000000002, 0},
                          {0, 0, 100.00000000000000}};

std::pair<double, double> calculateDistances(rgpot::AtomMatrix &positions,
                                             Eigen::VectorXi &atmtypes) {
  // Indices of Hydrogen atoms (Atomic number of Hydrogen is 1)
  std::vector<size_t> hIndices;
  // Positions of Copper atoms (Atomic number of Copper is 29)
  std::vector<Eigen::VectorXd> cuPositions;

  for (int i = 0; i < atmtypes.size(); ++i) {
    if (atmtypes[i] == 1) {
      hIndices.push_back(i);
    } else if (atmtypes[i] == 29) {
      cuPositions.push_back(positions.row(i));
    }
  }

  if (hIndices.size() != 2) {
    throw std::runtime_error("Expected exactly 2 Hydrogen atoms.");
  }

  // Calculate the distance between Hydrogen atoms
  double hDistance =
      (positions.row(hIndices[0]) - positions.row(hIndices[1])).norm();

  // Calculate the midpoint of Hydrogen atoms
  Eigen::VectorXd hMidpoint =
      (positions.row(hIndices[0]) + positions.row(hIndices[1])) / 2;

  // Calculate the minimum distance from the midpoint of H2 to each Copper atom
  double minCuDistance = std::numeric_limits<double>::max();
  for (const auto &cuPos : cuPositions) {
    double distance = (hMidpoint - cuPos).norm();
    if (distance < minCuDistance) {
      minCuDistance = distance;
    }
  }

  return std::make_pair(hDistance, minCuDistance);
}

[[cpp11::register]] writable::list cuh2pot_single_con(std::string fname) {
  std::vector<std::string> fconts =
      yodecon::helpers::file::read_con_file(fname);
  auto singleCon = yodecon::create_single_con<ConFrameVec>(fconts);

  // Populate atom positions matrix
  // Convert std::vectors to Eigen::VectorXd
  const size_t framesize = singleCon.x.size();
  Eigen::VectorXd xVec =
      Eigen::Map<Eigen::VectorXd>(singleCon.x.data(), framesize);
  Eigen::VectorXd yVec =
      Eigen::Map<Eigen::VectorXd>(singleCon.y.data(), framesize);
  Eigen::VectorXd zVec =
      Eigen::Map<Eigen::VectorXd>(singleCon.z.data(), framesize);
  std::vector<int> atmnums =
      yodecon::symbols_to_atomic_numbers(singleCon.symbol);
  // BUG: This shouldn't be necessary, and is a hacky fix
  atmnums.erase(std::remove(atmnums.begin(), atmnums.end(), 0), atmnums.end());
  Eigen::VectorXi atmtypes = Eigen::Map<VectorXi>((atmnums).data(), framesize);

  // Stack them into a matrix
  rgpot::AtomMatrix positions(framesize, 3);
  positions << xVec, yVec, zVec;

  // Compute the energy and forces
  auto cuh2pot = rgpot::CuH2Pot();
  auto [energy, forces] = cuh2pot(positions, atmtypes, BOX);

  // Prepare forces output matrix
  cpp11::writable::doubles_matrix<cpp11::by_row> forces_matrix(forces.rows(),
                                                               3);
  for (size_t idx{0}; idx < forces.rows(); ++idx) {
    if (!singleCon.is_fixed[idx]) {
      forces_matrix(idx, 0) = forces(idx, 0);
      forces_matrix(idx, 1) = forces(idx, 1);
      forces_matrix(idx, 2) = forces(idx, 2);
    } else {
      forces_matrix(idx, 0) = 0;
      forces_matrix(idx, 1) = 0;
      forces_matrix(idx, 2) = 0;
    }
  }

  // Return a named List with the energy and forces Matrix
  cpp11::writable::list result;
  result.push_back("energy"_nm = cpp11::writable::doubles({energy}));
  result.push_back("forces"_nm = forces_matrix);

  auto [hDistance, minCuDistance] = calculateDistances(positions, atmtypes);
  result.push_back("hDistance"_nm = cpp11::writable::doubles({hDistance}));
  result.push_back("minCuDistance"_nm =
                       cpp11::writable::doubles({minCuDistance}));

  return result;
}

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
  auto [energy, forces] = cuh2pot(positions, atmNumVec, BOX);

  // Calculate distances
  auto [hDistance, minCuDistance] = calculateDistances(positions, atmNumVec);

  // Return a named List with the energy, hDistance and minCuDistance
  cpp11::writable::list result;
  result.push_back("energy"_nm = cpp11::writable::doubles({energy}));
  result.push_back("hDistance"_nm = cpp11::writable::doubles({hDistance}));
  result.push_back("minCuDistance"_nm = cpp11::writable::doubles({minCuDistance}));

  return result;
}

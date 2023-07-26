#include "helpers.hpp"

namespace cuh2vizR::helpers {
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
} // namespace cuh2vizR::helpers

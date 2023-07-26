#include "helpers.hpp"

namespace cuh2vizR::helpers {
std::pair<double, double> calculateDistances(const rgpot::AtomMatrix &positions,
                                             Eigen::VectorXi &atmNumVec) {
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

  // Calculate the distance between Hydrogen atoms
  double hDistance =
      (positions.row(hIndices[0]) - positions.row(hIndices[1])).norm();

  // Calculate the midpoint of Hydrogen atoms
  Eigen::VectorXd hMidpoint =
      (positions.row(hIndices[0]) + positions.row(hIndices[1])) / 2;

  // Find the z-coordinate of the topmost Cu layer
  double maxCuZ = std::numeric_limits<double>::lowest();
  for (int cuIndex : cuIndices) {
    maxCuZ = std::max(maxCuZ, positions(cuIndex, 2));
  }

  double cuSlabDist = positions(hIndices[0], 2) - maxCuZ;

  return std::make_pair(hDistance, cuSlabDist);
}
} // namespace cuh2vizR::helpers

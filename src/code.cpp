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

  return result;
}

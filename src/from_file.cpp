#include "helpers.hpp"

[[cpp11::register]] cpp11::writable::list
cuh2pot_single_con(std::string fname) {
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
  auto [energy, forces] =
      cuh2pot(positions, atmtypes, cuh2vizR::helpers::DEFAULT_BOX);
  energy -= -697.311695; // 0 of the CuH2 system

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
  using namespace cpp11::literals;
  cpp11::writable::list result;
  result.push_back("energy"_nm = cpp11::writable::doubles({energy}));
  result.push_back("forces"_nm = forces_matrix);

  auto [hDistance, minCuDistance] =
      cuh2vizR::helpers::calculateDistances(positions, atmtypes);
  result.push_back("hDistance"_nm = cpp11::writable::doubles({hDistance}));
  result.push_back("minCuDistance"_nm =
                       cpp11::writable::doubles({minCuDistance}));

  return result;
}

[[cpp11::register]] cpp11::writable::data_frame
cuh2_pdat_con(std::string fname) {
  std::vector<std::string> fconts =
      yodecon::helpers::file::read_con_file(fname);
  std::vector<ConFrameVec> conList =
      yodecon::create_multi_con<ConFrameVec>(fconts);

  std::vector<double> energyVector;
  std::vector<double> hh_distVector;
  std::vector<double> hcu_distVector;

  auto cuh2pot = rgpot::CuH2Pot();
  for (auto &singleCon : conList) {
    // Populate atom positions matrix
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
    atmnums.erase(std::remove(atmnums.begin(), atmnums.end(), 0),
                  atmnums.end());
    Eigen::VectorXi atmtypes =
        Eigen::Map<VectorXi>((atmnums).data(), framesize);

    // Stack them into a matrix
    rgpot::AtomMatrix positions(framesize, 3);
    positions << xVec, yVec, zVec;

    // Compute the energy and forces
    auto [energy, forces] =
        cuh2pot(positions, atmtypes, cuh2vizR::helpers::DEFAULT_BOX);
    energy -= -697.311695; // 0 of the CuH2 system
    energyVector.push_back(energy);

    auto [hDistance, minCuDistance] =
        cuh2vizR::helpers::calculateDistances(positions, atmtypes);
    hh_distVector.push_back(hDistance);
    hcu_distVector.push_back(minCuDistance);
  }

  // Explicitly copy vectors into cpp11::writable::doubles
  cpp11::writable::doubles energyDoubles(energyVector.size());
  cpp11::writable::doubles hh_distDoubles(hh_distVector.size());
  cpp11::writable::doubles hcu_distDoubles(hcu_distVector.size());

  for (size_t i = 0; i < energyVector.size(); ++i) {
    energyDoubles[i] = energyVector[i];
    hh_distDoubles[i] = hh_distVector[i];
    hcu_distDoubles[i] = hcu_distVector[i];
  }

  using namespace cpp11::literals;
  cpp11::writable::data_frame result({"energy"_nm = energyDoubles,
                                      "hh_dist"_nm = hh_distDoubles,
                                      "hcu_dist"_nm = hcu_distDoubles});

  return result;
}

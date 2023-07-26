#include "../subprojects/potlib/CppCore/src/CuH2/CuH2Pot.hpp"
#include "../subprojects/potlib/CppCore/src/base_types.hpp"
#include "../subprojects/potlib/CppCore/src/pot_types.hpp"
#include "../subprojects/readCon/CppCore/include/BaseTypes.hpp"
#include "../subprojects/readCon/CppCore/include/ReadCon.hpp"
#include "../subprojects/readCon/CppCore/include/helpers/StringHelpers.hpp"
#include "cpp11.hpp"

using namespace yodecon::types;

namespace cuh2vizR::helpers {
// Define box (assuming it's constant for now)
const Eigen::Matrix3d DEFAULT_BOX{{15.345599999999999, 0, 0},
                                  {0, 21.702000000000002, 0},
                                  {0, 0, 100.00000000000000}};

std::pair<double, double> calculateDistances(rgpot::AtomMatrix &positions,
                                             Eigen::VectorXi &atmtypes);
}

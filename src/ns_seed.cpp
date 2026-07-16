/// @file ns_seed.cpp
/// @brief Implementation of random seed management utilities
#include "ns_seed.hpp"

/// @brief Anonymous namespace for internal state
namespace {
  /// @brief Global seed storage
  /// @details Default value is 42. Not protected by mutex - external
  ///          synchronization required for thread-safe access.
  unsigned global_seed = 42;
}

namespace seed {
  void set_global_seed(unsigned seed) {
    global_seed = seed;
  }

  unsigned get_global_seed() {
    return global_seed;
  }
}

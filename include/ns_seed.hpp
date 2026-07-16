/// @file ns_seed.hpp
/// @brief Random seed management utilities
/// @details
/// Provides a global seed management interface for random number generators.
/// The global seed can be set and retrieved for reproducible random sequences
/// across the application.
///
/// Thread-safety: Not thread-safe. Caller must ensure synchronization when
/// calling set_global_seed() from multiple threads.
#pragma once

namespace seed {
  /// @brief Set the global random seed
  /// @param[in] seed The seed value to set
  /// @note Not thread-safe. Ensure external synchronization if called concurrently.
  void set_global_seed(unsigned seed);

  /// @brief Get the current global random seed
  /// @return The current global seed value
  /// @note Not thread-safe. May return inconsistent values if set_global_seed()
  ///       is called concurrently.
  unsigned get_global_seed();
};

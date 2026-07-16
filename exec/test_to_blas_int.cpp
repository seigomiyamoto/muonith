#include <cstdint>
#include <iostream>
#include <limits>
#include <stdexcept>
#include <string>

#include "blas_backend.hpp"

int n_pass = 0;
int n_fail = 0;

void test_normal(const char* name, std::int64_t input, blas_int expected) {
  try {
    blas_int result = to_blas_int(input, name);
    if (result == expected) {
      std::cerr << "[PASS] " << name << ": to_blas_int(" << input << ") = " << result << "\n";
      n_pass++;
    } else {
      std::cerr << "[FAIL] " << name << ": to_blas_int(" << input << ") = " << result
                << ", expected " << expected << "\n";
      n_fail++;
    }
  } catch (const std::exception& e) {
    std::cerr << "[FAIL] " << name << ": unexpected exception: " << e.what() << "\n";
    n_fail++;
  }
}

void test_overflow(const char* name, std::int64_t input) {
  try {
    blas_int result = to_blas_int(input, name);
    std::cerr << "[FAIL] " << name << ": to_blas_int(" << input
              << ") returned " << result << ", expected std::overflow_error\n";
    n_fail++;
  } catch (const std::overflow_error& e) {
    std::cerr << "[PASS] " << name << ": caught overflow_error: " << e.what() << "\n";
    n_pass++;
  } catch (const std::exception& e) {
    std::cerr << "[FAIL] " << name << ": wrong exception type: " << e.what() << "\n";
    n_fail++;
  }
}

int main() {
  constexpr std::int64_t INT_MAX_64 = static_cast<std::int64_t>(std::numeric_limits<int>::max());
  constexpr std::int64_t INT_MIN_64 = static_cast<std::int64_t>(std::numeric_limits<int>::min());

  // Normal cases
  test_normal("positive",  42,         42);
  test_normal("zero",      0,          0);
  test_normal("INT_MAX",   INT_MAX_64, std::numeric_limits<int>::max());

  // Overflow / invalid cases
  test_overflow("negative",     -100);
  test_overflow("INT_MIN",      INT_MIN_64);
  test_overflow("INT_MAX+1",    INT_MAX_64 + 1);
  test_overflow("INT_MIN-1",    INT_MIN_64 - 1);
  test_overflow("INT64_MAX",    std::numeric_limits<std::int64_t>::max());

  // Summary
  std::cerr << "\n=== Summary: " << n_pass << " passed, " << n_fail << " failed ===\n";

  return (n_fail == 0) ? 0 : 1;
}

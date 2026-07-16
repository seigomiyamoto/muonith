/**
 * @file blas_backend.hpp
 * @brief Platform-agnostic BLAS/LAPACK backend selection header.
 *
 * @details Provides a unified interface for BLAS/LAPACK operations across
 * different backends:
 * - Linux/other: OpenBLAS (via cblas.h and lapacke.h)
 * - macOS: Apple Accelerate framework (with LAPACKE compatibility wrappers)
 *
 * Defines `blas_int` as the integer type for BLAS/LAPACK calls (LP64, 32-bit).
 *
 * @note All backends use LP64 (32-bit integer) interface.
 */
#pragma once

#include <cstdint>
#include <limits>
#include <stdexcept>
#include <string>

/// Integer type for BLAS/LAPACK calls (LP64: 32-bit)
using blas_int = int;

/**
 * @brief Cast a 64-bit integer to blas_int with overflow check.
 * @param[in] value Value to cast (e.g. Eigen::Index, size_t).
 * @param[in] context Caller name for error messages.
 * @return value as blas_int.
 * @throws std::overflow_error If value exceeds blas_int (int32) range.
 */
inline blas_int to_blas_int(std::int64_t value, const char* context) {
  if (value < 0) {
    throw std::overflow_error(
        std::string(context) + ": negative value " + std::to_string(value) +
        " is invalid for blas_int");
  }
  if (value > static_cast<std::int64_t>(std::numeric_limits<blas_int>::max())) {
    throw std::overflow_error(
        std::string(context) + ": value " + std::to_string(value) +
        " exceeds blas_int range (LP64 int32)");
  }
  return static_cast<blas_int>(value);
}

#if defined(__APPLE__)
  // --- Apple Accelerate Framework ---
  #define ACCELERATE_NEW_LAPACK
  #include <Accelerate/Accelerate.h>
  #include <vector>

  // Accelerate does not provide LAPACKE C interface.
  // Define thin wrappers around the Fortran-style API for compatibility.
  #define LAPACK_COL_MAJOR 102
  /**
   * @brief LAPACKE-compatible wrapper for sgetrf_ (float LU factorization).
   * @param[in] matrix_layout Must be LAPACK_COL_MAJOR (102).
   * @param[in] m Number of rows.
   * @param[in] n Number of columns.
   * @param[in,out] a Matrix to factorize (column-major, overwritten).
   * @param[in] lda Leading dimension of a.
   * @param[out] ipiv Pivot indices (size min(m,n)).
   * @return 0 on success, >0 if singular, <0 if invalid argument.
   */
  inline blas_int LAPACKE_sgetrf(int /*matrix_layout*/, blas_int m, blas_int n,
                                  float* a, blas_int lda, blas_int* ipiv) {
    blas_int info = 0;
    sgetrf_(&m, &n, a, &lda, ipiv, &info);
    return info;
  }

  /**
   * @brief LAPACKE-compatible wrapper for sgetri_ (float matrix inversion).
   * @param[in] matrix_layout Must be LAPACK_COL_MAJOR (102).
   * @param[in] n Order of the matrix.
   * @param[in,out] a LU-factorized matrix (overwritten with inverse).
   * @param[in] lda Leading dimension of a.
   * @param[in] ipiv Pivot indices from sgetrf.
   * @return 0 on success, >0 if singular, <0 if invalid argument.
   */
  inline blas_int LAPACKE_sgetri(int /*matrix_layout*/, blas_int n,
                                  float* a, blas_int lda, const blas_int* ipiv) {
    blas_int info = 0;
    blas_int lwork = -1;
    float work_query = 0.0f;
    // Workspace query
    sgetri_(&n, a, &lda, const_cast<blas_int*>(ipiv), &work_query, &lwork, &info);
    lwork = to_blas_int(static_cast<std::int64_t>(work_query), "LAPACKE_sgetri");
    std::vector<float> work(static_cast<size_t>(lwork));
    // Actual inversion
    sgetri_(&n, a, &lda, const_cast<blas_int*>(ipiv), work.data(), &lwork, &info);
    return info;
  }

  /**
   * @brief LAPACKE-compatible wrapper for dgetrf_ (double LU factorization).
   * @param[in] matrix_layout Must be LAPACK_COL_MAJOR (102).
   * @param[in] m Number of rows.
   * @param[in] n Number of columns.
   * @param[in,out] a Matrix to factorize (column-major, overwritten).
   * @param[in] lda Leading dimension of a.
   * @param[out] ipiv Pivot indices (size min(m,n)).
   * @return 0 on success, >0 if singular, <0 if invalid argument.
   */
  inline blas_int LAPACKE_dgetrf(int /*matrix_layout*/, blas_int m, blas_int n,
                                  double* a, blas_int lda, blas_int* ipiv) {
    blas_int info = 0;
    dgetrf_(&m, &n, a, &lda, ipiv, &info);
    return info;
  }

  /**
   * @brief LAPACKE-compatible wrapper for dgetri_ (double matrix inversion).
   * @param[in] matrix_layout Must be LAPACK_COL_MAJOR (102).
   * @param[in] n Order of the matrix.
   * @param[in,out] a LU-factorized matrix (overwritten with inverse).
   * @param[in] lda Leading dimension of a.
   * @param[in] ipiv Pivot indices from dgetrf.
   * @return 0 on success, >0 if singular, <0 if invalid argument.
   */
  inline blas_int LAPACKE_dgetri(int /*matrix_layout*/, blas_int n,
                                  double* a, blas_int lda, const blas_int* ipiv) {
    blas_int info = 0;
    blas_int lwork = -1;
    double work_query = 0.0;
    // Workspace query
    dgetri_(&n, a, &lda, const_cast<blas_int*>(ipiv), &work_query, &lwork, &info);
    lwork = to_blas_int(static_cast<std::int64_t>(work_query), "LAPACKE_dgetri");
    std::vector<double> work(static_cast<size_t>(lwork));
    // Actual inversion
    dgetri_(&n, a, &lda, const_cast<blas_int*>(ipiv), work.data(), &lwork, &info);
    return info;
  }

#else
  // --- OpenBLAS (Linux / other platforms) ---
  #include <cblas.h>
  #include <lapacke.h>
#endif

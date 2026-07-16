/// @file ns_eigen_blas.hpp
/// @brief Eigen and BLAS/LAPACK integration utilities for high-performance linear algebra operations
/// @details
/// This module provides a comprehensive bridge between Eigen linear algebra library and BLAS/LAPACK
/// optimized routines. It enables high-performance matrix operations,
/// conversions, and utilities for scientific computing applications.
///
/// **Workflow:**
/// - Matrix operations use optimized BLAS/LAPACK routines for improved performance
/// - Automatic thread management via OpenMP (functions marked with "Uses OpenMP")
/// - Conversion utilities between Eigen types and raw C arrays for BLAS/LAPACK interoperability
/// - I/O functions for matrix persistence and debugging
///
/// **Memory Layout:**
/// - All Eigen matrices must be ColMajor (Eigen's default) for BLAS/LAPACK compatibility
/// - Use check_colmajor_* functions to verify memory layout before BLAS/LAPACK calls
/// - Storage order: column-major (outer index = column, inner index = row)
///
/// **Thread-safety:**
/// - Functions using OpenMP are marked with "Uses OpenMP" and require external synchronization
/// - Thread count controlled via omp_get_max_threads() and set_threads_Eigen()
/// - Non-OpenMP functions are thread-safe for read-only operations
///
/// **Units:**
/// - Matrix/vector elements are dimensionless unless specified in application context
/// - Epsilon parameters use machine precision (std::numeric_limits<T>::epsilon())
///
/// **Coordinate System:**
/// - No specific coordinate system; application-dependent
#pragma once
#include <cstdio>
#include <string>
#include <cmath>
#include <iostream>
#include <fstream>
#include <chrono>
#include <thread>

#define _USE_MATH_DEFINES // for use M_PI, this order is required
#include <cmath>

#include "blas_backend.hpp"

// spdlog logger
#include "spdlog_pch.hpp"
#include "ns_mylogger.hpp"
#include "ns_type_definitions.hpp"

#include <Eigen/Dense>
#include <Eigen/Sparse>

/// @brief for type definitions
using namespace index_type_definitions;

/// @brief Type definition for sparse matrix
using SpMatf = Eigen::SparseMatrix<float>;

/// @brief Namespace for Eigen and BLAS/LAPACK related functions
namespace eigen_blas {

  //=================================================================
  /// @name blas_lapack_matrix 
  /// @details Functions for handling matrices with BLAS/LAPACK
  /// @{
  
  /// @brief Computes the inverse matrix using LU decomposition with float precision
  /// @param[in] mat Input square matrix (must be non-singular and ColMajor)
  /// @return Inverse of the input matrix
  /// @throws std::runtime_error If mat is not square (nrow!=ncol)
  /// @throws std::runtime_error If LU decomposition fails (matrix is singular)
  /// @throws std::runtime_error If inverse calculation fails
  /// @throws std::runtime_error If internal array comparison fails
  /// @note Memory layout: Input and output matrices are ColMajor
  /// @note Thread-safety: Uses omp_get_max_threads() but does not spawn threads
  /// @note Complexity: O(n^3) where n is the matrix dimension
  /// @note Uses OpenMP
  Eigen::MatrixXf calc_inverse_float( const Eigen::MatrixXf &mat );

  /// @brief Computes the inverse matrix using LAPACKE with float precision and alignment checks
  /// @param[in] inputMatrix Input square matrix (must be aligned and ColMajor)
  /// @return Inverse of the input matrix
  /// @throws std::runtime_error If inputMatrix data pointer is not aligned
  /// @throws std::runtime_error If inputMatrix data pointer is nullptr
  /// @throws std::runtime_error If matrix dimensions are non-positive (rows <= 0 or cols <= 0)
  /// @throws std::runtime_error If LU decomposition fails (LAPACKE_sgetrf returns non-zero)
  /// @throws std::runtime_error If inverse calculation fails (LAPACKE_sgetri returns non-zero)
  /// @note Memory layout: Input must be ColMajor with proper alignment (EIGEN_MAX_ALIGN_BYTES)
  /// @note Thread-safety: Thread-safe for concurrent calls with different matrices
  /// @note Complexity: O(n^3) where n is the matrix dimension
  Eigen::MatrixXf getInverseMatrixFloat(const Eigen::MatrixXf& inputMatrix);

  /// @brief Computes the inverse matrix using LAPACK with double precision, returning float result
  /// @param[in] mat Input square matrix (float precision, must be ColMajor)
  /// @return Inverse of the input matrix (converted back to float precision)
  /// @throws std::runtime_error If mat is not square (nrow!=ncol)
  /// @throws std::runtime_error If LU decomposition fails (LAPACKE_dgetrf returns non-zero)
  /// @throws std::runtime_error If inverse calculation fails (LAPACKE_dgetri returns non-zero)
  /// @note Memory layout: Input and output are ColMajor; internal computation uses double precision
  /// @note Thread-safety: Uses omp_get_max_threads() but does not spawn threads
  /// @note Complexity: O(n^3) where n is the matrix dimension
  /// @note Higher precision: Internally converts to double for computation, then back to float
  /// @note Uses OpenMP
  Eigen::MatrixXf getInverseMatrixDouble( const Eigen::MatrixXf &mat );

  /// @brief Computes the inverse matrix using LAPACK with double precision
  /// @param[in] mat Input square matrix (double precision, must be ColMajor)
  /// @return Inverse of the input matrix (double precision)
  /// @throws std::runtime_error If mat is not square (nrow!=ncol)
  /// @throws std::runtime_error If LU decomposition fails (LAPACKE_dgetrf returns non-zero)
  /// @throws std::runtime_error If inverse calculation fails (LAPACKE_dgetri returns non-zero)
  /// @note Memory layout: Input and output are ColMajor
  /// @note Thread-safety: Uses omp_get_max_threads() but does not spawn threads
  /// @note Complexity: O(n^3) where n is the matrix dimension
  /// @note Uses OpenMP
  Eigen::MatrixXd getInverseMatrixDouble( const Eigen::MatrixXd &mat );

  /// @brief Multiplies two matrices using cblas_sgemm with float precision
  /// @param[in] mat_A Left matrix (must be ColMajor)
  /// @param[in] mat_B Right matrix (must be ColMajor, cols_A must equal rows_B)
  /// @return Result matrix mat_C = mat_A * mat_B (ColMajor)
  /// @throws std::runtime_error If mat_A is not ColMajor
  /// @throws std::runtime_error If mat_B is not ColMajor
  /// @throws std::runtime_error If cols_A != rows_B (incompatible dimensions)
  /// @note Memory layout: All matrices (mat_A, mat_B, result) are ColMajor
  /// @note Thread-safety: Thread-safe for concurrent calls with different matrices
  /// @note Complexity: O(m*n*k) where mat_A is m×k and mat_B is k×n
  /// @note Performance: Uses optimized BLAS routine cblas_sgemm
  Eigen::MatrixXf multiplyMatrixColMajor(
    const Eigen::MatrixXf& mat_A, const Eigen::MatrixXf& mat_B);

  /// @brief Multiplies two matrices using cblas_dgemm with double precision
  /// @param[in] mat_A Left matrix (must be ColMajor)
  /// @param[in] mat_B Right matrix (must be ColMajor, cols_A must equal rows_B)
  /// @return Result matrix mat_C = mat_A * mat_B (ColMajor)
  /// @throws std::runtime_error If mat_A is not ColMajor
  /// @throws std::runtime_error If mat_B is not ColMajor
  /// @throws std::runtime_error If cols_A != rows_B (incompatible dimensions)
  /// @note Memory layout: All matrices (mat_A, mat_B, result) are ColMajor
  /// @note Thread-safety: Thread-safe for concurrent calls with different matrices
  /// @note Complexity: O(m*n*k) where mat_A is m×k and mat_B is k×n
  /// @note Performance: Uses optimized BLAS routine cblas_dgemm
  Eigen::MatrixXd multiplyMatrixColMajor(
    const Eigen::MatrixXd& mat_A, const Eigen::MatrixXd& mat_B);

  /// @brief Converts double-precision array to single-precision array
  /// @param[in] input Source array (double precision, must be non-null with at least size elements)
  /// @param[out] output Destination array (float precision, must be pre-allocated with at least size elements)
  /// @param[in] size Number of elements to convert (must be non-negative)
  /// @note Memory layout: Both arrays are contiguous in memory
  /// @note Thread-safety: Thread-safe if input and output do not overlap
  /// @note Complexity: O(size)
  /// @note Precision loss: Conversion from double to float may lose precision
  void convertDoubleToFloat(double* input, float* output, size_t size);

  /// @brief Outputs float C-array matrix to ASCII file in ColMajor order
  /// @param[in] mat Input matrix as C-array (ColMajor storage, must be non-null)
  /// @param[in] nrow Number of rows (must be positive)
  /// @param[in] ncol Number of columns (must be positive)
  /// @param[in] pathout Output file path
  /// @note Memory layout: Input is ColMajor (index = icol * nrow + irow)
  /// @note Thread-safety: Not thread-safe; external synchronization required for file I/O
  /// @note Complexity: O(nrow * ncol)
  /// @note Format: Each line contains "irow icol value" in scientific notation
  void out_matrix_array( const float *mat
  , const int nrow, const int ncol, const std::filesystem::path& pathout );

  /// @brief Outputs double C-array matrix to ASCII file in ColMajor order
  /// @param[in] mat Input matrix as C-array (ColMajor storage, must be non-null)
  /// @param[in] nrow Number of rows (must be positive)
  /// @param[in] ncol Number of columns (must be positive)
  /// @param[in] pathout Output file path
  /// @note Memory layout: Input is ColMajor (index = icol * nrow + irow)
  /// @note Thread-safety: Not thread-safe; external synchronization required for file I/O
  /// @note Complexity: O(nrow * ncol)
  /// @note Format: Each line contains "irow icol value" in scientific notation
  void out_matrix_array( const double *mat
  , const int nrow, const int ncol, const std::filesystem::path& pathout );
  
  /// @brief Outputs vector of matrices to ASCII file, skipping zero elements
  /// @param[in] pathout Output file path
  /// @param[in] vec_mat Vector of matrices to output
  /// @note Memory layout: Each matrix is ColMajor
  /// @note Thread-safety: Not thread-safe; external synchronization required for file I/O
  /// @note Complexity: O(total_elements) where total_elements is sum of all matrix sizes
  /// @note Format: Header "det_id = i nrow ncol" followed by non-zero elements "irow icol value"
  void out_vec_mat( const std::filesystem::path& pathout, const std::vector<Eigen::MatrixXf> &vec_mat );

  /// @brief Compares float C-array (ColMajor) with Eigen::MatrixXf element-wise
  /// @param[in] mat C-array matrix in ColMajor storage (must be non-null)
  /// @param[in] matxf Eigen matrix to compare against
  /// @return true if all elements match within tolerance (10.0 * epsilon), false otherwise
  /// @note Memory layout: Input array mat must be ColMajor (index = icol * nrow + irow)
  /// @note Thread-safety: Thread-safe for read-only comparison
  /// @note Complexity: O(nrow * ncol)
  /// @note Tolerance: Absolute difference must be <= 10.0 * std::numeric_limits<float>::epsilon()
  bool is_same_array_matrixXf( const float* mat, const Eigen::MatrixXf &matxf );

  /// @brief Compares double C-array (ColMajor) with Eigen::MatrixXd element-wise
  /// @param[in] mat C-array matrix in ColMajor storage (must be non-null)
  /// @param[in] matxd Eigen matrix to compare against
  /// @return true if all elements match within tolerance (10.0 * epsilon), false otherwise
  /// @note Memory layout: Input array mat must be ColMajor (index = icol * nrow + irow)
  /// @note Thread-safety: Thread-safe for read-only comparison
  /// @note Complexity: O(nrow * ncol)
  /// @note Tolerance: Absolute difference must be <= 10.0 * std::numeric_limits<double>::epsilon()
  bool is_same_array_matrixXd( const double* mat, const Eigen::MatrixXd &matxd );

  ///@} ------------------------------------------------------------------

  // ======================================================================
  /// @name Eigen_functions
  /// @details Functions related to Eigen
  /// @{

  /// @brief Validates Eigen::MatrixXf for non-empty dimensions and non-null data
  /// @param[in] mat Matrix to validate
  /// @return true if valid (rows > 0, cols > 0, data() != nullptr), false otherwise
  /// @note Thread-safety: Thread-safe for read-only validation
  /// @note Complexity: O(1)
  /// @note Side effects: Logs error messages if validation fails
  bool check_MatrixXf( const Eigen::MatrixXf &mat );

  /// @brief Computes element-wise division of two vectors with zero-handling
  /// @param[in] vec1 Numerator vector
  /// @param[in] vec2 Denominator vector (must have same size as vec1)
  /// @return Result vector where result(i) = vec1(i) / vec2(i), or -931931.0 if vec2(i) == 0.0
  /// @throws std::runtime_error If vec1.rows() != vec2.rows()
  /// @note Thread-safety: Thread-safe for concurrent calls with different vectors
  /// @note Complexity: O(n) where n is vector size
  /// @note Special value: Division by zero yields -931931.0 as sentinel value
  Eigen::VectorXf get_vecxf_div( const Eigen::VectorXf &vec1, const Eigen::VectorXf &vec2 );

  /// @brief Converts std::vector<double> to Eigen::VectorXf with precision downcast
  /// @param[in] vec Input vector (double precision)
  /// @return Output vector (single precision)
  /// @note Thread-safety: Thread-safe for concurrent calls with different vectors
  /// @note Complexity: O(n) where n is vec.size()
  /// @note Precision loss: Conversion from double to float may lose precision
  Eigen::VectorXf get_vecxf(const std::vector<double>& vec);


  /// @brief Reads sparse matrix data from ASCII file (format produced by out_matrix_non_zero)
  /// @param[in] path_in Input file path (must exist)
  /// @return Dense matrix with non-zero elements restored
  /// @throws std::runtime_error If path_in does not exist
  /// @throws std::runtime_error If file format is invalid (first line must be "# nrow ncol")
  /// @throws std::runtime_error If any data line does not have exactly 3 elements
  /// @note Thread-safety: Not thread-safe; external synchronization required for file I/O
  /// @note Complexity: O(nrow * ncol) for allocation + O(nnz) for reading non-zeros
  /// @note Format: First line "# nrow ncol", subsequent lines "irow icol value"
  Eigen::MatrixXf read_matrix_non_zero( const std::filesystem::path &path_in );

  /// @brief Computes row-wise sum to create vector from matrix
  /// @param[in] mat Input matrix
  /// @return Vector where each element is the sum of corresponding matrix row
  /// @note Thread-safety: Thread-safe for read-only operation
  /// @note Complexity: O(nrow * ncol)
  /// @note Precision: Returns double precision for increased accuracy of sums
  std::vector<double> mat_to_vec( const Eigen::MatrixXf &mat );

  /// @brief Computes matrix-vector multiplication with detailed per-element logging to file
  /// @param[in,out] fout Output file handle (must be open for writing)
  /// @param[in] matxf Matrix to multiply from the left
  /// @param[in] vecxf Vector to multiply from the right (size must equal matxf.cols())
  /// @return Result vector = matxf * vecxf
  /// @throws std::runtime_error If ncol != vecxf.rows()
  /// @throws std::runtime_error If result vector contains non-finite values
  /// @note Thread-safety: Not thread-safe; external synchronization required for file I/O
  /// @note Complexity: O(nrow * ncol)
  /// @note Format: Each line contains "irow icol factor_mat factor_vec product"
  /// @note Purpose: Debugging tool for verifying matrix-vector multiplication element-by-element
  Eigen::VectorXf check_mat_by_vec( FILE *fout
    , const Eigen::MatrixXf &matxf, const Eigen::VectorXf &vecxf );
  
  /// @brief Computes matrix-vector multiplication with detailed per-element logging to file
  /// @param[in] pathout Output file path
  /// @param[in] matxf Matrix to multiply from the left
  /// @param[in] vecxf Vector to multiply from the right (size must equal matxf.cols())
  /// @return Result vector = matxf * vecxf
  /// @throws std::runtime_error If ncol != vecxf.rows()
  /// @throws std::runtime_error If result vector contains non-finite values
  /// @note Thread-safety: Not thread-safe; external synchronization required for file I/O
  /// @note Complexity: O(nrow * ncol)
  /// @note Format: Each line contains "irow icol factor_mat factor_vec product"
  /// @note Purpose: Debugging tool for verifying matrix-vector multiplication element-by-element
  Eigen::VectorXf check_mat_by_vec( const std::filesystem::path& pathout
    , const Eigen::MatrixXf &matxf, const Eigen::VectorXf &vecxf );

  /// @brief Computes element-wise natural logarithm of vector using OpenMP
  /// @param[in] vec Input vector (all elements must be positive)
  /// @return Vector with natural logarithm of each element
  /// @throws std::runtime_error If any element of vec is <= 0.0
  /// @note Thread-safety: Uses OpenMP for parallelization; not thread-safe if called concurrently
  /// @note Complexity: O(n) with OpenMP parallelization
  /// @note Uses OpenMP: Parallelizes loop over vector elements
  Eigen::VectorXf mp_loge_vecxf( const Eigen::VectorXf &vec );

  /// @brief Checks if Eigen::VectorXf has contiguous storage (ColMajor equivalent for vectors)
  /// @param[in] vec Vector to check
  /// @return true if vec.innerStride() == 1 (contiguous), false otherwise
  /// @note Thread-safety: Thread-safe for read-only check
  /// @note Complexity: O(1)
  /// @note Memory layout: Returns true if elements are stored contiguously
  bool check_colmajor_vecxf(const Eigen::VectorXf& vec);

  /// @brief Checks if Eigen::VectorXd has contiguous storage (ColMajor equivalent for vectors)
  /// @param[in] vec Vector to check
  /// @return true if vec.innerStride() == 1 (contiguous), false otherwise
  /// @note Thread-safety: Thread-safe for read-only check
  /// @note Complexity: O(1)
  /// @note Memory layout: Returns true if elements are stored contiguously
  bool check_colmajor_vecxd(const Eigen::VectorXd& vec);

  /// @brief Validates that Eigen::MatrixXf is ColMajor and throws if not
  /// @param[in] mat Matrix to validate
  /// @throws std::runtime_error If mat is not ColMajor (outerStride != rows)
  /// @note Thread-safety: Thread-safe for read-only validation
  /// @note Complexity: O(1)
  /// @note Memory layout: Checks that outerStride == nrow (ColMajor requirement)
  void check_colmajor_matxf(const Eigen::MatrixXf& mat);

  /// @brief Validates that Eigen::MatrixXd is ColMajor and throws if not
  /// @param[in] mat Matrix to validate
  /// @throws std::runtime_error If mat is not ColMajor (outerStride != rows)
  /// @note Thread-safety: Thread-safe for read-only validation
  /// @note Complexity: O(1)
  /// @note Memory layout: Checks that outerStride == nrow (ColMajor requirement)
  void check_colmajor_matxd(const Eigen::MatrixXd& mat);  

  /// @brief Multiplies Eigen::MatrixXf and Eigen::VectorXf using cblas_sgemv with float precision
  /// @note mat, vec, and mat_ret must all be ColMajor
  /// ! Currently always causes Segmentation fault, so do not use
  // Eigen::VectorXf multiplyMatrixToVector(
  //   const Eigen::MatrixXf& mat, const Eigen::VectorXf& vec);

  /// @brief Multiplies Eigen::MatrixXd and Eigen::VectorXd using cblas_dgemv with double precision
  /// @note mat, vec, and mat_ret must all be ColMajor
  /// ! Currently always causes Segmentation fault, so do not use
  // Eigen::VectorXd multiplyMatrixToVector(
  //   const Eigen::MatrixXd& mat, const Eigen::VectorXd& vec);

  /// @brief Concatenates two matrices vertically (along rows)
  /// @param[in] mat1 Top matrix (must have positive rows and same cols as mat2)
  /// @param[in] mat2 Bottom matrix (must have positive rows and same cols as mat1)
  /// @return Merged matrix with mat1 on top and mat2 on bottom
  /// @throws std::runtime_error If mat1.cols() != mat2.cols()
  /// @throws std::runtime_error If mat1.rows() <= 0 or mat2.rows() <= 0
  /// @note Thread-safety: Thread-safe for concurrent calls with different matrices
  /// @note Complexity: O((rows1 + rows2) * cols)
  /// @note Memory layout: Result is ColMajor
  Eigen::MatrixXf get_merged_matrix_rows(
    const Eigen::MatrixXf& mat1, const Eigen::MatrixXf& mat2);

  /// @brief Appends matrix rows in-place (vertical concatenation)
  /// @param[in,out] mat_large Target matrix to expand (modified in-place)
  /// @param[in] mat_tobe_added Matrix to append at bottom (must have positive rows)
  /// @throws std::runtime_error If mat_large.cols() != mat_tobe_added.cols() (when mat_large is non-empty)
  /// @throws std::runtime_error If mat_tobe_added.rows() <= 0
  /// @note Thread-safety: Not thread-safe; modifies mat_large in-place
  /// @note Complexity: O((rows_org + rows_add) * cols) due to memory reallocation
  /// @note Special case: If mat_large is empty (rows == 0 or cols == 0), simply copies mat_tobe_added
  /// @note Memory layout: Result is ColMajor
  void append_matrix_rows(
    Eigen::MatrixXf& mat_large, const Eigen::MatrixXf& mat_tobe_added);

  /// @brief Logs dimensions of each matrix in vector to logger
  /// @param[in] vec_mat Vector of matrices to display
  /// @param[in] level Logging level (default: spdlog::level::debug)
  /// @note Thread-safety: Thread-safe if logger is thread-safe
  /// @note Complexity: O(vec_mat.size())
  /// @note Format: Logs "vec_mat[i].rows()=X, vec_mat[i].cols()=Y" for each matrix
  void disp_vec_matxf(
    const std::vector<Eigen::MatrixXf>& vec_mat
  , const spdlog::level::level_enum& level = spdlog::level::debug);

  /// @brief Logs dimensions of each sparse matrix in vector to logger
  /// @param[in] vec_mat Vector of sparse matrices to display
  /// @param[in] level Logging level (default: spdlog::level::debug)
  /// @note Thread-safety: Thread-safe if logger is thread-safe
  /// @note Complexity: O(vec_mat.size())
  /// @note Format: Logs "vec_mat[i].rows()=X, vec_mat[i].cols()=Y" for each sparse matrix
  void disp_vec_spmatf(
    const std::vector<SpMatf>& vec_mat
  , const spdlog::level::level_enum& level = spdlog::level::debug);

  /// @brief Compares two vectors of matrices for approximate equality
  /// @param[in] vec_mat1 First vector of matrices
  /// @param[in] vec_mat2 Second vector of matrices
  /// @return true if vectors have same size and all corresponding matrices are approximately equal, false otherwise
  /// @note Thread-safety: Thread-safe for read-only comparison
  /// @note Complexity: O(total_elements) where total_elements is sum of all matrix sizes
  /// @note Tolerance: Uses Eigen::MatrixXf::isApprox() default tolerance
  bool isApprox(const std::vector<Eigen::MatrixXf>& vec_mat1,
                     const std::vector<Eigen::MatrixXf>& vec_mat2);

  /// @brief Configures Eigen thread count for parallel operations
  /// @param[in] n_threads_in Requested number of threads (clamped to n_core if exceeded)
  /// @note Thread-safety: Not thread-safe; should be called during initialization only
  /// @note Complexity: O(1)
  /// @note Clamping: If n_threads_in > n_core (hardware_concurrency/2), uses n_core instead
  /// @note Side effects: Calls Eigen::setNbThreads() and Eigen::initParallel(); outputs to std::cout
  void set_threads_Eigen( const int n_threads_in );

  /// @brief Converts dense matrix to sparse matrix by pruning small elements
  /// @param[in] dense Input dense matrix
  /// @param[in] reference Reference magnitude for threshold calculation
  /// @param[in] epsilon Relative tolerance (default: machine epsilon for float)
  /// @return Sparse matrix containing only elements with |value| >= reference * epsilon
  /// @note Thread-safety: Thread-safe for concurrent calls with different matrices
  /// @note Complexity: O(nrow * ncol)
  /// @note Memory layout: Result is ColMajor sparse matrix
  /// @note Threshold: Elements with absolute value < reference * epsilon are considered zero
  inline SpMatf convertDenseToSparse(
    const Eigen::MatrixXf& dense, const float reference
    , const float epsilon = std::numeric_limits<float>::epsilon() ){
    return dense.sparseView(reference, epsilon);
  };

  /// @brief Converts sparse matrix to dense matrix
  /// @param[in] sparse Input sparse matrix
  /// @return Dense matrix with all elements (zeros filled in)
  /// @note Thread-safety: Thread-safe for concurrent calls with different matrices
  /// @note Complexity: O(nrow * ncol)
  /// @note Memory layout: Result is ColMajor dense matrix
  inline Eigen::MatrixXf convertSparseToDense(const SpMatf& sparse)
  {
    return Eigen::MatrixXf(sparse);
  };

  /// @brief Compares sparse matrix with dense matrix for approximate equality
  /// @param[in] mat_sps Sparse matrix
  /// @param[in] mat_den Dense matrix
  /// @return true if matrices are approximately equal, false otherwise
  /// @note Thread-safety: Thread-safe for read-only comparison
  /// @note Complexity: O(nnz) where nnz is number of non-zeros in sparse matrix
  /// @note Tolerance: Uses Eigen's default isApprox() tolerance
  inline bool isApprox(const SpMatf& mat_sps, const Eigen::MatrixXf& mat_den)
  {
    return mat_sps.isApprox(mat_den);
  };

  /// @brief Checks if Eigen matrix data pointer is properly aligned for SIMD operations
  /// @tparam Derived Eigen matrix type (auto-deduced)
  /// @param[in] matrix Matrix to check alignment
  /// @return true if data pointer is aligned to EIGEN_MAX_ALIGN_BYTES boundary, false otherwise
  /// @note Thread-safety: Thread-safe for read-only check
  /// @note Complexity: O(1)
  /// @note Alignment: Checks alignment to EIGEN_MAX_ALIGN_BYTES (typically 16 or 32 bytes for AVX)
  /// @note Performance: Aligned matrices enable vectorized SIMD operations
  template <typename Derived>
  bool isAlignedMatrix(const Eigen::MatrixBase<Derived>& matrix) {
    return (reinterpret_cast<std::uintptr_t>(matrix.derived().data()) % EIGEN_MAX_ALIGN_BYTES) == 0;
  };

  ///@} ------------------------------------------------------------------


  // ======================================================================
  /// @name Eigen_output_functions
  /// @details Functions related to Eigen output
  /// @{

  /// @brief Outputs Eigen::VectorXf to open file handle
  /// @param[in,out] fout Output file handle (must be open for writing)
  /// @param[in] vec Vector to output
  /// @note Thread-safety: Not thread-safe; external synchronization required for file I/O
  /// @note Complexity: O(vec.size())
  /// @note Format: Each line contains "irow value" in scientific notation
  void out_vecxf( FILE *fout, const Eigen::VectorXf &vec );

  /// @brief Outputs Eigen::VectorXf to ASCII file
  /// @param[in] pathout Output file path
  /// @param[in] vec Vector to output
  /// @note Thread-safety: Not thread-safe; external synchronization required for file I/O
  /// @note Complexity: O(vec.size())
  /// @note Format: Each line contains "irow value" in scientific notation
  void out_vecxf( const std::filesystem::path& pathout, const Eigen::VectorXf &vec );

  /// @brief Outputs element-wise difference of two vectors to file handle, skipping zero pairs
  /// @param[in,out] fout Output file handle (must be open for writing)
  /// @param[in] vec0 First vector
  /// @param[in] vec1 Second vector (must have same size as vec0)
  /// @throws std::runtime_error If vec0.rows() != vec1.rows()
  /// @note Thread-safety: Not thread-safe; external synchronization required for file I/O
  /// @note Complexity: O(vec0.size())
  /// @note Format: Each line contains "irow val0 val1 diff" where diff = val1 - val0
  /// @note Filtering: Skips rows where both vec0(i) == 0.0 and vec1(i) == 0.0
  void out_vecxf_diff( FILE *fout
  , const Eigen::VectorXf &vec0, const Eigen::VectorXf &vec1 );

  /// @brief Outputs element-wise difference of two vectors to ASCII file, skipping zero pairs
  /// @param[in] pathout Output file path
  /// @param[in] vec0 First vector
  /// @param[in] vec1 Second vector (must have same size as vec0)
  /// @throws std::runtime_error If vec0.rows() != vec1.rows()
  /// @note Thread-safety: Not thread-safe; external synchronization required for file I/O
  /// @note Complexity: O(vec0.size())
  /// @note Format: Each line contains "irow val0 val1 diff" where diff = val1 - val0
  /// @note Filtering: Skips rows where both vec0(i) == 0.0 and vec1(i) == 0.0
  void out_vecxf_diff( const std::filesystem::path& pathout
  , const Eigen::VectorXf &vec0, const Eigen::VectorXf &vec1 );

  /// @brief Outputs all elements of Eigen::MatrixXf to file handle in ColMajor order
  /// @param[in,out] fout Output file handle (must be open for writing)
  /// @param[in] mat Matrix to output
  /// @note Thread-safety: Not thread-safe; external synchronization required for file I/O
  /// @note Complexity: O(nrow * ncol)
  /// @note Format: Each line contains "irow icol value" in scientific notation
  /// @note Memory layout: Outputs in ColMajor order (iterates columns first)
  void out_matrix( FILE *fout, const Eigen::MatrixXf &mat );

  /// @brief Outputs all elements of Eigen::MatrixXf to ASCII file in ColMajor order
  /// @param[in] pathout Output file path
  /// @param[in] mat Matrix to output
  /// @note Thread-safety: Not thread-safe; external synchronization required for file I/O
  /// @note Complexity: O(nrow * ncol)
  /// @note Format: Each line contains "irow icol value" in scientific notation
  /// @note Memory layout: Outputs in ColMajor order (iterates columns first)
  void out_matrix( const std::filesystem::path& pathout, const Eigen::MatrixXf &mat );

  /// @brief Outputs all elements of Eigen::MatrixXd to file handle in ColMajor order
  /// @param[in,out] fout Output file handle (must be open for writing)
  /// @param[in] mat Matrix to output
  /// @note Thread-safety: Not thread-safe; external synchronization required for file I/O
  /// @note Complexity: O(nrow * ncol)
  /// @note Format: Each line contains "irow icol value" in scientific notation
  /// @note Memory layout: Outputs in ColMajor order (iterates columns first)
  void out_matrix( FILE *fout, const Eigen::MatrixXd &mat );

  /// @brief Outputs all elements of Eigen::MatrixXd to ASCII file in ColMajor order
  /// @param[in] pathout Output file path
  /// @param[in] mat Matrix to output
  /// @note Thread-safety: Not thread-safe; external synchronization required for file I/O
  /// @note Complexity: O(nrow * ncol)
  /// @note Format: Each line contains "irow icol value" in scientific notation
  /// @note Memory layout: Outputs in ColMajor order (iterates columns first)
  void out_matrix( const std::filesystem::path& pathout, const Eigen::MatrixXd &mat );

  /// @brief Outputs non-zero elements of Eigen::MatrixXf to file handle with dimension header
  /// @param[in,out] fout Output file handle (must be open for writing)
  /// @param[in] mat Matrix to output
  /// @param[in] epsilon Threshold for zero detection (elements with |value| < epsilon are skipped)
  /// @note Thread-safety: Not thread-safe; external synchronization required for file I/O
  /// @note Complexity: O(nrow * ncol)
  /// @note Format: First line "# nrow ncol", subsequent lines "irow icol value" for non-zeros only
  /// @note Memory layout: Outputs in ColMajor order (iterates columns first)
  void out_matrix_non_zero(
    FILE *fout, const Eigen::MatrixXf &mat, const float epsilon );

  /// @brief Outputs non-zero elements of Eigen::MatrixXf to ASCII file with dimension header
  /// @param[in] pathout Output file path
  /// @param[in] mat Matrix to output
  /// @param[in] epsilon Threshold for zero detection (default: std::numeric_limits<float>::min())
  /// @note Thread-safety: Not thread-safe; external synchronization required for file I/O
  /// @note Complexity: O(nrow * ncol)
  /// @note Format: First line "# nrow ncol", subsequent lines "irow icol value" for non-zeros only
  /// @note Memory layout: Outputs in ColMajor order (iterates columns first)
  void out_matrix_non_zero(
      const std::filesystem::path& pathout, const Eigen::MatrixXf &mat
    , const float epsilon=std::numeric_limits<float>::min() );

  /// @brief Outputs non-zero elements of Eigen::MatrixXd to file handle with dimension header
  /// @param[in,out] fout Output file handle (must be open for writing)
  /// @param[in] mat Matrix to output
  /// @param[in] epsilon Threshold for zero detection (elements with |value| < epsilon are skipped)
  /// @note Thread-safety: Not thread-safe; external synchronization required for file I/O
  /// @note Complexity: O(nrow * ncol)
  /// @note Format: First line "# nrow ncol", subsequent lines "irow icol value" for non-zeros only
  /// @note Memory layout: Outputs in ColMajor order (iterates columns first)
  void out_matrix_non_zero(
    FILE *fout, const Eigen::MatrixXd &mat, const double epsilon );

  /// @brief Outputs non-zero elements of Eigen::MatrixXf to ASCII file (convenience wrapper)
  /// @param[in] pathout Output file path
  /// @param[in] mat Matrix to output
  /// @note Thread-safety: Not thread-safe; external synchronization required for file I/O
  /// @note Complexity: O(nrow * ncol)
  /// @note Format: First line "# nrow ncol", subsequent lines "irow icol value" for non-zeros only
  /// @note Wrapper: Calls out_matrix_non_zero() internally
  void out_matxf(
    const std::filesystem::path& pathout, const Eigen::MatrixXf &mat);

  /// @brief Outputs non-zero elements of Eigen::MatrixXf to ASCII file in RowMajor order
  /// @param[in] pathout Output file path
  /// @param[in] mat Matrix to output
  /// @note Thread-safety: Not thread-safe; external synchronization required for file I/O
  /// @note Complexity: O(nrow * ncol)
  /// @note Format: First line "# nrow ncol output_row_major", subsequent lines "irow icol value"
  /// @note Memory layout: Outputs in RowMajor order (iterates rows first) unlike other output functions
  /// @note Filtering: Skips elements with |value| < std::numeric_limits<float>::min()
  void out_matxf_rowmajor(
    const std::filesystem::path& pathout, const Eigen::MatrixXf &mat);

  /// @brief Outputs non-zero elements of Eigen::MatrixXd to ASCII file with dimension header
  /// @param[in] pathout Output file path
  /// @param[in] mat Matrix to output
  /// @param[in] epsilon Threshold for zero detection (default: std::numeric_limits<double>::epsilon())
  /// @note Thread-safety: Not thread-safe; external synchronization required for file I/O
  /// @note Complexity: O(nrow * ncol)
  /// @note Format: First line "# nrow ncol", subsequent lines "irow icol value" for non-zeros only
  /// @note Memory layout: Outputs in ColMajor order (iterates columns first)
  void out_matrix_non_zero(
    const std::filesystem::path& pathout, const Eigen::MatrixXd &mat
  , const double epsilon=std::numeric_limits<double>::epsilon() );

  /// @brief Outputs non-zero elements of sparse matrix to file handle with dimension header
  /// @param[in,out] fout Output file handle (must be open for writing)
  /// @param[in] mat Sparse matrix to output
  /// @param[in] epsilon Threshold for zero detection (elements with |value| < epsilon are skipped)
  /// @note Thread-safety: Not thread-safe; external synchronization required for file I/O
  /// @note Complexity: O(nnz) where nnz is number of non-zeros
  /// @note Format: First line "# nrow ncol", subsequent lines "irow icol value" for non-zeros only
  /// @note Memory layout: Iterates in ColMajor order using SpMatf::InnerIterator
  void out_matrix_non_zero(
    FILE *fout, const SpMatf &mat, const float epsilon );
  
  /// @brief Outputs non-zero elements of sparse matrix to ASCII file with dimension header
  /// @param[in] pathout Output file path
  /// @param[in] mat Sparse matrix to output
  /// @param[in] epsilon Threshold for zero detection (default: std::numeric_limits<float>::epsilon())
  /// @note Thread-safety: Not thread-safe; external synchronization required for file I/O
  /// @note Complexity: O(nnz) where nnz is number of non-zeros
  /// @note Format: First line "# nrow ncol", subsequent lines "irow icol value" for non-zeros only
  /// @note Memory layout: Iterates in ColMajor order using SpMatf::InnerIterator
  void out_matrix_non_zero(
    const std::filesystem::path& pathout, const SpMatf &mat
  , const float epsilon=std::numeric_limits<float>::epsilon() );

  /// @brief Outputs each matrix in vector to separate files with indexed filenames
  /// @tparam MatrixType Matrix type (must be Eigen::MatrixXf, Eigen::MatrixXd, or SpMatf)
  /// @param[in] suffix Suffix for output filenames (generates "out_mat_detXX_suffix.tmp")
  /// @param[in] vec_mat Vector of matrices to output
  /// @param[in] epsilon Threshold for zero detection (default: type-specific machine epsilon)
  /// @note Thread-safety: Not thread-safe; external synchronization required for file I/O
  /// @note Complexity: O(total_nnz) where total_nnz is sum of non-zeros in all matrices
  /// @note Format: Each file named "out_mat_detXX_suffix.tmp" where XX is zero-padded index
  /// @note Type safety: static_assert ensures MatrixType is one of the allowed types
  template<typename MatrixType>
  void out_vec_mat_non_zero(
      const std::string& suffix,
      const std::vector<MatrixType>& vec_mat,
      typename MatrixType::Scalar epsilon = std::numeric_limits<typename MatrixType::Scalar>::epsilon())
  {
    // Verify at compile time that template type is an allowed type
    static_assert(
      std::is_same_v<MatrixType, Eigen::MatrixXf> ||
      std::is_same_v<MatrixType, Eigen::MatrixXd> ||
      std::is_same_v<MatrixType, SpMatf>,
      "MatrixType must be either Eigen::MatrixXf, Eigen::MatrixXd, or SpMatf"
    );
    
    char cfname[512];

    for (size_t i = 0; i < vec_mat.size(); ++i) {
      // Generate output filename
      sprintf(cfname, "out_mat_det%02zu_%s.tmp", i, suffix.c_str());
      std::filesystem::path pathout(cfname);

      // Call existing function to output non-zero elements of matrix
      out_matrix_non_zero(pathout, vec_mat[i], epsilon);
    }
  };
  ///@} ------------------------------------------------------------------


};

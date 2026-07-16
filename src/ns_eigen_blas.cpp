// src/ns_eigen_blas.cpp
#include "ns_eigen_blas.hpp"
#include "ns_mylogger.hpp"
#include "ns_type_definitions.hpp"
#include "ns_myapp.hpp"

// Inverse matrix function with float precision
Eigen::MatrixXf eigen_blas::calc_inverse_float( const Eigen::MatrixXf &matxf )
{
  // Specify the number of threads.
  // If 0 is specified, BLAS automatically selects the optimal number of threads.
  const int n_threads = omp_get_max_threads();

  // Specify the number of OpenMP threads
  // omp_set_num_threads(n_threads);

  // Specify the number of Eigen threads
  // eigen_blas::set_threads_Eigen(n_threads);

  // Get number of rows and columns
  const int nrow = matxf.rows();
  const int ncol = matxf.cols();

  // Check if matrix is square
  if(nrow!=ncol) THROW_ERROR("nrow!=ncol");

  // Allocate memory for inverse matrix array
  float* mat_inv = new float[nrow * ncol];

  // Array for LU decomposition
  std::vector<blas_int> ipiv(nrow); // Use blas_int as an alias for long long type

  // Convert Eigen::MatrixXf to float*
  Eigen::Map<Eigen::MatrixXf>(mat_inv, nrow, ncol) = matxf;
  bool tf_same1 = eigen_blas::is_same_array_matrixXf(mat_inv,matxf);
  if(tf_same1==false) THROW_ERROR("mat_inv is not equal to matxf");

  // LU decomposition (float)
  const blas_int lu_result = LAPACKE_sgetrf(LAPACK_COL_MAJOR, nrow, ncol, mat_inv, ncol, ipiv.data());
  if( lu_result > 0 ) THROW_ERROR2("eigen_blas::calc_inverse_float: The k-th diagonal element is zero. Matrix is singular and inverse does not exist. k={}",lu_result);
  if( lu_result < 0 ) THROW_ERROR("eigen_blas::calc_inverse_float: Function was not called correctly. One or more arguments are invalid.");

  // Inverse matrix calculation (float)
  const blas_int inv_result = LAPACKE_sgetri(LAPACK_COL_MAJOR, ncol, mat_inv, ncol, ipiv.data());
  if( inv_result > 0 ) THROW_ERROR2("eigen_blas::calc_inverse_float: The k-th diagonal element is zero. Matrix is singular and inverse does not exist. k={}",inv_result);
  if( inv_result < 0 ) THROW_ERROR("eigen_blas::calc_inverse_float: Function was not called correctly. One or more arguments are invalid.");

  // Convert float* to Eigen::MatrixXf
  // Eigen::Map<Eigen::MatrixXf> matxf_inv(mat_inv, nrow, ncol); // ! mistake. it's not deep copy.
  Eigen::MatrixXf matxf_inv = Eigen::Map<Eigen::MatrixXf>(mat_inv, nrow, ncol);

  bool tf_same2 = eigen_blas::is_same_array_matrixXf(mat_inv,matxf_inv);
  if(tf_same2==false) THROW_ERROR("eigen_blas::calc_inverse_float: mat_inv is not equal to matxf_inv");

  // Free memory
  delete[] mat_inv;

  // Return inverse matrix
  return matxf_inv;
}

// Returns the inverse matrix of Eigen::MatrixXf
Eigen::MatrixXf eigen_blas::getInverseMatrixFloat(const Eigen::MatrixXf& inputMatrix)
{

  // Alignment check for MatrixXf
  const bool is_aligned_input = eigen_blas::isAlignedMatrix(inputMatrix);
  LOG_DEBUG("is_aligned_input = {}", is_aligned_input ? "true" : "false");
  if( !is_aligned_input ) {
    LOG_ERROR("inputMatrix.data() is not aligned");
    THROW_ERROR("inputMatrix.data() is not aligned");
  }

  // Create a copy of input matrix (used for inverse matrix calculation)
  Eigen::MatrixXf matrix = inputMatrix.eval(); // copy with safety

  // Get matrix size and cast to blas_int type
  const blas_int rows = to_blas_int(matrix.rows(), "eigen_blas::getInverseMatrixFloat");
  const blas_int cols = to_blas_int(matrix.cols(), "eigen_blas::getInverseMatrixFloat");
  
  // check matrix.data() pointer
  if( matrix.data()==nullptr) THROW_ERROR("matrix.data()==nullptr");
  
  // check rows and cols
  if( rows <= 0 || cols <= 0 ) THROW_ERROR3("rows <= 0 || cols <= 0",rows,cols);

  // Pivot array for LU decomposition (declared as blas_int type)
  std::vector<blas_int> pivot(std::min(rows, cols));

  // Execute LU decomposition
  // blas_int info = LAPACKE_sgetrf(LAPACK_COL_MAJOR, rows, cols, matrix.data(), cols, pivot.data());
  blas_int info = LAPACKE_sgetrf(LAPACK_COL_MAJOR, rows, cols, matrix.data(), rows, pivot.data());

  if (info != 0) {
    THROW_ERROR2("eigen_blas::getInverseMatrixFloat: LU decomposition failed with info", info);
  }

  // Calculate inverse matrix
  // info = LAPACKE_sgetri(LAPACK_COL_MAJOR, cols, matrix.data(), cols, pivot.data());
  info = LAPACKE_sgetri(LAPACK_COL_MAJOR, cols, matrix.data(), rows, pivot.data());

  if (info != 0) {
    THROW_ERROR2("eigen_blas::getInverseMatrixFloat: Inverse matrix calculation failed with info", info);
  }

  // alignment check of output matrix
  const bool is_aligned_output = eigen_blas::isAlignedMatrix(matrix);
  LOG_DEBUG("is_aligned_output = {}", is_aligned_output ? "true" : "false");
  if( !is_aligned_output ) {
    LOG_ERROR("matrix.data() is not aligned");
    THROW_ERROR("matrix.data() is not aligned");
  }

  // ! return matrix;
  return Eigen::MatrixXf(matrix);
}

// Inverse matrix function with double precision
Eigen::MatrixXf eigen_blas::getInverseMatrixDouble( const Eigen::MatrixXf &matxf)
{
  const int n_threads = omp_get_max_threads();

  const int nrow = matxf.rows();
  const int ncol = matxf.cols();

  // Check if matrix is square
  if(nrow!=ncol) THROW_ERROR("nrow!=ncol");

  // Array for inverse matrix
  double* mat_inv = new double[nrow * ncol];

  // Array for LU decomposition
  std::vector<blas_int> ipiv(nrow); // Use blas_int as an alias for long long type

  // Convert Eigen::MatrixXf to double*
  Eigen::MatrixXd matxd = matxf.cast<double>();
  Eigen::Map<Eigen::MatrixXd>(mat_inv, nrow, ncol) = matxd;

  const blas_int lu_result = LAPACKE_dgetrf(LAPACK_COL_MAJOR, nrow, ncol, mat_inv, ncol, ipiv.data());
  if( lu_result > 0 ) THROW_ERROR2("eigen_blas::getInverseMatrixDouble: The k-th diagonal element is zero. Matrix is singular and inverse does not exist. k={}",lu_result);
  if( lu_result < 0 ) THROW_ERROR("eigen_blas::getInverseMatrixDouble: Function was not called correctly. One or more arguments are invalid.");

  // Inverse matrix calculation (double)
  const blas_int inv_result = LAPACKE_dgetri(LAPACK_COL_MAJOR, ncol, mat_inv, ncol, ipiv.data());
  if( inv_result > 0 ) THROW_ERROR2("eigen_blas::getInverseMatrixDouble: The k-th diagonal element is zero. Matrix is singular and inverse does not exist. k={}",inv_result);
  if( inv_result < 0 ) THROW_ERROR("eigen_blas::getInverseMatrixDouble: Function was not called correctly. One or more arguments are invalid.");

  // Convert double* to float*
  float* mat_inv_f = new float[nrow * ncol];
  eigen_blas::convertDoubleToFloat(mat_inv, mat_inv_f, nrow * ncol);

  // Convert float* to Eigen::MatrixXf
  // Eigen::Map<Eigen::MatrixXf> matxf_inv(mat_inv_f, nrow, ncol); // ! mistake. it's not deep copy.
  Eigen::MatrixXf matxf_inv = Eigen::Map<Eigen::MatrixXf>(mat_inv_f, nrow, ncol);

  // Free memory
  delete[] mat_inv;
  delete[] mat_inv_f;

  // Return inverse matrix
  return matxf_inv;
}

// Inverse matrix function with double precision
Eigen::MatrixXd eigen_blas::getInverseMatrixDouble( const Eigen::MatrixXd &matxd )
{
  // Specify the number of threads.
  // If 0 is specified, BLAS automatically selects the optimal number of threads.
  const int n_threads = omp_get_max_threads();

  const int nrow = matxd.rows();
  const int ncol = matxd.cols();

  // Check if matrix is square
  if(nrow!=ncol) THROW_ERROR("nrow!=ncol");

  // Array for inverse matrix
  double* mat_inv = new double[nrow * ncol];

  // Array for LU decomposition
  std::vector<blas_int> ipiv(nrow); // Use blas_int as an alias for long long type

  // Convert Eigen::MatrixXd to double*
  Eigen::Map<Eigen::MatrixXd>(mat_inv, nrow, ncol) = matxd;

  const blas_int lu_result = LAPACKE_dgetrf(LAPACK_COL_MAJOR, nrow, ncol, mat_inv, ncol, ipiv.data());
  if( lu_result > 0 ) THROW_ERROR2("eigen_blas::getInverseMatrixDouble: The k-th diagonal element is zero. Matrix is singular and inverse does not exist. k={}",lu_result);
  if( lu_result < 0 ) THROW_ERROR("eigen_blas::getInverseMatrixDouble: Function was not called correctly. One or more arguments are invalid.");

  // Inverse matrix calculation (double)
  const blas_int inv_result = LAPACKE_dgetri(LAPACK_COL_MAJOR, ncol, mat_inv, ncol, ipiv.data());
  if( inv_result > 0 ) THROW_ERROR2("eigen_blas::getInverseMatrixDouble: The k-th diagonal element is zero. Matrix is singular and inverse does not exist. k={}",inv_result);
  if( inv_result < 0 ) THROW_ERROR("eigen_blas::getInverseMatrixDouble: Function was not called correctly. One or more arguments are invalid.");

  // Convert double* to Eigen::MatrixXd
  // Eigen::Map<Eigen::MatrixXd> matxd_inv(mat_inv, nrow, ncol); // ! mistake. it's not deep copy.
  Eigen::MatrixXd matxd_inv = Eigen::Map<Eigen::MatrixXd>(mat_inv, nrow, ncol);

  // Free memory
  delete[] mat_inv;

  // Return inverse matrix
  return matxd_inv;
}

// @brief Multiplies two Eigen::MatrixXf matrices using BLAS with float precision
// @note  mat_A, mat_B, and mat_C must all be ColMajor
Eigen::MatrixXf eigen_blas::multiplyMatrixColMajor(
  const Eigen::MatrixXf& mat_A, const Eigen::MatrixXf& mat_B)
{
  LOG_INFO("");
  LOG_INFO("Notice that mat_A, mat_B, and return mat_C should be colmajor");

  eigen_blas::check_colmajor_matxf(mat_A);
  eigen_blas::check_colmajor_matxf(mat_B);

  // Get matrix sizes
  const blas_int rows_A = to_blas_int(mat_A.rows(), "eigen_blas::multiplyMatrixColMajor(float)");
  const blas_int cols_A = to_blas_int(mat_A.cols(), "eigen_blas::multiplyMatrixColMajor(float)");
  const blas_int rows_B = to_blas_int(mat_B.rows(), "eigen_blas::multiplyMatrixColMajor(float)");
  const blas_int cols_B = to_blas_int(mat_B.cols(), "eigen_blas::multiplyMatrixColMajor(float)");

  // Check matrix sizes
  if (cols_A != rows_B) {
    THROW_ERROR3("eigen_blas::multiplyMatrixColMajor - cols_A != rows_B, can't multiply", cols_A, rows_B);
  }

  const blas_int rows_C = rows_A;
  const blas_int cols_C = cols_B;

  // Create matrix for storing result
  Eigen::MatrixXf mat_C(rows_C, cols_C);
  eigen_blas::check_colmajor_matxf(mat_C);

  // Set constexpr constant values
  constexpr float alpha = 1.0f;
  constexpr float beta = 0.0f;

  // Matrix multiplication
  cblas_sgemm(
    CblasColMajor, CblasNoTrans, CblasNoTrans,
    rows_A, // Number of rows in A
    cols_B, // Number of columns in B
    cols_A, // Number of columns in A
    alpha, // Scalar alpha
    mat_A.data(), // A
    rows_A, // lda (rows_A for column-major)
    mat_B.data(), // B
    rows_B, // ldb (rows_B for column-major)
    beta, // Scalar beta
    mat_C.data(), // C
    rows_A // ldc (rows_C for column-major)
  );

  return mat_C;
}

// @brief Multiplies two Eigen::MatrixXd matrices using BLAS with double precision
//        mat_A, mat_B, and mat_C must all be ColMajor
Eigen::MatrixXd eigen_blas::multiplyMatrixColMajor(
  const Eigen::MatrixXd& mat_A, const Eigen::MatrixXd& mat_B)
{
  LOG_INFO("");
  LOG_INFO("Notice that mat_A, mat_B, and return mat_C should be colmajor");

  eigen_blas::check_colmajor_matxd(mat_A);
  eigen_blas::check_colmajor_matxd(mat_B);

  // Get matrix sizes
  const blas_int rows_A = to_blas_int(mat_A.rows(), "eigen_blas::multiplyMatrixColMajor(double)");
  const blas_int cols_A = to_blas_int(mat_A.cols(), "eigen_blas::multiplyMatrixColMajor(double)");
  const blas_int rows_B = to_blas_int(mat_B.rows(), "eigen_blas::multiplyMatrixColMajor(double)");
  const blas_int cols_B = to_blas_int(mat_B.cols(), "eigen_blas::multiplyMatrixColMajor(double)");

  // Check matrix sizes
  if (cols_A != rows_B) {
    THROW_ERROR3("eigen_blas::multiplyMatrixColMajor - cols_A != rows_B, can't multiply", cols_A, rows_B);
  }

  const blas_int rows_C = rows_A;
  const blas_int cols_C = cols_B;

  // Create matrix for storing result
  Eigen::MatrixXd mat_C(rows_C, cols_C);
  eigen_blas::check_colmajor_matxd(mat_C);

  // Set constexpr constant values
  constexpr double alpha = 1.0;
  constexpr double beta = 0.0;

  // Matrix multiplication
  cblas_dgemm(
    CblasColMajor, CblasNoTrans, CblasNoTrans,
    rows_A, // Number of rows in A
    cols_B, // Number of columns in B
    cols_A, // Number of columns in A
    alpha, // Scalar alpha
    mat_A.data(), // A
    rows_A, // lda (rows_A for column-major)
    mat_B.data(), // B
    rows_B, // ldb (rows_B for column-major)
    beta, // Scalar beta
    mat_C.data(), // C
    rows_A // ldc (rows_C for column-major)
  );

  return mat_C;
}

// Convert double* to float*
void eigen_blas::convertDoubleToFloat(double* input, float* output, size_t size)
{
  for (size_t i = 0; i < size; i++) {
    output[i] = static_cast<float>(input[i]);
  }
}

// out c_array float
// row
void eigen_blas::out_matrix_array( const float *mat
, const int nrow, const int ncol, const std::filesystem::path& pathout )
{
  LOG_INFO("to pathout={}", pathout.string());
  FILE *fout = myapp::get_fout(pathout);
  // colmajor output
  for(int icol = 0; icol < ncol; icol++){
    for(int irow = 0; irow < nrow; irow++){
      const int colmajor_index = icol * nrow + irow;
      fprintf(fout,"%d %d %E\n",irow,icol,mat[colmajor_index]);
    }
  }
  myapp::close(fout,pathout);
}

// out c_array double
void eigen_blas::out_matrix_array( const double *mat
, const int nrow, const int ncol, const std::filesystem::path& pathout )
{
  LOG_INFO("to pathout={}", pathout.string());
  FILE *fout = myapp::get_fout(pathout);
  // colmajor output
  for(int icol = 0; icol < ncol; icol++){
    for(int irow = 0; irow < nrow; irow++){
      const int colmajor_index = icol * nrow + irow;
      fprintf(fout,"%d %d %E\n",irow,icol,mat[colmajor_index]);
    }
  }
  myapp::close(fout,pathout);
}

/// @brief output vec_mat to ascii file
/// @note zero value is not output
void eigen_blas::out_vec_mat(
  const std::filesystem::path& pathout, const std::vector<Eigen::MatrixXf> &vec_mat )
{
  const int n_det = vec_mat.size();
  FILE *fout = myapp::get_fout(pathout);
  for(int i=0;i<n_det;++i){
    const Eigen::MatrixXf &mat = vec_mat.at(i);
    const int nrow = mat.rows();
    const int ncol = mat.cols();
    fprintf(fout,"det_id = %d %d %d\n",i,nrow,ncol);
    for(int icol = 0; icol < ncol; icol++){
      for(int irow = 0; irow < nrow; irow++){
        const int colmajor_index = icol * nrow + irow;
        float value = mat(colmajor_index);
        if( value == 0.0 ) continue; // zero value is not output
        fprintf(fout,"%d %d %E\n",irow,icol,mat(colmajor_index));
      }
    }
  }
  myapp::close(fout,pathout);
}

// Compare array and EigenMatrix (float)
bool eigen_blas::is_same_array_matrixXf( const float* mat, const Eigen::MatrixXf &matxf )
{
  const int nrow = matxf.rows();
  const int ncol = matxf.cols();
  for(int icol=0;icol<ncol;icol++){
    for(int irow=0;irow<nrow;irow++){
      const int colmajor_index = icol * nrow + irow;
      const float value1 = mat[colmajor_index];
      const float value2 = matxf(irow,icol);
      if( std::abs(value1-value2) > 10.0*std::numeric_limits<float>::epsilon() ){
        return false;
      }
    }
  }
  return true;
}

// Compare array and EigenMatrix (double)
bool eigen_blas::is_same_array_matrixXd( const double* mat, const Eigen::MatrixXd &matxd )
{
  const int nrow = matxd.rows();
  const int ncol = matxd.cols();
  for(int icol=0;icol<ncol;icol++){
    for(int irow=0;irow<nrow;irow++){
      const int colmajor_index = icol * nrow + irow;
      const double value1 = mat[colmajor_index];
      const double value2 = matxd(irow,icol);
      if( std::abs(value1-value2) > 10.0*std::numeric_limits<double>::epsilon() ){
        return false;
      }
    }
  }
  return true;
}


// check matrixxf
bool eigen_blas::check_MatrixXf( const Eigen::MatrixXf &mat ){
  const int rows = mat.rows();
  const int cols = mat.cols();
  if(rows<=0){
    LOG_ERROR("rows={}", rows);
    return false;
  }
  if(cols<=0){
    LOG_ERROR("cols={}", cols);
    return false;
  }
  if( mat.data()==nullptr ){
    LOG_ERROR("mat.data()==nullptr");
    return false;
  }
  return true;
}

// Convert std::vector<double> to Eigen::VectorXf
Eigen::VectorXf eigen_blas::get_vecxf(const std::vector<double>& vec)
{
  Eigen::VectorXf eigen_density(vec.size());
  for (size_t i = 0; i < vec.size(); ++i) {
    eigen_density(i) = static_cast<float>(vec[i]);
  }
  return eigen_density;
}

// @brief get Eigen::VectorXf from Eigen::VectorXf divided Eigen::VectorXf
Eigen::VectorXf eigen_blas::get_vecxf_div( const Eigen::VectorXf &vec1, const Eigen::VectorXf &vec2 )
{
  const int nrow1 = vec1.rows();
  const int nrow2 = vec2.rows();
  if( nrow1 != nrow2 ) THROW_ERROR("vec1.rows() != vec2.rows()");
  Eigen::VectorXf vec3 = Eigen::VectorXf::Zero(nrow1);
  for(int irow=0;irow<nrow1;irow++){
    const float value1 = vec1(irow);
    const float value2 = vec2(irow);
    if (value2 == 0.0) {
      vec3(irow) = -931931.0;
    } else {
      vec3(irow) = value1 / value2;
    }
  }
  return vec3;
}


// Read Eigen::MatrixXf non-zero elements
Eigen::MatrixXf eigen_blas::read_matrix_non_zero( const std::filesystem::path &path_in )
{
  fprintf(stderr,"eigen_blas::read_matrix_non_zero reading file %s ...\n",path_in.c_str());
  std::ifstream reading_file;
  if(!std::filesystem::exists(path_in)) THROW_ERROR( path_in.string() + " does not exist");
  reading_file.open(path_in,std::ios::in);
  std::string reading_line_str;
  std::vector<std::string> vec_str;
  
  // Read 1st line: nrows, ncols
  std::getline(reading_file, reading_line_str);
  constexpr char delimiter = myapp::char_delim_default;
  constexpr char comment = 'a';
  // Cannot use constexpr char comment = myapp::char_commentout_default;
  // Because the first character of the first line of path_in has '#'.
  // It's a bit tricky, but we read rows and columns by explicitly not treating '#' as a comment on the first line.
  vec_str = myapp::split(reading_line_str,delimiter,comment);
  if( vec_str.size()!=3 ) THROW_ERROR2("eigen_blas::read_matrix_non_zero vec_str.size()!=3",vec_str.size());

  // Check that first character is '#'
  if( vec_str.at(0).c_str()[0]!=myapp::char_commentout_default ) THROW_ERROR("eigen_blas::read_matrix_non_zero: vec_str.at(0).c_str()[0]!=myapp::char_commentout_default");
  const int nrow = std::stoi( vec_str.at(1) );
  const int ncol = std::stoi( vec_str.at(2) );

  // memory allocation of matrix
  Eigen::MatrixXf mat = Eigen::MatrixXf::Zero(nrow,ncol);
  fprintf(stderr,"********************\n");
  fprintf(stderr,"mat.rows=%d, cols=%d\n",nrow,ncol);
  fprintf(stderr,"********************\n");

  // Start loop to read data from line 2 onwards
  while( std::getline(reading_file, reading_line_str) ){
    // Read data
    vec_str = myapp::split(reading_line_str);
    if( vec_str.size()!=3 ) THROW_ERROR("eigen_blas::read_matrix_non_zero vec_str.size()!=3");
    const int irow = std::stoi( vec_str.at(0) );
    const int icol = std::stoi( vec_str.at(1) );
    const float value = std::stof(vec_str.at(2));
    mat(irow,icol) = value;
  }

  fprintf(stderr,"eigen_blas::read_matrix_non_zero done.\n");
  return mat;
}



// sum up cols and make vector
std::vector<double> eigen_blas::mat_to_vec( const Eigen::MatrixXf &mat )
{
  const size_t nrow = mat.rows();
  const size_t ncol = mat.cols();
  std::vector<double> ret_vec;
  // ret_vec.reserve(nrow); // error " ret_vec.at(irow) = sum_cols; "
  ret_vec.resize(nrow);

  for(int irow=0; irow<nrow; irow++ ){
    // for(int icol=0; icol<ncol; icol++ ) sum_cols += (double)( mat(irow,icol) );
    ret_vec.at(irow) = (double) mat.row(irow).sum();
  }
  return ret_vec;
}

// @brief Computes matrix-vector multiplication (inner product) in the most primitive way
// @param[in] matxf Matrix to multiply from the left
// @param[in] vecxf Vector to multiply from the right
// @return matxf*vecxf
// @details Outputs the multiplication result of each element obtained during computation to FILE*
Eigen::VectorXf eigen_blas::check_mat_by_vec( FILE *fout
  , const Eigen::MatrixXf &matxf, const Eigen::VectorXf &vecxf )
{
  const int nrow = matxf.rows();
  Eigen::VectorXf vec_ret(nrow);
  const int ncol = matxf.cols();

  LOG_INFO("nrow={} ncol={} vecxf.rows()={}",nrow,ncol,vecxf.rows());

  if(ncol!=vecxf.rows()) THROW_ERROR("eigen_blas::check_mat_by_vec: ncol!=vec.rows()");
  for(int irow=0;irow<nrow;irow++){
    const Eigen::VectorXf vecxf_matxf_irow = matxf.row(irow);
    float sum = 0.0;
    for(int icol=0;icol<ncol;icol++){
      const float factor_mat = vecxf_matxf_irow(icol);
      const float factor_vec = vecxf(icol);
      const float product = factor_mat * factor_vec;
      fprintf(fout,"%d %d %E %E %E\n",irow,icol,factor_mat,factor_vec,product);
      sum += product;
    }
    vec_ret(irow) = sum;
  }
  if(vec_ret.allFinite()==false) THROW_ERROR("vec_ret.allFinite()==false");
  return vec_ret;
}

Eigen::VectorXf eigen_blas::check_mat_by_vec( const std::filesystem::path& pathout
  , const Eigen::MatrixXf &matxf, const Eigen::VectorXf &vecxf )
{
  FILE *fout = myapp::get_fout(pathout);
  Eigen::VectorXf vec_ret = check_mat_by_vec(fout,matxf,vecxf);
  myapp::close(fout,pathout);
  return vec_ret;
}

/// @brief Computes the logarithm of Eigen::VectorXf
/// @param[in] vec Vector to compute logarithm
/// @return Vector with logarithm computed
/// @note Throws error if any element of vec is less than or equal to 0
Eigen::VectorXf eigen_blas::mp_loge_vecxf( const Eigen::VectorXf &vec )
{
  const int nrow = vec.rows();
  Eigen::VectorXf vec_log(nrow);
  #pragma omp parallel for
  for(int irow=0;irow<nrow;irow++){
    const float value = vec(irow);
    if( value <= 0.0 ) THROW_ERROR2("value <= 0.0, value=%E",value);
    vec_log(irow) = std::log(value);
  }
  return vec_log;
}

bool eigen_blas::check_colmajor_vecxf(const Eigen::VectorXf& vec)
{
  return vec.innerStride() == 1;
}

bool eigen_blas::check_colmajor_vecxd(const Eigen::VectorXd& vec)
{
  return vec.innerStride() == 1;
}

/// @brief Checks if Eigen::MatrixXf is ColMajor
void eigen_blas::check_colmajor_matxf(const Eigen::MatrixXf& mat)
{
  const int nrow = mat.rows();
  const int outer = mat.outerStride();
  const bool is_colmajor = (outer == nrow);
  if(  !is_colmajor ) {
    LOG_ERROR("is not ColMajor, outer={}, nrow={}", outer, nrow);
    THROW_ERROR("Eigen::MatrixXf is not ColMajor");
  }
}

/// @brief Checks if Eigen::MatrixXd is ColMajor
void eigen_blas::check_colmajor_matxd(const Eigen::MatrixXd& mat)
{
  const int nrow = mat.rows();
  const int outer = mat.outerStride();
  const bool is_colmajor = (outer == nrow);
  if(  !is_colmajor ) {
    LOG_ERROR("is not ColMajor, outer={}, nrow={}", outer, nrow);
    THROW_ERROR("Eigen::MatrixXf is not ColMajor");
  }
}

// Returns Eigen::MatrixXf created by concatenating two matrices (mat1, mat2) along rows
Eigen::MatrixXf eigen_blas::get_merged_matrix_rows(
  const Eigen::MatrixXf& mat1, const Eigen::MatrixXf& mat2)
{
  // Get matrix sizes
  const int rows1 = mat1.rows();
  const int cols1 = mat1.cols();
  const int rows2 = mat2.rows();
  const int cols2 = mat2.cols();

  // check mat1 and mat2
  if( cols1 != cols2 ) THROW_ERROR("mat1.cols() != mat2.cols()");
  if( rows1 <=0 ) THROW_ERROR("mat1.rows() <= 0");
  if( rows2 <=0 ) THROW_ERROR("mat2.rows() <= 0");

  // Create merged matrix. Expand in row direction
  Eigen::MatrixXf mat(rows1 + rows2, cols1);

  // Merge mat1 and mat2
  mat << mat1, mat2;

  return mat;
}

// Expands the rows of mat_large and adds mat_tobe_added
// If mat_large.rows()==0 or mat_large.cols()==0, simply copies mat_tobe_added
void eigen_blas::append_matrix_rows(Eigen::MatrixXf& mat_large, const Eigen::MatrixXf& mat_tobe_added)
{
  // Get matrix sizes
  const int rows_org = mat_large.rows();
  const int cols_org = mat_large.cols();
  const int rows_add = mat_tobe_added.rows();
  const int cols_add = mat_tobe_added.cols();

  // Size check
  if (cols_org != cols_add && rows_org != 0 && cols_org != 0) THROW_ERROR("mat_large.cols() != mat_tobe_added.cols()");
  if (rows_add <= 0) THROW_ERROR("mat_tobe_added.rows() <= 0");

  // If any dimension of mat_large is 0, simply copy mat_tobe_added
  if (rows_org == 0 || cols_org == 0) {
    mat_large = mat_tobe_added;
  } else {
    // Temporarily save existing matrix
    const Eigen::MatrixXf temp = mat_large;

    // Expand in row direction and merge
    mat_large.resize(rows_org + rows_add, cols_org);
    mat_large << temp, mat_tobe_added;
  }
  LOG_INFO("mat_large.rows()={}, mat_large.cols()={}", mat_large.rows(), mat_large.cols());
}

// Displays the number of rows and columns of each matrix in std::vector<Eigen::MatrixXf>
void eigen_blas::disp_vec_matxf(
    const std::vector<Eigen::MatrixXf>& vec_mat
  , const spdlog::level::level_enum& level)
{
  for(int i=0;i<vec_mat.size();++i){
    mylogger::g_logger->log(level
      ,"vec_mat[{}].rows()={}, vec_mat[{}].cols()={}"
      ,i,vec_mat.at(i).rows(),i,vec_mat.at(i).cols());
  }
}

// Displays the number of rows and columns of each matrix in std::vector<SpMatf>
void eigen_blas::disp_vec_spmatf(
  const std::vector<SpMatf>& vec_mat
, const spdlog::level::level_enum& level)
{
  for(int i=0;i<vec_mat.size();++i){
    mylogger::g_logger->log(level
      ,"vec_mat[{}].rows()={}, vec_mat[{}].cols()={}"
      ,i,vec_mat.at(i).rows(),i,vec_mat.at(i).cols());
  }
}

/// @brief Checks if two std::vector<Eigen::MatrixXf> are equal
bool eigen_blas::isApprox(const std::vector<Eigen::MatrixXf>& vec_mat1,
                          const std::vector<Eigen::MatrixXf>& vec_mat2)
{
  // If sizes differ, vectors do not match
  if (vec_mat1.size() != vec_mat2.size()) {
    return false;
  }

  // Check if each element of vec_mat1 matches vec_mat2
  for (int i = 0; i < vec_mat1.size(); ++i) {
    if ( !vec_mat1.at(i).isApprox(vec_mat2.at(i)) ) return false;
  }

  // If all above conditions are passed, consider equal
  return true;
}

// set Eigen::nbThreads() and Eigen::initParallel()
// if n_threads_in > n_core, n_core is used.
void eigen_blas::set_threads_Eigen( const int n_threads_in )
{
  const int n_core = std::thread::hardware_concurrency() / 2;
  int number_of_threads = n_threads_in;
  if( n_threads_in > n_core ) number_of_threads = n_core;
  Eigen::setNbThreads(  number_of_threads  );
  std::cout << "Eigen::nbThreads() = " << Eigen::nbThreads() << std::endl;
  Eigen::initParallel();
}

// output Eigen::VectorXf to FILE *
void eigen_blas::out_vecxf( FILE *fout, const Eigen::VectorXf &vec ){
  const size_t nrow = vec.rows();
  // colmajor output
  for(int irow=0; irow<nrow; irow++ ){
    fprintf(fout,"%d %7.4E\n",irow,vec(irow));
  }
}

// output Eigen::VectorXf to ascii file
void eigen_blas::out_vecxf(
  const std::filesystem::path& pathout, const Eigen::VectorXf &vec )
{
  FILE *fout = myapp::get_fout(pathout);
  out_vecxf(fout,vec);
  myapp::close(fout,pathout);
}

// output 2 vecxf to FILE *
void eigen_blas::out_vecxf_diff( FILE *fout
, const Eigen::VectorXf &vec0, const Eigen::VectorXf &vec1 )
{
  const int nrow = vec0.rows();
  if(  nrow != vec1.rows() )
    THROW_ERROR("eigen_blas::out_vecxf_diff: vec0.rows() != vec1.rows()");

  for(int irow=0; irow<nrow; irow++ ){
    const float val0 = vec0(irow);
    const float val1 = vec1(irow);
    const float diff = val1 - val0;
    if(val0==0.0f && val1==0.0f) continue;
    fprintf(fout,"%5d %7.4E %7.4E %7.4E\n",irow,val0,val1,diff);
  }
}

// output 2 vecxf to ascii file
void eigen_blas::out_vecxf_diff( const std::filesystem::path& pathout
, const Eigen::VectorXf &vec0, const Eigen::VectorXf &vec1 )
{
  LOG_INFO("for {}", pathout.string());
  FILE *fout = myapp::get_fout(pathout);
  out_vecxf_diff(fout,vec0,vec1);
  myapp::close(fout,pathout);
}


// output Eigen::MatrixXf
void eigen_blas::out_matrix( FILE *fout, const Eigen::MatrixXf &mat ){
  const size_t nrow = mat.rows();
  const size_t ncol = mat.cols();
  // colmajor output
  for(int icol=0; icol<ncol; icol++ ){
    for(int irow=0; irow<nrow; irow++ ){
      fprintf(fout,"%d %d %E\n",irow,icol,mat(irow,icol));
    }
  }
}

void eigen_blas::out_matrix(
  const std::filesystem::path& pathout, const Eigen::MatrixXf &mat )
{
  FILE *fout = myapp::get_fout(pathout);
  out_matrix(fout,mat);
  myapp::close(fout,pathout);
}

// output Eigen::MatrixXd
void eigen_blas::out_matrix(FILE *fout, const Eigen::MatrixXd &mat) {
  const size_t nrow = mat.rows();
  const size_t ncol = mat.cols();
  // colmajor output
  for (int icol = 0; icol < ncol; icol++) {
    for (int irow = 0; irow < nrow; irow++) {
      fprintf(fout, "%d %d %E\n", irow, icol, mat(irow, icol));
    }
  }
}

void eigen_blas::out_matrix(
  const std::filesystem::path& pathout, const Eigen::MatrixXd &mat)
{
  FILE *fout = myapp::get_fout(pathout);
  out_matrix(fout, mat);
  myapp::close(fout, pathout);
}

/// @brief output non zero element of Eigen::MatrixXf to ascii file
void eigen_blas::out_matxf(
 const std::filesystem::path& pathout, const Eigen::MatrixXf &mat)
{
  LOG_INFO("for {}", pathout.string());
  out_matrix_non_zero(pathout, mat);
}

/// @brief output non zero element of Eigen::MatrixXf to ascii file
/// @note this function output with row major
void eigen_blas::out_matxf_rowmajor(
  const std::filesystem::path& pathout, const Eigen::MatrixXf &mat)
{
  LOG_INFO("for {}", pathout.string());
  FILE *fout = myapp::get_fout(pathout);
  const int nrow = mat.rows();
  const int ncol = mat.cols();
  // output nrow and ncols
  fprintf(fout, "# %d %d output_row_major\n", nrow, ncol);
  // output non zero elements
  for (int irow = 0; irow < nrow; irow++) {
    Eigen::VectorXf vecxf_irow = mat.row(irow);
    for (int icol = 0; icol < ncol; icol++) {
      const float value = vecxf_irow(icol);
      if (fabs(value) < std::numeric_limits<float>::min()) continue;
      fprintf(fout, "%d %d %E\n", irow, icol, value);
    }
  }
  myapp::close(fout, pathout);
}

void eigen_blas::out_matrix_non_zero(
  FILE *fout, const Eigen::MatrixXf &mat, const float epsilon )
{
  // constexpr double epsilon_double = std::numeric_limits<double>::epsilon();
  const int nrow = mat.rows();
  const int ncol = mat.cols();
  // output nrow and ncols
  fprintf(fout,"# %d %d\n",nrow,ncol);
  // output non zero elements
  // colmajor output
  for(int icol=0; icol<ncol; icol++ ){
    Eigen::VectorXf vecxf_icol = mat.col(icol);
    for(int irow=0; irow<nrow; irow++ ){
      const float value = vecxf_icol(irow);
      if( fabs(value) < epsilon ) continue;
      fprintf(fout,"%d %d %E\n",irow,icol,value);
    }
  }
}

// output Eigen::MatrixXf
void eigen_blas::out_matrix_non_zero(
  const std::filesystem::path& pathout, const Eigen::MatrixXf &mat
  , const float epsilon )
{
  LOG_INFO("for {}", pathout.string());
  FILE *fout = myapp::get_fout(pathout);
  out_matrix_non_zero(fout, mat, epsilon);
  myapp::close(fout, pathout);
}

void eigen_blas::out_matrix_non_zero(
  FILE *fout, const Eigen::MatrixXd &mat, const double epsilon )
{
  const int nrow = mat.rows();
  const int ncol = mat.cols();
  // output nrow and ncols
  fprintf(fout, "# %d %d\n", nrow, ncol);
  // output non zero elements
  for (int icol = 0; icol < ncol; icol++) {
    Eigen::VectorXd vecxd_icol = mat.col(icol);
    for (int irow = 0; irow < nrow; irow++) {
      const double value = vecxd_icol(irow);
      if (fabs(value) < epsilon) continue;
      fprintf(fout, "%d %d %E\n", irow, icol, value);
    }
  }
}

void eigen_blas::out_matrix_non_zero(
  const std::filesystem::path& pathout, const Eigen::MatrixXd &mat
  , const double epsilon )
{
  FILE *fout = myapp::get_fout(pathout);
  out_matrix_non_zero(fout, mat, epsilon);
  myapp::close(fout, pathout);
  fprintf(stderr, "eigen_blas::out_matrix_non_zero ... done.\n");
}

// output non zero element of SpMatf to FILE*
void eigen_blas::out_matrix_non_zero(
  FILE *fout, const SpMatf &mat, const float epsilon )
{
  const int nrow = mat.rows();
  const int ncol = mat.cols();
  // output nrow and ncols
  fprintf(fout, "# %d %d\n", nrow, ncol);
  // output non zero elements with colmajor
  for (int icol = 0; icol < ncol; icol++) {
    for (SpMatf::InnerIterator it(mat, icol); it; ++it) {
      const int irow = it.row();
      const float value = it.value();
      if (fabs(value) < epsilon) continue;
      fprintf(fout, "%d %d %E\n", irow, icol, value);
    }
  }
}

// output non zero element of SpMatf to ascii file
void eigen_blas::out_matrix_non_zero(
  const std::filesystem::path& pathout, const SpMatf &mat, const float epsilon)
{
  LOG_INFO("for {}", pathout.string());
  FILE *fout = myapp::get_fout(pathout);
  out_matrix_non_zero(fout, mat, epsilon);
  myapp::close(fout, pathout);
  fprintf(stderr, "eigen_blas::out_matrix_non_zero ... done.\n");
}
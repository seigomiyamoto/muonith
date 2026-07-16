// this is io_binary.cpp
#include "ns_io_binary.hpp"
#include "ns_mymacro.hpp"
#include "ns_mylogger.hpp"
#include "spdlog_pch.hpp"
#include <cstdint>
#include <cstring>
#ifdef _MSC_VER
#include <intrin.h>
#endif

constexpr size_t MAX_STRING_SIZE = 1024 * 1024; // Set appropriate upper limit, e.g., 1MB

/// @brief Binary save/load function group
namespace io_binary {

// Determine if a file is binary or text
bool is_binary_file(const std::filesystem::path& path_in)
{
  // 1. File opening: open in binary mode
  std::ifstream file(path_in, std::ios::in | std::ios::binary);
  if (!file.is_open()) {
    return false;
  }

  // 2. Check file size and data consistency:
  //    Get file size, and if size is 0 (empty file), consider it unreadable
  file.seekg(0, std::ios::end);
  std::streampos fileSize = file.tellg();
  if (fileSize == 0) {
    return false;  // Empty file is not treated as binary (or handle as needed)
  }
  file.seekg(0, std::ios::beg);

  // 3. Adjust read position:
  //    Read required bytes (here max 512 bytes) and check current read position and bytes read
  char buffer[512];
  file.read(buffer, sizeof(buffer));
  std::streamsize bytesRead = file.gcount();
  // Get current read position (for debugging, can be adjusted if needed)
  std::streampos currentPos = file.tellg();
  (void)currentPos;  // Not currently used, but kept for verification

  // If read data contains control characters (other than newline/tab), judge as binary
  for (std::streamsize i = 0; i < bytesRead; ++i) {
    if (buffer[i] < 0x20 &&
        buffer[i] != '\r' &&
        buffer[i] != '\n' &&
        buffer[i] != '\t')
    {
      return true;
    }
  }
  return false;
}

std::ifstream open_ifstream(const std::filesystem::path& path)
{
  // Open file in binary mode
  std::ifstream ifs(path, std::ios::in | std::ios::binary);
  if (!ifs.is_open()) {
    THROW_ERROR("io_binary::open_ifstream: Cannot open file. path={}", path.string());
  }

  // Check file size (error if empty)
  ifs.seekg(0, std::ios::end);
  std::streampos fileSize = ifs.tellg();
  if (fileSize == 0) {
    THROW_ERROR("io_binary::open_ifstream: File is empty. path={}", path.string());
  }
  ifs.seekg(0, std::ios::beg);  // Reset read position to beginning

  return ifs; // Can be returned via move semantics in C++11 and later
}

std::ofstream open_ofstream(const std::filesystem::path& path) {
  std::ofstream ofs(path, std::ios::out | std::ios::binary);
  if (!ofs.is_open()) {
    THROW_ERROR("io_binary::open_ofstream: Cannot open file. path={}", path.string());
  }
  return ofs;
}

// Special type handling
void write_bool(std::ostream& os, bool value)
{
  char byte_value = value ? 1 : 0;
  write_binary(os, byte_value);
}

bool read_bool(std::istream& is)
{
  char byte_value = read_binary<char>(is);
  return byte_value != 0;
}

void write_vec_bool(std::ostream& os, const std::vector<bool>& vec)
{
  uint32_t size = static_cast<uint32_t>(vec.size());
  write_binary(os, size);
  static_assert(std::is_same<decltype(vec.size()), size_t>::value, "vec.size() must be size_t");
  for (bool value : vec) {
    write_binary(os, value);
  }
}

std::vector<bool> read_vec_bool(std::istream& is)
{
  // LOG_DEBUG("Before reading vec_bool size, stream pos = {}", static_cast<long long>(is.tellg()));
  uint32_t size = read_binary<uint32_t>(is);
  // LOG_DEBUG("Read vec_bool size = {}", size);
  std::vector<bool> vec(size);
  static_assert(std::is_same<decltype(vec.size()), size_t>::value, "vec.size() must be size_t");
  for (uint32_t i = 0; i < size; ++i) {
    vec.at(i) = read_bool(is);
  }
  return vec;
}

std::vector<std::vector<bool>> read_vec_vec_bool(std::istream& is)
{
  size_t outer_size = read_binary<size_t>(is);
  std::vector<std::vector<bool>> vec(outer_size);
  for (size_t i = 0; i < outer_size; ++i) {
    vec[i] = read_vec_bool(is);
  }
  return vec;
}

void write_vec_vec_bool(std::ostream& os, const std::vector<std::vector<bool>>& vec)
{
  write_binary(os, vec.size());
  for (const auto& inner_vec : vec) {
    write_vec_bool(os, inner_vec);
  }
}

void write_string(std::ostream& os, const std::string& str) {
  if (str.size() > MAX_STRING_SIZE) {
    THROW_ERROR("io_binary::write_string: String size too large. size={}", str.size());
  }
  write_binary(os, str.size());
  os.write(str.data(), str.size());
  if (!os) THROW_ERROR("io_binary::write_string: Failed to write string to binary.");
}

std::string read_string(std::istream& is) {
  size_t size = read_binary<size_t>(is);
  if (size > MAX_STRING_SIZE) {
    THROW_ERROR("io_binary::read_string: String size is too large. size={}", size);
  }
  std::string str(size, '\0');
  if (!is.read(&str[0], size)) {
    THROW_ERROR("io_binary::read_string: Failed to read string from binary.");
  }
  return str;
}

// write std::filesystem::path to binary std::ostream& os
void write_path(std::ostream& os, const std::filesystem::path& path) {
  write_string(os, path.string());
}

// read std::filesystem::path from binary std::istream& ifs
std::filesystem::path read_path(std::istream& is) {
  return std::filesystem::path(read_string(is));
}

// write std::map<int, std::array<int,2>> to binary std::ostream& os
void write_multimap_int_int2(std::ostream& os, const std::multimap<int, std::array<int, 2>>& map)
{
  write_binary(os, map.size());
  for (const auto& [key, array] : map) {
    write_binary(os, key);
    write_binary(os, array[0]);
    write_binary(os, array[1]);
  }
}

// read std::map<int, std::array<int,2>> from binary std::istream& is
std::multimap<int, std::array<int, 2>> read_multimap_int_int2(std::istream& is)
{
  size_t size = read_binary<size_t>(is);
  std::multimap<int, std::array<int, 2>> map;
  for (size_t i = 0; i < size; ++i) {
    int key = read_binary<int>(is);
    int value0 = read_binary<int>(is);
    int value1 = read_binary<int>(is);
    map.insert({key, {value0, value1}});
  }
  return map;
}

// write std::map<std::array<int,2>, int> to binary std::ostream& os
void write_map_int2_int(std::ostream& os, const std::map<std::array<int, 2>, int>& map)
{
  write_binary(os, map.size());
  for (const auto& [array, value] : map) {
    write_binary(os, array[0]);
    write_binary(os, array[1]);
    write_binary(os, value);
  }
}

// read std::map<std::array<int,2>, int> from binary std::istream& is
std::map<std::array<int, 2>, int> read_map_int2_int(std::istream& is)
{
  size_t size = read_binary<size_t>(is);
  std::map<std::array<int, 2>, int> map;
  for (size_t i = 0; i < size; ++i) {
    int key0 = read_binary<int>(is);
    int key1 = read_binary<int>(is);
    int value = read_binary<int>(is);
    map[{key0, key1}] = value;
  }
  return map;
}

// Eigen data structure handling
// write Eigen::Vector3d to std::ostream& os
void write_vec3d(std::ostream& os, const Eigen::Vector3d& vec)
{
  write_tuple(os, std::make_tuple(vec.x(), vec.y(), vec.z()));
  if (os.fail()) {
    THROW_ERROR("io_binary::write_vec3d: Failed to write Eigen::Vector3d to binary file.");
  }
}

// read Eigen::Vector3d from std::istream& ifs
Eigen::Vector3d read_vec3d(std::istream& is)
{
  auto [x, y, z] = read_tuple<double, double, double>(is);
  return Eigen::Vector3d(x, y, z);
}

// Generic read of Eigen::VectorXf from file stream
Eigen::VectorXf read_vecxf_stream(std::istream& is)
{
  size_t size = read_binary<size_t>(is);
  Eigen::VectorXf vec(size);
  if (!is.read(reinterpret_cast<char*>(vec.data()), sizeof(float) * size)) {
    THROW_ERROR("io_binary::read_vecxf_stream: Failed to read Eigen::VectorXf from stream.");
  }
  return vec;
}

// Read Eigen::VectorXf from file path
Eigen::VectorXf read_vecxf_bin(const std::filesystem::path& path_in)
{
  std::ifstream ifs = open_ifstream(path_in);  
  return read_vecxf_stream(ifs);
}

// Write Eigen::VectorXf to file stream
void write_vecxf_stream(std::ostream& os, const Eigen::VectorXf& vec) {
  write_binary(os, static_cast<size_t>(vec.size())); // Write size
  os.write(reinterpret_cast<const char*>(vec.data()), sizeof(float) * vec.size());
  if (!os) {
    THROW_ERROR("io_binary::write_vecxf_stream: Failed to write Eigen::VectorXf to stream.");
  }
}

// Write Eigen::VectorXf to file path
void out_vecxf_bin(const std::filesystem::path& pathout, const Eigen::VectorXf& vecxf)
{
  std::ofstream ofs = open_ofstream(pathout);
  write_vecxf_stream(ofs, vecxf);
}

// Generic write of Eigen::MatrixXf to binary file
void write_matxf_stream(std::ostream& os, const Eigen::MatrixXf& mat)
{
  write_binary(os, static_cast<size_t>(mat.rows())); // Write row count
  write_binary(os, static_cast<size_t>(mat.cols())); // Write column count
  os.write(reinterpret_cast<const char*>(mat.data()), sizeof(float) * mat.size());
  if (!os) THROW_ERROR("Failed to write Eigen::MatrixXf to stream.");
}

// Generic read of Eigen::MatrixXf from binary file
Eigen::MatrixXf read_matxf_stream(std::istream& is)
{
  std::size_t rows = read_binary<std::size_t>(is); // Read row count
  std::size_t cols = read_binary<std::size_t>(is); // Read column count

  Eigen::MatrixXf mat(rows, cols);
  is.read(reinterpret_cast<char*>(mat.data()), sizeof(float) * mat.size());
  if (!is) THROW_ERROR("Failed to read Eigen::MatrixXf from stream.");
  return mat;
}

// Write Eigen::MatrixXf to binary file with filepath
void out_matxf_bin(const std::filesystem::path& pathout, const Eigen::MatrixXf& mat)
{
  std::ofstream ofs = open_ofstream(pathout);
  write_matxf_stream(ofs, mat);
}

// Read Eigen::MatrixXf from binary file
Eigen::MatrixXf read_matxf_bin(const std::filesystem::path& path_in)
{
  std::ifstream ifs = open_ifstream(path_in);
  return read_matxf_stream(ifs);
}

// Write std::vector<Eigen::MatrixXf> to stream
void write_vec_matxf(std::ostream& os, const std::vector<Eigen::MatrixXf>& vec_mat)
{
  size_t size = vec_mat.size();
  if (size == 0) {
    THROW_ERROR("io_binary::write_vec_matxf: Vector of Eigen::MatrixXf is empty.");
  }

  // Write vector size
  io_binary::write_binary(os, size);

  // Write each matrix
  for (const auto& mat : vec_mat) {
    io_binary::write_matxf_stream(os, mat);
  }

  if (!os) {
    THROW_ERROR("io_binary::write_vec_matxf: Failed to write std::vector<Eigen::MatrixXf> to stream.");
  }
}

// Read std::vector<Eigen::MatrixXf> from stream
std::vector<Eigen::MatrixXf> read_vec_matxf(std::istream& is)
{
  size_t size = io_binary::read_binary<size_t>(is);

  std::vector<Eigen::MatrixXf> vec_mat;
  vec_mat.reserve(size);

  for (size_t i = 0; i < size; ++i) {
    vec_mat.push_back(io_binary::read_matxf_stream(is));
  }

  if (!is) {
    THROW_ERROR("io_binary::read_vec_matxf: Failed to read std::vector<Eigen::MatrixXf> from stream.");
  }

  return vec_mat;
}

// Write std::vector<Eigen::MatrixXf> to file path
void write_vec_matxf(const std::filesystem::path& path,
                     const std::vector<Eigen::MatrixXf>& vec_mat)
{
  std::ofstream ofs = open_ofstream(path);
  io_binary::write_vec_matxf(ofs, vec_mat);

  ofs.close();
  if (!ofs) {
    THROW_ERROR("io_binary::write_vec_matxf: Failed to write std::vector<Eigen::MatrixXf> to file. path={}", path.string());
  }
}

// Read std::vector<Eigen::MatrixXf> from file path
std::vector<Eigen::MatrixXf> read_vec_matxf(const std::filesystem::path& path)
{
  std::ifstream ifs = open_ifstream(path);
  std::vector<Eigen::MatrixXf> vec_mat = io_binary::read_vec_matxf(ifs);

  ifs.close();
  if (!ifs) {
    THROW_ERROR("io_binary::read_vec_matxf: Failed to read std::vector<Eigen::MatrixXf> from file. path={}", path.string());
  }

  return vec_mat;
}

void write_spmatf_stream(std::ostream& os, const SpMatf& mat)
{
  // Write row count, column count, and number of non-zero elements
  write_binary(os, static_cast<size_t>(mat.rows()));
  write_binary(os, static_cast<size_t>(mat.cols()));
  write_binary(os, static_cast<size_t>(mat.nonZeros()));

  // Write (row, col, value) for each non-zero element
  for (int k = 0; k < mat.outerSize(); ++k) {
    for (Eigen::SparseMatrix<float>::InnerIterator it(mat, k); it; ++it) {
      int r = it.row();
      int c = it.col();
      float v = it.value();
      write_binary(os, r);
      write_binary(os, c);
      write_binary(os, v);
    }
  }

  if (!os) THROW_ERROR("Failed to write Eigen::SparseMatrix<float> to stream.");
}

SpMatf read_spmatf_stream(std::istream& is)
{
  std::size_t rows = read_binary<std::size_t>(is);
  std::size_t cols = read_binary<std::size_t>(is);
  std::size_t nnz  = read_binary<std::size_t>(is);

  std::vector<Eigen::Triplet<float>> triplets;
  triplets.reserve(nnz);

  for (std::size_t i = 0; i < nnz; ++i) {
    int r = read_binary<int>(is);
    int c = read_binary<int>(is);
    float v = read_binary<float>(is);
    triplets.emplace_back(r, c, v);
  }

  SpMatf mat(rows, cols);
  mat.setFromTriplets(triplets.begin(), triplets.end());
  mat.makeCompressed();

  if (!is) THROW_ERROR("Failed to read Eigen::SparseMatrix<float> from stream.");
  return mat;
}

void write_vec_spmatf(std::ostream& os, const std::vector<SpMatf>& vec_mat)
{
  // First write the vector size
  write_binary(os, static_cast<size_t>(vec_mat.size()));
  // Write each sparse matrix in sequence
  for (const auto& mat : vec_mat) {
    write_spmatf_stream(os, mat);
  }
}

/// @brief Reads a std::vector<SpMatf> from a general-purpose binary file
/// @param is Input stream
/// @return Vector of sparse matrices that was read
std::vector<SpMatf> read_vec_spmatf(std::istream& is)
{
  std::size_t vec_size = read_binary<size_t>(is);
  std::vector<SpMatf> vec_mat;
  vec_mat.reserve(vec_size);
  for (std::size_t i = 0; i < vec_size; ++i) {
    vec_mat.push_back(read_spmatf_stream(is));
  }
  return vec_mat;
}

/// @brief Function that writes a std::vector<SpMatf> to a file path
/// @param path Output file path
/// @param vec_mat Vector of sparse matrices to write
void write_vec_spmatf(const std::filesystem::path& path,
                      const std::vector<SpMatf>& vec_mat)
{
  std::ofstream ofs = open_ofstream(path);
  io_binary::write_vec_spmatf(ofs, vec_mat);

  ofs.close();
  if (!ofs) {
    THROW_ERROR("io_binary::write_vec_spmatf: Failed to write std::vector<Eigen::SparseMatrix<float>> to file. path={}", path.string());
  }
}

std::vector<SpMatf> read_vec_spmatf(const std::filesystem::path& path)
{
  std::ifstream ifs = open_ifstream(path);
  std::vector<SpMatf> vec_mat = io_binary::read_vec_spmatf(ifs);
  ifs.close();
  if (!ifs) {
    THROW_ERROR("io_binary::read_vec_spmatf: Failed to read std::vector<Eigen::SparseMatrix<float>> from file. path={}", path.string());
  }
  return vec_mat;
}

// @brief write std::vector<std::tuple<int,double>> to std::ostream& os
void write_vec_tp_int_double(
  std::ostream& os, const std::vector<std::tuple<int,double>>& vec)
{
  write_binary(os, vec.size());
  for (const auto& tp : vec) {
    write_tuple(os, tp);
  }
}

std::vector<std::tuple<int,double>>
  read_vec_tp_int_double(std::istream& is)
{
  size_t size = read_binary<size_t>(is);
  std::vector<std::tuple<int,double>> vec(size);
  for (size_t i = 0; i < size; ++i) {
    vec.at(i) = read_tuple<int,double>(is);
  }
  return vec;
}

void write_uomap_int_double(
  std::ostream& os, const std::unordered_map<int,double>& map)
{
  write_binary(os, map.size());
  for (const auto& [key, value] : map) {
    write_binary(os, key);
    write_binary(os, value);
  }
}

// @brief read std::unordered_map<int,double> from std::istream& is
std::unordered_map<int,double> read_uomap_int_double(std::istream& is)
{
  size_t size = read_binary<size_t>(is);
  std::unordered_map<int,double> map;
  for (size_t i = 0; i < size; ++i) {
    int key = read_binary<int>(is);
    double value = read_binary<double>(is);
    map[key] = value;
  }
  return map;
}

// write set_int to ostream 
void write_set_int(std::ostream& os, const std::set<int>& set_int) {
  write_binary(os, set_int.size()); // Write element count
  for (const auto& value : set_int) {
    write_binary(os, value); // Write each element
  }
}

// read set_int from istream
std::set<int> read_set_int(std::istream& is) {
  size_t size = read_binary<size_t>(is); // Read element count
  std::set<int> set_int;
  for (size_t i = 0; i < size; ++i) {
    set_int.insert(read_binary<int>(is)); // Read and insert element
  }
  return set_int;
}

void write_vec_tp_bool_double(std::ostream& os, const std::vector<std::tuple<bool,double>>& vec)
{
  write_binary(os, vec.size());
  for (const auto& t : vec) {
    write_bool(os, std::get<0>(t));
    write_binary(os, std::get<1>(t));
  }
}

std::vector<std::tuple<bool,double>> read_vec_tp_bool_double(std::istream& is)
{
  // Read saved size
  size_t size = read_binary<size_t>(is);
  std::vector<std::tuple<bool,double>> vec;
  vec.reserve(size);
  // Read each tuple in sequence
  for (size_t i = 0; i < size; ++i) {
    bool b = read_bool(is);         // Read bool
    double d = read_binary<double>(is);  // Read double
    vec.emplace_back(b, d);
  }
  return vec;
}

// unnamed namespace
namespace {
  /// @brief Determine if system is little endian
  bool is_little_endian() {
    /// Simple implementation
    union {
      uint32_t i;
      char c[4];
    } test = { 0x01020304 };

    return (test.c[0] == 0x04);
  }

  /// @brief Infer CPU architecture name
  std::string get_architecture_name() {
#if defined(__x86_64__) || defined(_M_X64)
    return "x86_64";
#elif defined(__i386) || defined(_M_IX86)
    return "x86";
#elif defined(__arm__) || defined(_M_ARM)
    return "arm";
#elif defined(__aarch64__) || defined(_M_ARM64)
    return "arm64";
#elif defined(__powerpc64__) || defined(_M_PPC)
    return "powerpc64";
#else
    return "unknown";
#endif
  }

  /// @brief Pointer size in bytes
  uint16_t get_pointer_size() {
    // For example, 8 for 64-bit, 4 for 32-bit
    // Generally determined by sizeof(void*)
    return static_cast<uint16_t>(sizeof(void*));
  }
} // unnamed namespace

//==========================================================
// 1) Get ArchitectureInfo of the current runtime environment
ArchitectureInfo get_current_architecture_info()
{
  ArchitectureInfo info;
  info.architectureName = get_architecture_name();
  info.isLittleEndian   = is_little_endian();
  info.pointerSize      = get_pointer_size();
  return info;
}

//==========================================================
// 2) Write ArchitectureInfo to binary
void write_architecture_info(std::ostream& os, const ArchitectureInfo& info)
{
  // architectureName ( std::string )
  write_string(os, info.architectureName);
  
  // isLittleEndian ( bool )
  write_bool(os, info.isLittleEndian);
  
  // pointerSize ( uint16_t )
  write_binary(os, info.pointerSize);
}

//==========================================================
// 3) Read ArchitectureInfo from binary
ArchitectureInfo read_architecture_info(std::istream& is)
{
  ArchitectureInfo info;
  info.architectureName = read_string(is);
  info.isLittleEndian   = read_bool(is);
  info.pointerSize      = read_binary<uint16_t>(is);
  return info;
}

//==========================================================
// 4) Compare read ArchitectureInfo with current runtime environment,
//    throw exception if mismatch
void check_architecture_compatibility_or_throw(const ArchitectureInfo& fileInfo)
{
  ArchitectureInfo current = get_current_architecture_info();

  // Check if architecture name, endianness, and pointer size all match
  if (fileInfo != current) {
    std::string err_msg =
      "Incompatible architecture:\n"
      "  - file:    " + fileInfo.architectureName
                    + ", " + (fileInfo.isLittleEndian ? "LE" : "BE")
                    + ", pointerSize=" + std::to_string(fileInfo.pointerSize) + "\n"
      "  - current: " + current.architectureName
                    + ", " + (current.isLittleEndian ? "LE" : "BE")
                    + ", pointerSize=" + std::to_string(current.pointerSize) + "\n";
    LOG_ERROR(err_msg);
    THROW_ERROR("io_binary::check_architecture_compatibility_or_throw: {}", err_msg);
  }
}

//==========================================================
// Write std::vector<Eigen::VectorXf> to binary stream
void write_vec_vecxf(std::ostream& os, const std::vector<Eigen::VectorXf>& vec)
{
  write_binary(os, static_cast<uint64_t>(vec.size()));
  for (const auto& v : vec) {
    write_vecxf_stream(os, v);
  }
}

//==========================================================
// Read std::vector<Eigen::VectorXf> from binary stream
std::vector<Eigen::VectorXf> read_vec_vecxf(std::istream& is)
{
  uint64_t n = read_binary<uint64_t>(is);
  std::vector<Eigen::VectorXf> vec(n);
  for (auto& v : vec) {
    v = read_vecxf_stream(is);
  }
  return vec;
}

//==========================================================
// Write std::vector<Eigen::VectorXf> to file path
void write_vec_vecxf(const std::filesystem::path& path,
                     const std::vector<Eigen::VectorXf>& vec)
{
  std::ofstream ofs = open_ofstream(path);
  io_binary::write_vec_vecxf(ofs, vec);

  ofs.close();
  if (!ofs) {
    THROW_ERROR("io_binary::write_vec_vecxf: Failed to write std::vector<Eigen::VectorXf> to file. path={}", path.string());
  }
}

//==========================================================
// Read std::vector<Eigen::VectorXf> from file path
std::vector<Eigen::VectorXf> read_vec_vecxf(const std::filesystem::path& path_in)
{
  std::ifstream ifs = open_ifstream(path_in);
  std::vector<Eigen::VectorXf> vec = io_binary::read_vec_vecxf(ifs);
  ifs.close();
  if (!ifs) {
    THROW_ERROR("io_binary::read_vec_vecxf: Failed to read std::vector<Eigen::VectorXf> from file. path={}", path_in.string());
  }
  return vec;
}

} // namespace io_binary end

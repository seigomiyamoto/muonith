/// @file ns_io_binary.hpp
/// @brief Binary file I/O utilities
/// @details Namespace for binary serialization and deserialization of data structures.
#pragma once

#include <iostream>
#include <set>
#include <fstream>
#include <stdexcept>
#include <vector>
#include <tuple>
#include <string>
#include <filesystem>
#include <Eigen/Dense>
#include <Eigen/Sparse>
#include "ns_mymacro.hpp"

namespace io_binary {
  /// @brief Type alias for Eigen::SparseMatrix<float>
  using SpMatf = Eigen::SparseMatrix<float>;

  //=================================================================
  /// @name binary_io "Binary I/O function group"
  /// @{

  /// @brief Determine if a file is binary or text
  /// @param path_in Input file path
  /// @return true if binary file, false if text file
  bool is_binary_file(const std::filesystem::path& path_in);

  /// @brief Open a binary file for reading
  /// @param path File path
  /// @return std::ifstream Input file stream
  std::ifstream open_ifstream(const std::filesystem::path& path);

  /// @brief Open a binary file for writing
  /// @param path File path
  /// @return std::ofstream Output file stream
  std::ofstream open_ofstream(const std::filesystem::path& path);

  /// @brief Write generic binary data
  /// @tparam T Type of data to write
  /// @param os Output stream
  /// @param value Value to write
  /// @details usage example: write_binary(ofs, vec_xyz.size());
  template <typename T>
  void write_binary(std::ostream& os, const T& value) {
    os.write(reinterpret_cast<const char*>(&value), sizeof(T));
    if (!os) THROW_ERROR("io_binary::write_binary: Failed to write binary data.");
  };

  /// @brief Read generic binary data
  /// @tparam T Type of data to read
  /// @param is Input stream
  /// @return Read value
  /// @details usage example: size_t size = read_binary<size_t>(ifs);
  template <typename T>
  T read_binary(std::istream& is) {
    T value;
    if (!is.read(reinterpret_cast<char*>(&value), sizeof(T))) {
      fprintf(stderr, "Failed to read binary type %s data.\n", typeid(T).name());
      THROW_ERROR("io_binary::read_binary: Failed to read binary data.");
    }
    return value;
  };

  /// @brief Write bool value to binary file
  /// @param os Output stream
  /// @param value Bool value to write
  void write_bool(std::ostream& os, bool value);

  /// @brief Read bool value from binary file
  /// @param is Input stream
  /// @return Read bool value
  bool read_bool(std::istream& is);

  /// @brief Write string to binary file
  /// @param os Output stream
  /// @param str String to write
  void write_string(std::ostream& os, const std::string& str);

  /// @brief Read string from binary file
  /// @param is Input stream
  /// @return Read string
  std::string read_string(std::istream& is);

  /// @brief Write Eigen::Vector3d to output stream
  void write_vec3d(std::ostream& os, const Eigen::Vector3d& vec);

  /// @brief Read Eigen::Vector3d from binary file
  Eigen::Vector3d read_vec3d(std::istream& is);

  /// @brief Write Eigen::VectorXf to binary stream
  /// @param os Output stream
  /// @param vec Eigen::VectorXf to write
  void write_vecxf_stream(std::ostream& os, const Eigen::VectorXf& vec);

  /// @brief Read Eigen::VectorXf from binary stream
  /// @param is Input stream
  /// @return Read Eigen::VectorXf
  Eigen::VectorXf read_vecxf_stream(std::istream& is);

  /// @brief Write Eigen::MatrixXf to binary stream
  /// @param os Output stream
  /// @param mat Eigen::MatrixXf to write
  void write_matxf_stream(std::ostream& os, const Eigen::MatrixXf& mat);

  /// @brief Read Eigen::MatrixXf from binary stream
  /// @param is Input stream
  /// @return Read Eigen::MatrixXf
  Eigen::MatrixXf read_matxf_stream(std::istream& is);

  /// @brief Read Eigen::VectorXf from file path
  /// @param path_in Input file path
  /// @return Eigen::VectorXf
  Eigen::VectorXf read_vecxf_bin(const std::filesystem::path& path_in);

  /// @brief Write Eigen::VectorXf to file path
  /// @param pathout Output file path
  /// @param vec Vector to write
  void out_vecxf_bin(const std::filesystem::path& pathout, const Eigen::VectorXf& vec);

  /// @brief Read Eigen::MatrixXf from file path
  /// @param path_in Input file path
  /// @return Eigen::MatrixXf
  Eigen::MatrixXf read_matxf_bin(const std::filesystem::path& path_in);

  /// @brief Write Eigen::MatrixXf to file path
  void out_matxf_bin(const std::filesystem::path& pathout, const Eigen::MatrixXf& mat);

  /// @brief Write std::filesystem::path to binary file
  void write_path(std::ostream& os, const std::filesystem::path& path);

  /// @brief write std::map<int, std::array<int,2>> to binary std::ostream& os
  void write_multimap_int_int2(std::ostream& os
    , const std::multimap<int, std::array<int, 2>>& map);
  
  /// @brief Read std::multimap<int, std::array<int,2>> from binary stream
  std::multimap<int, std::array<int, 2>> read_multimap_int_int2(std::istream& is);

  /// @brief Write std::map<std::array<int,2>, int> to binary stream
  void write_map_int2_int(std::ostream& os
    , const std::map<std::array<int, 2>, int>& map);

  /// @brief Read std::map<std::array<int,2>, int> from binary stream
  std::map<std::array<int, 2>, int> read_map_int2_int(std::istream& is);

  /// @brief Read std::filesystem::path from binary file
  std::filesystem::path read_path(std::istream& is);

  /// @brief Write std::vector to binary file
  /// @tparam T Element type of vector
  /// @param os Output stream
  /// @param vec std::vector to write
  template <typename T>
  void write_vec(std::ostream& os, const std::vector<T>& vec) {
    write_binary(os, vec.size());
    for (const auto& element : vec) {
      write_binary(os, element);
    }
  };

  /// @brief Read std::vector from binary file
  /// @tparam T Element type of vector
  /// @param is Input stream
  /// @return Read std::vector
  template <typename T>
  std::vector<T> read_vec(std::istream& is) {
    size_t size = read_binary<size_t>(is);
    std::vector<T> vec(size);
    for (auto& element : vec) {
      element = read_binary<T>(is);
    }
    return vec;
  };

  /// @brief Write std::vector<bool> to binary file
  /// @param os Output stream
  /// @param vec std::vector<bool> to write
  void write_vec_bool(std::ostream& os, const std::vector<bool>& vec);

  /// @brief Read std::vector<bool> from binary file
  /// @param is Input stream
  /// @return Read std::vector<bool>
  std::vector<bool> read_vec_bool(std::istream& is);

  /// @brief Read 2D std::vector<bool> from binary file
  /// @param is Input stream
  /// @return Read 2D std::vector<bool>
  std::vector<std::vector<bool>> read_vec_vec_bool(std::istream& is);

  /// @brief Write 2D std::vector<bool> to binary file
  /// @param os Output stream
  /// @param vec 2D std::vector<bool> to write
  void write_vec_vec_bool(std::ostream& os, const std::vector<std::vector<bool>>& vec);

  /// @brief Write 2D std::vector to binary file
  /// @tparam T Element type of vector
  /// @param os Output stream
  /// @param vec 2D std::vector to write
  template <typename T>
  void write_vec_vec(std::ostream& os, const std::vector<std::vector<T>>& vec) {
    write_binary(os, vec.size());
    for (const auto& inner_vec : vec) {
      write_vec(os, inner_vec);
    }
  };

  /// @brief Read 2D std::vector from binary file
  /// @tparam T Element type of vector
  /// @param is Input stream
  /// @return Read 2D std::vector
  template <typename T>
  std::vector<std::vector<T>> read_vec_vec(std::istream& is) {
    size_t outer_size = read_binary<size_t>(is);
    std::vector<std::vector<T>> vec(outer_size);
    for (auto& inner_vec : vec) {
      inner_vec = read_vec<T>(is);
    }
    return vec;
  };

  /// @brief Write tuple to binary file
  /// @tparam Args Tuple element types
  /// @param os Output stream
  /// @param tuple Tuple to write
  template <typename... Args>
  void write_tuple(std::ostream& os, const std::tuple<Args...>& tuple) {
    std::apply([&os](const auto&... args) {
      (..., write_binary(os, args));
    }, tuple);
  };

  /// @brief Read tuple from binary file
  /// @tparam Args Tuple element types
  /// @param is Input stream
  /// @return Read tuple
  template <typename... Args>
  std::tuple<Args...> read_tuple(std::istream& is) {
    std::tuple<Args...> tuple;
    std::apply([&is](auto&... args) {
      (..., (args = read_binary<std::decay_t<decltype(args)>>(is)));
    }, tuple);
    return tuple;
  };


  /// @brief Write std::vector<Eigen::MatrixXf> to binary stream
  /// @param os Output stream
  /// @param vec_mat std::vector<Eigen::MatrixXf> to write
  void write_vec_matxf(std::ostream& os, const std::vector<Eigen::MatrixXf>& vec_mat);

  /// @brief Read std::vector<Eigen::MatrixXf> from binary stream
  /// @param is Input stream
  /// @return Read std::vector<Eigen::MatrixXf>
  std::vector<Eigen::MatrixXf> read_vec_matxf(std::istream& is);

  /// @brief Write std::vector<Eigen::VectorXf> to binary stream
  /// @param os Output stream
  /// @param vec std::vector<Eigen::VectorXf> to write
  void write_vec_vecxf(std::ostream& os, const std::vector<Eigen::VectorXf>& vec);

  /// @brief Read std::vector<Eigen::VectorXf> from binary stream
  /// @param is Input stream
  /// @return Read std::vector<Eigen::VectorXf>
  std::vector<Eigen::VectorXf> read_vec_vecxf(std::istream& is);

  /// @brief Write std::vector<Eigen::VectorXf> to file path.
  /// @param path Output file path.
  /// @param vec Vector of VectorXf to write.
  void write_vec_vecxf(const std::filesystem::path& path,
                      const std::vector<Eigen::VectorXf>& vec);

  /// @brief Read std::vector<Eigen::VectorXf> from file path.
  /// @param path_in Input file path.
  /// @return std::vector<Eigen::VectorXf>
  std::vector<Eigen::VectorXf> read_vec_vecxf(const std::filesystem::path& path_in);

  /// @brief Write std::vector<Eigen::MatrixXf> to file path
  /// @param path Output file path
  /// @param vec_mat Vector of matrices to write
  void write_vec_matxf(const std::filesystem::path& path,
                      const std::vector<Eigen::MatrixXf>& vec_mat);

  /// @brief Read std::vector<Eigen::MatrixXf> from file path
  /// @param path_in Input file path
  /// @return std::vector<Eigen::MatrixXf>
  std::vector<Eigen::MatrixXf> read_vec_matxf(const std::filesystem::path& path_in);

  /// @brief Write Eigen::SparseMatrix<float> to binary stream
  void write_spmatf_stream(std::ostream& os, const SpMatf& mat);

  /// @brief Read Eigen::SparseMatrix<float> from binary stream
  SpMatf read_spmatf_stream(std::istream& is);

  /// @brief Write std::vector<SpMatf> to binary stream
  /// @param os Output stream
  /// @param vec_mat Vector of sparse matrices to write
  void write_vec_spmatf(std::ostream& os, const std::vector<SpMatf>& vec_mat);

  /// @brief Read std::vector<SpMatf> from binary stream
  /// @param is Input stream
  /// @return Read vector of sparse matrices
  std::vector<SpMatf> read_vec_spmatf(std::istream& is);

  /// @brief Write std::vector<SpMatf> to file path
  /// @param path Output file path
  /// @param vec_mat Vector of sparse matrices to write
  void write_vec_spmatf(const std::filesystem::path& path,
    const std::vector<SpMatf>& vec_mat);

  /// @brief Read std::vector<SpMatf> from file path
  /// @param path Input file path
  /// @return Read vector of sparse matrices
  std::vector<SpMatf> read_vec_spmatf(const std::filesystem::path& path);

  /// @brief write std::vector<std::tuple<int,double>> to std::ostream& os
  void write_vec_tp_int_double(
    std::ostream& os, const std::vector<std::tuple<int,double>>& vec);
  
  /// @brief read std::vector<std::tuple<int,double>> from std::istream& is
  std::vector<std::tuple<int,double>>
    read_vec_tp_int_double(std::istream& is);

  /// @brief write std::unordered_map<int,double>>& to std::ostream& os
  void write_uomap_int_double(
    std::ostream& os, const std::unordered_map<int,double>& map);

  /// @brief read std::unordered_map<int,double> from std::istream& is
  std::unordered_map<int,double> read_uomap_int_double(std::istream& is);

  /// @brief Write std::set<int> to binary file
  /// @param os Output stream
  /// @param set_int std::set<int> to write
  void write_set_int(std::ostream& os, const std::set<int>& set_int);

  /// @brief Read std::set<int> from binary file
  /// @param is Input stream
  /// @return Read std::set<int>
  std::set<int> read_set_int(std::istream& is);

  /// @brief Write std::vector<std::tuple<bool,double>> to binary file
  void write_vec_tp_bool_double(std::ostream& os, const std::vector<std::tuple<bool,double>>& vec);

  /// @brief Read std::vector<std::tuple<bool,double>> from binary file
  std::vector<std::tuple<bool,double>> read_vec_tp_bool_double(std::istream& is);

  /// @} ------------------------------------------------------------------

  //==========================================================
  /// @brief Architecture information about the environment where the file was written
  struct ArchitectureInfo {
    std::string architectureName; ///< "x86", "x86_64", "arm", "arm64", etc.
    bool        isLittleEndian;   ///< true if little endian, false if big endian
    uint16_t    pointerSize;      ///< Pointer size in bytes. Example: 8 (for 64-bit)

    // Utility: equality check
    bool operator==(const ArchitectureInfo& other) const {
      return (architectureName == other.architectureName) &&
             (isLittleEndian   == other.isLittleEndian  ) &&
             (pointerSize      == other.pointerSize     );
    }

    bool operator!=(const ArchitectureInfo& other) const {
      return !(*this == other);
    }
  };

  //==========================================================
  /// @brief Get ArchitectureInfo of the current runtime environment
  ArchitectureInfo get_current_architecture_info();

  //==========================================================
  /// @brief Write ArchitectureInfo to binary file
  void write_architecture_info(std::ostream& os, const ArchitectureInfo& info);

  //==========================================================
  /// @brief Read ArchitectureInfo from binary file
  ArchitectureInfo read_architecture_info(std::istream& is);

  //==========================================================
  /// @brief Check if the ArchitectureInfo written in file matches current runtime environment
  /// @throws std::runtime_error if architectures are incompatible
  void check_architecture_compatibility_or_throw(const ArchitectureInfo& fileInfo);

} // namespace io_binary

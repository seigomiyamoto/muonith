/// @file ns_myapp.hpp
/// @brief Application-wide utility functions
/// @details Namespace for common application utilities including random number generation and statistics.
#pragma once

// eigen library

#include <cstdio>
#include <string>
#include <cmath>
#include <iostream>
#include <fstream>
#include <chrono>
#include <thread>

#define _USE_MATH_DEFINES // for use M_PI, this order is required.
#include <cmath>

// std::map<>
#include <map>
#include <unordered_map>
#include <algorithm>
#include <cassert>
#include <vector>
#include <filesystem> // for std::filesystem::exists
#include <limits>
#include <tuple>

#include <stdexcept>  // required for std::runtime_error
#include <type_traits> // for std::is_same_v

// spdlog logger
#include "spdlog_pch.hpp"
#include "ns_mylogger.hpp"
#include <stdexcept>

// json
#include <nlohmann/json.hpp>
using json = nlohmann::json;

#include "ns_io_binary.hpp"
#include <Eigen/Dense>
#include <Eigen/Sparse>
#include "cls_Angle.hpp"
#include "cls_AABB.hpp"
#include "cls_Ray.hpp"
#include "ns_tuple_int.hpp"
#include "ns_mymacro.hpp"
#include "ns_mylogger.hpp"
#include "ns_type_definitions.hpp"
#include "ns_angle_util.hpp"

using namespace index_type_definitions;


// Static assertion to verify that Eigen::MatrixXf is ColumnMajor
static_assert((Eigen::MatrixXf::Flags & Eigen::RowMajorBit) == 0, "MatrixXf must be ColumnMajor");

/// @brief Enum class representing the order of a vector.
enum class VectorOrder {
  NotOrdered,
  Increasing,
  Decreasing
};

/// @brief Enum class representing monotonic status of a dataset
enum class MonotonicStatus {
    Increasing,
    Decreasing,
    NotMonotonic
};

/// @brief Type definition for sparse matrix
using SpMatf = Eigen::SparseMatrix<float>;

//#####################################################################################
//#####################################################################################
/// @namespace myapp
/// @brief Namespace for miscellaneous utility functions
/// @ingroup basicTools
//#####################################################################################
//#####################################################################################
namespace myapp {
  //=======================================================================
  /// @name constant_myapp
  /// @brief Constants
  /// @{
  
  /// @brief default delimiter for split 
  static constexpr char char_delim_default = ' ';
  
  /// @brief default commentout character when read ascii file
  static constexpr char char_commentout_default = '#';
  

  ///@} ------------------------------------------------------------------
  

  //=======================================================================
  /// @name basic_myapp
  /// @brief Basic utility functions
  /// @{

  /// @brief open file and get fout with safety
  /// @ingroup basicTools
  FILE* get_fout( const std::filesystem::path& pathout );

  /// @brief open file for binary output with safety checks
  FILE* get_fout_binary( const std::filesystem::path &path_out );

  /// @brief close file with safety
  /// @ingroup basicTools
  /// @hidecallgraph
  /// @hidecallergraph
  void close( FILE* fout, const std::filesystem::path& pathout );

  /// @brief check file exists or not
  /// @ingroup basicTools
  /// @hidecallgraph
  /// @hidecallergraph
  void filecheck(const std::filesystem::path& pathout );

  /// @brief Round to n decimal places
  float  round_n( float  number, const double n);
  /// @brief Round to n decimal places
  double round_n( double number, const double n);

  /// @brief Get random value from Poisson distribution
  int poisson(const double mean);

  // Count digits
  // https://www.delftstack.com/ja/howto/cpp/number-of-digits-in-a-number-cpp/
  // size_t count_digits(const double number);

  // chatGPT4 2023-04-10 16:43:06
  /// @brief Count decimal places
  int countDecimalPlaces(double number, const int seido_keta=8);

  /// @brief log10 that returns error if x<0
  double log10(const double x);

  /// @brief Convert number to string with specified digit format
  std::string formatNumber(int digits, int number);

  /// @brief Split long string containing char_delim into multiple strings
  /// @details https://qiita.com/iseki-masaya/items/70b4ee6e0877d12dafa8
  std::vector<std::string> split_naive(const std::string &s
  , const char char_delim=char_delim_default);

  /// @brief Avoid returning empty strings even if multiple char_delim are consecutive. Higher performance than split_naive.
  std::vector< std::string > split(const std::string &s
  , const char char_delim=char_delim_default
  , const char char_commentout=char_commentout_default);

  /// @brief bool_zero or not for double
  bool is_zero( const double value, const double factor );

  // @brief return distance between 2 Eigen::Vector3d
  // double distance( const Eigen::Vector3d v3_pos1, const Eigen::Vector3d v3_pos2 );

  /// @brief Determine whether given dataset is monotonically increasing or decreasing
  /// @details "Monotonic" here means values change in only one direction. \n
  /// enum class MonotonicStatus { Increasing, Decreasing, NotMonotonic};
  MonotonicStatus getMonotonicStatus(
    const std::vector<std::pair<double, double>>& data);

  ///@brief Enum class representing the order of a vector.
  ///@details The possible values are: \n
  ///- Increasing: the vector is in increasing order. \n
  ///- Decreasing: the vector is in decreasing order. \n
  ///- NotOrdered: the vector is not ordered. 
  VectorOrder getVectorOrder(const std::vector<double>& vec);

  /// @brief ump search function
  template<typename K, typename V>
  inline std::pair<K,V> getMaximumKey(const std::unordered_map<K,V> &ump) {
    return *std::max_element(ump.begin(), ump.end(), [](std::pair<K,V> const &x, std::pair<K,V> const &y) {
      return x.first < y.first;
    });
  };

  /// @brief Check continuity of vector<int>.
  /// @return true:continuous false:not continuous
  /// @note Empty vector is considered continuous.
  bool sort_and_check_vec_int_continuity( std::vector<int>& vec_int);

  ///@brief Function to release vector memory.
  /// @note Usage example \n
  /// std::vector<int> myVector = ... \n
  /// clearVectorMemory(myVector);
  template <typename T>
  void clearVectorMemory(std::vector<T>& v) {
    v.clear();  // Clear vector elements (unnecessary?)
    v.shrink_to_fit();  // Attempt to release memory?
    std::vector<T>().swap(v);  // Ensure memory release by swapping with empty vector
  };

  ///@} ------------------------------------------------------------------


  //=======================================================================
  /// @name read_data_myapp
  /// @details Data reading functions
  /// @{

  /// @brief read (double,double) array from file
  std::vector<std::pair<double,double>> read_double2( const std::filesystem::path &path_file );

  /// @brief read 3clm ascii data using FILE* and sscanf
  std::vector<Eigen::Vector3d> read_xyz( const std::filesystem::path &path_file );

  /// @brief read 3clm ascii data using std::ifstream and myapp::split
  std::vector<Eigen::Vector3d> read_xyz_stream( const std::filesystem::path &path_file );

  /// @brief read 3-column ASCII xyz data. Binary point-cloud input is no longer supported (throws std::runtime_error).
  std::vector< std::array<double,3> >
    read_vec_xyz( const std::filesystem::path &path_in );

  /// @brief read 3clm ascii data
  std::vector< std::array<double,3> >
    read_vec_xyz_txt( const std::filesystem::path &path_file );
  
  /// @brief read 3double + 3int ascii data
  std::vector< std::tuple<double,double,double,int,int,int> >
    read_vec_xyzrgb( const std::filesystem::path &path_file );

  ///@} ------------------------------------------------------------------

  //=======================================================================
  /// @name output_myapp
  /// @details Miscellaneous output functions
  /// @{

  /// @brief output std::vector<Eigen::Vector3d>
  void out_vec_v3( const std::filesystem::path& pathout
  , const std::vector<Eigen::Vector3d> &vec_v3 );

  /// @brief out_vec_double to ascii file
  void out( const std::filesystem::path& pathout, const std::vector<double> &vec );

  /// @brief output vecxf to ascii file
  void out( const std::filesystem::path& pathout, const Eigen::VectorXf &vecxf );

  /// @brief output std::vector< std::array<double,3> > to ascii
  void out( const std::filesystem::path& pathout
    , const std::vector< std::array<double,3>> &vec_tp3 );

  /// @brief output std::vector< std::array<double,4> > to ascii
  void out( const std::filesystem::path& pathout
    , const std::vector< std::array<double,4>> &vec_tp4 );

  /// @brief output std::vector< std::tuple<double,double,double,double,int,int,int> > to ascii
  void out_vec_xyzdzrgb( const std::filesystem::path& pathout
    , const std::vector< std::tuple<double,double,double,double,int,int,int>> &vec_xyzdzrgb );

  /// @brief print all std::vector<std::pair<int,int>> to FILE*
  void out_vec_pair( FILE *fout, const std::vector<std::pair<int,int>> &vec_ipair );
  /// @brief print all std::vector<std::pair<int,int>> to ascii file
  void out_vec_pair( const std::filesystem::path& pathout
        , const std::vector<std::pair<int,int>> &vec_ipair );

  /// @brief print all std::vector<std::pair<double,double>> to FILE*
  void out_vec_pair( FILE *fout, const std::vector<std::pair<double,double>> &vec_pair );
  /// @brief print all std::vector<std::pair<double,double>> to ascii file
  void out_vec_pair( const std::filesystem::path& pathout
        , const std::vector<std::pair<double,double>> &vec_pair );

  /// @brief for output std::vector<std::tuple<int,int,int>> to FILE*
  void out_int_triple(FILE *fout, const std::tuple<int,int,int> &tp);
  /// @brief for output std::vector<std::tuple<int,int,int>> to FILE*
  void out_vec_triple( FILE *fout
  , const std::vector<std::tuple<int,int,int>> &vec_triple );
  /// @brief for output std::vector<std::tuple<int,int,int>> to ascii file
  void out_vec_triple( const std::filesystem::path& pathout
  , const std::vector<std::tuple<int,int,int>> &vec_triple );

  /// @brief output (ix,iy,path_length) to FILE*
  void out_ix_iy_PL(FILE *fout, const std::tuple<int,int,double> &tp);
  /// @brief output all (ix,iy,path_length) to FILE*
  void out_vec_ix_iy_PL( FILE *fout
  , const std::vector<std::tuple<int,int,double>> &vec_ix_iy_PL );
  /// @brief output all (ix,iy,path_length) to ascii file
  void out_vec_ix_iy_PL( const std::filesystem::path& pathout
  , const std::vector<std::tuple<int,int,double>> &vec_ix_iy_PL );

  /// @brief for output (ix,iy,iz,path_length) to FILE*
  void out_ix_iy_iz_path(FILE *fout, const std::tuple<int,int,int,double> &tp);
  /// @brief for output all (ix,iy,iz,path_length) to FILE*
  void out_vec_ix_iy_iz_path( FILE *fout
  , const std::vector<std::tuple<int,int,int,double>> &vec_ix_iy_iz_path );
  /// @brief for output all (ix,iy,iz,path_length) to ascii file
  void out_vec_ix_iy_iz_path( const std::filesystem::path& pathout
  , const std::vector<std::tuple<int,int,int,double>> &vec_ix_iy_iz_path );

  /// @brief output all x, y, y1=ax+b, dy=y-y1 to FILE*
  void debug_y_ax_b(FILE *fout, const std::vector<std::pair<double,double>> &vec_xy, const double &a, const double &b );
  /// @brief output all x, y, y1=ax+b, dy=y-y1 to ascii file
  void debug_y_ax_b(const std::filesystem::path& pathout, const std::vector<std::pair<double,double>> &vec_xy, const double &a, const double &b );

  /// @brief Function to output std::set<int> to file with one element per line
  /// @param[in] data Set to output
  /// @param[in] outpath Output file path
  void out_set(
    const std::set<int>& data,
    const std::filesystem::path& pathout);  

  /// @brief for output std::map<int, bool> to file
  void out_map(const std::map<int, bool>& map_int_bool
    , const std::filesystem::path& pathout);

  /// @brief for output std::map<int, int> to file
  void out_map(const std::map<int, int>& map_int_int
    , const std::filesystem::path& pathout);

  // for output std::unordered_multimap<int, int>
  void out_uommp_int2(const std::unordered_multimap<int,int>& uommp
    , const std::filesystem::path& pathout);

  // for output std::unordered_multimap<int, int>
  void out_uomp_int2(const std::unordered_map<int,int>& uomp
    , const std::filesystem::path& pathout);

  /// @brief for output SortedDetidUqigSet to file
  void write_sorted_detid_uqigavail_to_file(
    const SortedDetidUqigSet& data,
    const std::filesystem::path& outpath);

  ///@} ------------------------------------------------------------------

  //=======================================================================
  /// @name file_operation_myapp
  /// @details File operation functions
  /// @{

  /// @brief erace all files in directory. directory is not erased.
  void clear_dir(const std::filesystem::path& dir_path);

  ///@} ------------------------------------------------------------------

  //======================================================================
  /// @name time_myapp
  /// @details Time measurement functions
  /// @{

  /// @brief display time (min) to FILE*
  void cast_time_min( FILE *fout, const std::string &msg
  , const std::chrono::system_clock::time_point &start
  , const std::chrono::system_clock::time_point &end );
  /// @brief display time (sec) to FILE*
  void cast_time_sec( FILE *fout, const std::string &msg
  , const std::chrono::system_clock::time_point &start
  , const std::chrono::system_clock::time_point &end );
  /// @brief display time (msec) to FILE*
  void cast_time_msec( FILE *fout, const std::string &msg
  , const std::chrono::system_clock::time_point &start
  , const std::chrono::system_clock::time_point &end );

  /// @brief display time (min) to spdlog
  void cast_time_min(
      const spdlog::level::level_enum log_level
    , const std::string &msg
    , const std::chrono::system_clock::time_point &start
    , const std::chrono::system_clock::time_point &end );
  /// @brief display time (sec) to spdlog
  void cast_time_sec(
      const spdlog::level::level_enum log_level
    , const std::string &msg
    , const std::chrono::system_clock::time_point &start
    , const std::chrono::system_clock::time_point &end );
  /// @brief display time (msec) to spdlog
  void cast_time_msec(
      const spdlog::level::level_enum log_level
    , const std::string &msg
    , const std::chrono::system_clock::time_point &start
    , const std::chrono::system_clock::time_point &end );
  ///@} ------------------------------------------------------------------


  //===============================================================================
  /// @name linear_interpolation_myapp
  /// @details Interpolation calculations
  /// @{

  // Interpolation calculation
  double linear_interpolation(double x0, double y0, double x1, double y1, double x);
  double bilinear_interpolation(double x1, double x2, double y1, double y2,
                                double f11, double f21, double f12, double f22,
                                double x, double y);

  ///@} ------------------------------------------------------------------


  //=================================================================
  /// @name compare_myapp
  /// @details Comparison functions
  /// @{

  /// @brief Comparison of primitive floating point numbers (specify allowable error with epsilon)
  bool is_eq( const double a, const double b, const double epsilon);

  /// @brief Comparison of std::tuple<double, double, double>
  bool is_eq(
      const std::tuple<double, double, double>& t1
    , const std::tuple<double, double, double>& t2
    , const double epsilon);

  /// @brief Comparison of std::vector< std::tuple<double, double, double> >
  bool is_eq(
      const std::vector< std::tuple<double, double, double> >& v1
    , const std::vector< std::tuple<double, double, double> >& v2
    , const double epsilon=1.0E-9);

  /// @brief compare unordered_multimap<int,int> is same or not
  bool is_eq(const std::unordered_multimap<int, int>& ummap1
           , const std::unordered_multimap<int, int>& ummap2);

  /// @brief compare unordered_multimap<int,int> is same or not \n
  /// Intended to consider stored order
  template <typename Key, typename Value>
  bool is_eq2(const std::unordered_multimap<Key, Value>& ummap1,
              const std::unordered_multimap<Key, Value>& ummap2)
  {
    if (ummap1.size() != ummap2.size()) {
      return false;
    }

    // Map to store frequency of values for each key
    std::map<Key, std::map<Value, int>> frequencyMap1, frequencyMap2;

    for (const auto& pair : ummap1) {
      frequencyMap1[pair.first][pair.second]++;
    }

    for (const auto& pair : ummap2) {
      frequencyMap2[pair.first][pair.second]++;
    }

    // Compare frequency maps
    return frequencyMap1 == frequencyMap2;
  }
  /// @brief Function to determine whether two tuple_int::UmpIntInt3 are equal
  bool are_ump_eq( const tuple_int::UmpIntInt3& map1,const tuple_int::UmpIntInt3& map2 );

  /// @brief Function to determine whether two tuple_int::UmpIntInt3 are equal
  /// @details pointer input version
  // inline bool are_ump_eq( const tuple_int::UmpIntInt3* map1, const tuple_int::UmpIntInt3* map2 ){
  //   return are_ump_eq(*map1, *map2);
  // };

  /// @brief Function to determine whether two std::unordered_map< std::string, int > are equal
  bool are_ump_eq(const tuple_int::UmpInt3Int& map1,const tuple_int::UmpInt3Int& map2 );


  /// @brief for sort. compare tuple<double,double,double> by yxz
  bool compare_tuple_yxz(
    const std::tuple<double, double, double>& a,
    const std::tuple<double, double, double>& b);

  /// @brief for sort. compare array<double,3> by yxz
  bool compare_array_yxz(
    const std::array<double, 3>& a,
    const std::array<double, 3>& b);


  /// @brief Determine whether x of vec_double2 is equally spaced within double_epsilon range
  /// @returns true:equally spaced false:not equally spaced
  bool is_same_interval_x(
      const std::vector<std::pair<double,double>> &vec_double2
    , const double double_epsilon );

  ///@} ------------------------------------------------------------------

  //=================================================================
  /// @name Grouping_myapp
  /// @details Functions for grouping std::vector<int,int> etc.
  /// @{

  /// @brief Split arr array into N parts. Assign remainder to targetGroup-th group. \n
  /// If targetGroup=0, remainder is not assigned.
  std::vector<std::vector<int>> splitArrayIntoN(
    const std::vector<int>& arr, int N, const int targetGroup);

  /// @brief Function to split 2D array into specified grid number
  ///
  /// @param[in] arr 2D array of type `std::vector<std::tuple<int, int>>` to be split
  /// @param[in] nx_div Number of divisions in x direction
  /// @param[in] ny_div Number of divisions in y direction
  /// @param[in] targetGroup_x Group number (1-based) in x direction to assign surplus elements
  /// @param[in] targetGroup_y Group number (1-based) in y direction to assign surplus elements
  /// @return 4D vector. Array split into `ny_div * nx_div` grid.
  ///
  /// @details
  /// This function splits array `arr` into specified `nx_div * ny_div` grid. \n
  /// Each grid holds a partial subset of the 2D array, with each element storing a `std::tuple<int, int>` pair. \n
  /// If the array is not divisible by `nx_div` or `ny_div`, the surplus is assigned to the group specified by `targetGroup_x`, `targetGroup_y`. \n
  /// Elements within each split group preserve the order of the original array.
  ///
  /// For example, if the array is 6x8 in size and `nx_div=3`, `ny_div=2`, 6 groups are created after splitting, and any surplus \n
  /// is assigned to the specified target group.
  ///
  /// @note
  /// `targetGroup_x` and `targetGroup_y` are 1-based indices.
  ///
  /// @attention
  /// If the size of array `arr` is not divisible by `nx_div` or `ny_div`, the remainder is assigned to the group specified by `targetGroup_x`, `targetGroup_y`, \n
  /// so each group may have different sizes.
  std::vector<std::vector<std::vector<std::vector<std::array<int,2>>>>>
    split2DArrayIntoGrid(
      const std::vector<std::vector<std::array<int,2>>>& arr
    , const int nx_div, const int ny_div
    , const int targetGroup_x, const int targetGroup_y);



  void printResult(
    const std::vector<std::vector<std::vector<std::vector<std::array<int,2>>>>>& result);

  void printGroups(
    const std::vector<std::vector<std::vector<std::array<int,2>>>>& groups);

  // Function to copy int corresponding to specified num from std::vector<std::array<int,3>>
  // to std::vector<int>
  std::vector<int> copyArrayElementToVector(
    const std::vector<std::array<int,3>>& input, const int num );

  ///@} ------------------------------------------------------------------


  //----------------------------------------------------------------------------------



  //----------------------------------------------------------------------------------
  
  // ======================================================================
  /// @name json_myapp
  /// @details JSON-related functions
  /// @{

  /// @brief Function to load json file ignoring comments
  /// @param pathin
  /// @return nlohmann::json
  nlohmann::json load_json(const std::filesystem::path &pathin);

  /// @brief Function to format JSON and output to spdlog
  /// @note usage1 : dump_with_section(js, {}, "address", spdlog::level::debug); --> output hierarchy up to "address" and "address" section \n
  /// usage2 : dump_with_section(js, {"profile", "address"}, "location", spdlog::level::debug); --> output only "location" section \n
  /// usage3 : dump_with_section(js, {}, "", spdlog::level::debug); --> output all hierarchy
  void dump_json(
      const nlohmann::json& js
    , const std::vector<std::string>& current_path
    , const std::string& target_section
    , const spdlog::level::level_enum& log_level
    // Do not specify arguments below normally
    , const int indent = 2
    , std::ostringstream* oss_ptr = nullptr
    , bool already_wrapped = false
    , bool inline_mode = false );

  /// @brief Function to format JSON and output to FILE *fout
  /// @note usage is same as above
  void dump_json(
      const nlohmann::json& js
    , const std::vector<std::string>& current_path
    , const std::string& target_section
    , FILE* fout
    // Do not specify arguments below normally
    , const int indent = 2
    , bool already_wrapped = false
    , bool inline_mode = false );


  ///@} ------------------------------------------------------------------
};

//#####################################################################################
/// @namespace mysort
/// @brief Sorting functions
/// @ingroup basicTools
//#####################################################################################
namespace mysort {

  ///@brief Sort std::vector<std::pair<double,double>> in ascending order by first
  void sort_vec_double2_x_increasing(
    std::vector<std::pair<double,double>> &vec_double2 );

  ///@brief Sort std::vector<std::pair<double,double>> in ascending order by second
  void sort_vec_double2_y_increasing(
    std::vector<std::pair<double,double>> &vec_double2 );

  // @brief Sort std::vector<std::tuple<x,y,z>> in ascending order by x
  // void sort_vec_double3_x_increasing(
  //   std::vector<std::array<double,3>> &vec_double3 );

  /// @brief sort Eigen_Vector3d in descending order with xy
  /// @details (x1,y1), (x2,y1), (x3,y1), .... , (x1,y2), (x2,y2), (x3,y2) .... \n
  // for sorting std::vector<Eigen::Vector3d> \n
  // usage : std::sort(vec_v3.begin(),vec_v3.end(),Eigen_Vector3d_sort_xy());
  struct Eigen_Vector3d_sort_xy{
    bool operator()(const Eigen::Vector3d& a, const Eigen::Vector3d& b){
      if (a.x() == b.x()){
        return a.y() < b.y();
      }
      else{
        return a.x() < b.x();
      }
    }
  };
  /// @brief for sorting std::vector<Eigen::Vector3d> in descending order yx
  /// @details (x1,y1), (x1,y2), (x1,y3), .... , (x2,y1), (x2,y2), (x2,y3) , ... \n
  /// usage : std::sort(vec_v3.begin(),vec_v3.end(),Eigen_Vector3d_sort_yx());
  struct Eigen_Vector3d_sort_yx{
    bool operator()(const Eigen::Vector3d& a, const Eigen::Vector3d& b)
    {
      if (a.y() == b.y()){
        return a.x() < b.x();
      }else{
        return a.y() < b.y();
      }
    }
  };

  /// @brief for sorting std::vector<std::pair<double,double>> in descending order with xy
  struct vec_double2_sort_x_increasing{
    bool operator()(const std::pair<double,double>& a, const std::pair<double,double>& b)
    {
      if (a.first == b.first) {
        return a.second < b.second;
      } else {
        return a.first < b.first;
      }
    }

  };
  /// @brief for sorting std::vector<std::pair<double,double>> in descending order with yx
  struct vec_double2_sort_y_increasing {
    bool operator()(const std::pair<double, double>& a, const std::pair<double, double>& b) {
      if (a.second == b.second) {
        return a.first < b.first;
      } else {
        return a.second < b.second;
      }
    }
  };
};


// struct vec_double3_sort_x_increasing{
//   bool operator()(const std::array<double,3>& a
//                 , const std::array<double,3>& b)
//   {
//     const auto&[x1,y1,z1] = a;
//     const auto&[x2,y2,z2] = b;
//     if (x1 == x2) {
//       return y1 < y2;
//     } else {
//       return x1 < x2;
//     }
//   }
// };

//===================================================

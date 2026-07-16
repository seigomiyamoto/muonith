/// @file ns_tuple_int.hpp
/// @brief Integer tuple type definitions and utilities for map-based indexing
/// @details
/// This module provides utilities for working with std::tuple<int,int> and std::tuple<int,int,int>
/// as keys in unordered maps. It includes:
/// - Custom hash functions (Tuple2Hash, Tuple3Hash) for tuple-based map keys
/// - Type aliases for common map patterns (UmpIntInt2, UmpInt2Int, UmpIntInt3, UmpInt3Int)
/// - Stream output operators for convenient tuple printing
/// - Conversion functions between tuples and strings
/// - Map lookup utilities with sentinel values for not-found cases
///
/// Typical usage:
/// @code
/// // Create a map from int2-tuple to int
/// tuple_int::UmpInt2Int edge_map;
/// edge_map[std::make_tuple(0,1)] = 100;
///
/// // Lookup with sentinel value
/// int result = tuple_int::get_value_from_key(edge_map, 0, 1); // returns 100
/// int not_found = tuple_int::get_value_from_key(edge_map, 5, 6); // returns -1
/// @endcode
///
/// @note Thread-safety: All functions are thread-safe for read-only operations on shared maps.
///       Concurrent modifications require external synchronization.
#pragma once

#include <string>
#include <cstdio>
#include <iostream>
#include <fstream>
#include <set>

#include <chrono>
#include <thread>

#define _USE_MATH_DEFINES // For M_PI constant, must be defined before <cmath>
#include <cmath>

// std::map<>
#include <map>
#include <unordered_map>
#include <algorithm>

#include <cassert>
#include <vector>
// for std::filesystem::path
#include <filesystem>

//#########################################################################
//#########################################################################
/// @namespace tuple_int
/// @brief Functions for handling std::tuple<int,int>, std::tuple<int,int,int>, and their associated maps
//#########################################################################
//#########################################################################
namespace tuple_int {
  //==================================================================
  /// @name type definitions and related function of std::unordered_map
  ///@{
  
  /// @brief Hash function for std::tuple<int,int>
  /// @details Uses XOR with bit-shift to combine hash values of tuple elements.
  ///          Provides reasonable hash distribution for unordered_map usage.
  /// @note Complexity: O(1)
  struct Tuple2Hash {
    /// @param tuple Tuple to hash
    /// @return Hash value
    size_t operator()(const std::tuple<int,int>& tuple) const {
      const auto& [a, b] = tuple;
      size_t hash_a = std::hash<int>()(a);
      size_t hash_b = std::hash<int>()(b);
      
      // Combine hash values
      return hash_a ^ (hash_b << 1);
    }
  };

  /// @brief Equality comparison for std::tuple<int,int>
  /// @details Delegates to tuple's built-in operator==
  /// @note Complexity: O(1)
  struct Tuple2Equal {
    /// @param t1 First tuple
    /// @param t2 Second tuple
    /// @return true if tuples are equal, false otherwise
    bool operator()(const std::tuple<int,int>& t1, const std::tuple<int,int>& t2) const {
      return t1 == t2;
    }
  };

  using UmpIntInt2 = std::unordered_map<int, std::tuple<int, int>>;
  using UmpInt2Int = std::unordered_map<std::tuple<int,int>, int, Tuple2Hash, Tuple2Equal>;

  /// @brief Hash function for std::tuple<int,int,int>
  /// @details Uses XOR with bit-shift to combine hash values of tuple elements.
  ///          Provides reasonable hash distribution for unordered_map usage.
  /// @note Complexity: O(1)
  struct Tuple3Hash {
    /// @param tuple Tuple to hash
    /// @return Hash value
    size_t operator()(const std::tuple<int, int, int>& tuple) const {
      const auto& [a, b, c] = tuple;
      size_t hash_a = std::hash<int>()(a);
      size_t hash_b = std::hash<int>()(b);
      size_t hash_c = std::hash<int>()(c);
      
      // Combine hash values
      return hash_a ^ (hash_b << 1) ^ (hash_c << 2);
    }
  };

  /// @brief Equality comparison for std::tuple<int,int,int>
  /// @details Delegates to tuple's built-in operator==
  /// @note Complexity: O(1)
  struct Tuple3Equal {
    /// @param t1 First tuple
    /// @param t2 Second tuple
    /// @return true if tuples are equal, false otherwise
    bool operator()(const std::tuple<int, int, int>& t1, const std::tuple<int, int, int>& t2) const {
      return t1 == t2;
    }
  };

  using UmpIntInt3 = std::unordered_map<int, std::tuple<int,int,int>>;
  using UmpInt3Int = std::unordered_map<std::tuple<int,int,int>, int, Tuple3Hash, Tuple3Equal>;


  ///@} ------------------------------------------------------------------

  //=======================================================================
  /// @name operators for std::tuple<int,int>, std::tuple<int,int,int>
  ///@{
  
  /// @brief Enables direct stream output for std::tuple<int,int> \n
  /// std::cout << "tuple=" << t << std::endl;
  inline std::ostream& operator<<(std::ostream& os, const std::tuple<int, int>& t) {
    os << "(" << std::get<0>(t) << "," << std::get<1>(t) << ")";
    return os;
  };

  /// @brief Enables direct stream output for std::tuple<int,int,int> \n
  /// std::cout << "tuple=" << t << std::endl;
  inline std::ostream& operator<<(std::ostream& os, const std::tuple<int,int,int>& t) {
    os << "(" << std::get<0>(t) << "," << std::get<1>(t) << "," << std::get<2>(t) << ")";
    return os;
  };
  ///@} ------------------------------------------------------------------

  //=======================================================================
  /// @name constants
  ///@{

  /// @brief const tuple<int,int> for not found
  static constexpr std::tuple<int,int> tp_int2_notfound = std::make_tuple(-1,-1);
  
  /// @brief const tuple<int,int,int> for not found
  static constexpr std::tuple<int,int,int> tp_int3_notfound = std::make_tuple(-1,-1,-1);
  
  /// @brief const value for not found
  static constexpr int int_notfound = -1;

  /// @brief char delimiter when reading from ascii file
  static constexpr char char_delim = ' ';
  ///@} ------------------------------------------------------------------

  //=======================================================================
  /// @name std::tuple<int,int>, UmpIntInt2, UmpInt2Int functions
  ///@{


  /// @brief Outputs a map of int to int2-tuple to a file
  /// @param[in] fout Output file pointer (must be valid and open for writing)
  /// @param[in] ump_int_int2 Map to output
  /// @note Format: Each line contains "key int0 int1"
  void out_ump_int_int2(FILE* fout, const UmpIntInt2& ump_int_int2);

  /// @brief Outputs a map of int2-tuple to int to a file
  /// @param[in] fout Output file pointer (must be valid and open for writing)
  /// @param[in] ump_int2_int Map to output
  /// @note Format: Each line contains "i0 i1 value"
  void out_ump_int2_int(FILE* fout, const UmpInt2Int& ump_int2_int);

  /// @brief Outputs a map of int to int2-tuple to a file
  /// @param[in] pathout Output file path
  /// @param[in] ump_int_int2 Map to output
  /// @note Format: Each line contains "key int0 int1"
  /// @throws std::runtime_error If file cannot be opened for writing
  void out_ump_int_int2( const std::filesystem::path& pathout, const UmpIntInt2 &ump_int_int2 );

  /// @brief Outputs a map of int2-tuple to int to a file
  /// @param[in] pathout Output file path
  /// @param[in] ump_int2_int Map to output
  /// @note Format: Each line contains "i0 i1 value"
  /// @throws std::runtime_error If file cannot be opened for writing
  void out_ump_int2_int( const std::filesystem::path& pathout, const UmpInt2Int &ump_int2_int );

  /// @brief Retrieves the key corresponding to an int2-tuple value from a map
  /// @param[in] ump Map to search
  /// @param[in] value Value (tuple) to search for
  /// @return Corresponding key
  /// @throws std::runtime_error If key is not found
  /// @note Complexity: O(n) where n is map size (performs linear search)
  int get_key_from_value( const UmpIntInt2 &ump, const std::tuple<int,int> &value );

  /// @brief Retrieves the int2-tuple value corresponding to a key from a map
  /// @param[in] ump Map to search
  /// @param[in] key Key to search for
  /// @return Corresponding value (tuple), or tp_int2_notfound (-1,-1) if not found
  /// @note Complexity: O(1) average case
  std::tuple<int,int> get_value_from_key(const UmpIntInt2 &ump, const int &key );

  /// @brief Retrieves the int value corresponding to a key from a map
  /// @param[in] ump Map to search
  /// @param[in] key Key to search for (std::tuple<int,int>)
  /// @return Corresponding value (int), or int_notfound (-1) if not found
  /// @note Complexity: O(1) average case
  int get_value_from_key( const UmpInt2Int &ump, const std::tuple<int,int> &key );

  /// @brief Retrieves the int value corresponding to a pair of ints from a map
  /// @param[in] ump Map to search
  /// @param[in] i0 First int of the key pair
  /// @param[in] i1 Second int of the key pair
  /// @return Corresponding value (int), or int_notfound (-1) if not found
  /// @note Complexity: O(1) average case
  int get_value_from_key( const UmpInt2Int &ump, const int i0, const int i1 );

  /// @brief Checks if two int2-tuples are equal
  /// @param tp0 First tuple
  /// @param tp1 Second tuple
  /// @return true if the two tuples are equal, false otherwise
  bool equal(const std::tuple<int,int> &tp0, const std::tuple<int,int> &tp1);

  /// @brief Converts an int2-tuple to a string
  /// @param tp Tuple to convert
  /// @return Converted string
  std::string to_string( const std::tuple<int,int> &tp );

  /// @brief Converts two int values to a string
  /// @param i0 First int value
  /// @param i1 Second int value
  /// @return Converted string
  std::string to_string( const int i0, const int i1 );

  /// @brief Converts a string to an int2-tuple
  /// @param str String to convert (format: "i0 i1" with space delimiter)
  /// @return Converted int2-tuple
  /// @throws std::runtime_error If string does not contain exactly 2 integers
  /// @throws std::invalid_argument If string contains non-integer values
  std::tuple<int,int> string_to_int2( const std::string &str );

  ///@} ------------------------------------------------------------------


  //=======================================================================
  /// @name std::tuple<int,int,int>, UmpIntInt3, UmpInt3Int functions
  ///@{


  /// @brief Outputs a map of int to int3-tuple to a file
  /// @param[in] fout Output file pointer (must be valid and open for writing)
  /// @param[in] ump_int_int3 Map to output
  /// @note Format: Each line contains "key i0 i1 i2"
  void out_ump_int_int3(FILE* fout
  , const std::unordered_map<int, std::tuple<int,int,int>>& ump_int_int3);

  /// @brief Outputs a map of int3-tuple to int to a file
  /// @param[in] fout Output file pointer (must be valid and open for writing)
  /// @param[in] ump_int3_int Map to output
  /// @note Format: Each line contains "i0 i1 i2 value"
  void out_ump_int3_int(FILE* fout, const UmpInt3Int& ump_int3_int);

  /// @brief Outputs a map of int to int3-tuple to a file
  /// @param[in] pathout Output file path
  /// @param[in] ump_int_int3 Map to output
  /// @note Format: Each line contains "key i0 i1 i2"
  /// @throws std::runtime_error If file cannot be opened for writing
  void out_ump_int_int3( const std::filesystem::path& pathout
  , const UmpIntInt3 &ump_int_int3 );

  /// @brief Outputs a map of int3-tuple to int to a file
  /// @param[in] pathout Output file path
  /// @param[in] ump_int3_int Map to output
  /// @note Format: Each line contains "i0 i1 i2 value"
  /// @throws std::runtime_error If file cannot be opened for writing
  void out_ump_int3_int( const std::filesystem::path& pathout
  , const UmpInt3Int &ump_int3_int );

  /// @brief Retrieves the key corresponding to an int3-tuple value from a map
  /// @param[in] ump Map to search
  /// @param[in] value Value (tuple) to search for
  /// @return Corresponding key
  /// @throws std::runtime_error If key is not found
  /// @note Complexity: O(n) where n is map size (performs linear search)
  int get_key_from_value( const UmpIntInt3 &ump, const std::tuple<int,int,int> &value);

  /// @brief Retrieves the int3-tuple value corresponding to a key from a map
  /// @param[in] ump Map to search
  /// @param[in] key Key to search for
  /// @return Corresponding value (tuple), or tp_int3_notfound (-1,-1,-1) if not found
  /// @note Complexity: O(1) average case
  std::tuple<int,int,int> get_value_from_key( const UmpIntInt3 &ump, const int &key );

  /// @brief Retrieves the int value corresponding to a string key from a map
  /// @param[in] ump Map to search
  /// @param[in] key Key to search for (string format: "i0 i1 i2")
  /// @return Corresponding value (int), or int_notfound (-1) if not found
  /// @note Complexity: O(1) average case for map lookup
  /// @note This function is declared but may not be implemented
  int get_value_from_key( const UmpInt3Int &ump, const std::string &key );

  /// @brief Retrieves the int value corresponding to a triple of ints from a map
  /// @param[in] ump Map to search
  /// @param[in] i0 First int of the key triple
  /// @param[in] i1 Second int of the key triple
  /// @param[in] i2 Third int of the key triple
  /// @return Corresponding value (int), or int_notfound (-1) if not found
  /// @note Complexity: O(1) average case
  int get_value_from_key( const UmpInt3Int &ump, const int i0, const int i1, const int i2 );

  /// @brief Retrieves the int value corresponding to an int3-tuple key from a map
  /// @param[in] ump Map to search
  /// @param[in] key Key to search for (int3-tuple)
  /// @return Corresponding value (int), or int_notfound (-1) if not found
  /// @note Complexity: O(1) average case
  int get_value_from_key( const UmpInt3Int &ump, const std::tuple<int,int,int> &key );

  /// @brief Checks if two int3-tuples are equal
  /// @param tp0 First tuple
  /// @param tp1 Second tuple
  /// @return true if the two tuples are equal, false otherwise
  bool equal( const std::tuple<int,int,int> &tp0, const std::tuple<int,int,int> &tp1);

  /// @brief Converts an int3-tuple to a string
  /// @param tp Tuple to convert
  /// @return Converted string
  std::string to_string( const std::tuple<int,int,int> &tp );

  /// @brief Converts three int values to a string
  /// @param i0 First int value
  /// @param i1 Second int value
  /// @param i2 Third int value
  /// @return Converted string
  std::string to_string( const int i0, const int i1, const int i2 );

  /// @brief Converts a string to an int3-tuple
  /// @param str String to convert (format: "i0 i1 i2" with space delimiter)
  /// @return Converted int3-tuple
  /// @throws std::runtime_error If string does not contain exactly 3 integers
  /// @throws std::invalid_argument If string contains non-integer values
  std::tuple<int,int,int> string_to_int3( const std::string &str );

  /// @brief Extracts the set of keys from a std::map<int,int>
  std::set<int> get_set_key( const std::map<int,int> &map );

  /// @brief Extracts the set of values from a std::map<int,int>
  std::set<int> get_set_value( const std::map<int,int> &map );

  /// @brief Extracts the vector of keys from a std::map<int,int>
  std::vector<int> get_vec_key( const std::map<int,int> &map );

  /// @brief Extracts the vector of values from a std::map<int,int>
  std::vector<int> get_vec_value( const std::map<int,int> &map );

  ///@} ------------------------------------------------------------------

}

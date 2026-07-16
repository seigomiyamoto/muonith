/// @file cls_Grid1d.hpp
/// @brief One-dimensional grid class for spatial discretization
/// @details
/// Provides a base class for 1D grid structures used in ray tracing and spatial calculations.
/// The Grid1d class represents a uniform grid with configurable bin count, bounds, and interval.
///
/// **Key Features:**
/// - Uniform bin spacing across the grid range
/// - Index-to-value and value-to-index conversion
/// - Range queries and bin extraction
/// - Grid merging (coarsening) and splitting (refinement)
/// - Binary serialization support
///
/// **Coordinate System:**
/// - Bins are defined by [min, max) with uniform interval
/// - Index range: [0, nbin)
/// - Closed interval semantics for range operations
///
/// **Thread-Safety:**
/// - Const methods are thread-safe for read-only operations
/// - Non-const methods require external synchronization
///
/// **Units:**
/// - All spatial values are in the same units as the problem domain (typically meters)
/// - Intervals and tolerances must match these units
#pragma once

#include <Eigen/Dense>
#include <limits>
#include <source_location>

#include <fstream>
#include <iostream>
#include <sstream>
#include <string>
#include <cstdio>
#include <cmath>
#include <functional>
#include <algorithm>
#include <vector>
#include "spdlog_pch.hpp"
#include "ns_mylogger.hpp"

using srcloc = std::source_location;

//####################################################################
//####################################################################
/// @class Grid1d
/// @brief Class representing a one-dimensional grid
/// @ingroup basicGridClasses
///
/// @details
/// **Responsibilities:**
/// - Maintain a uniform 1D grid with nbin bins spanning [min, max)
/// - Provide bidirectional mapping between indices and spatial values
/// - Support grid transformation operations (merge, split, cut)
/// - Validate grid consistency and range queries
///
/// **Invariants:**
/// - nbin >= 1
/// - min < max
/// - interval > 0.0
/// - interval == (max - min) / nbin (within tolerance)
/// - std::remainder(max - min, interval) ~= 0 (within tolerance)
///
/// @note
/// min/max always hold the canonical bin-edge coordinates.  Any half-interval
/// shift (bin-center interpretation of raw input data) is applied exactly once
/// at grid construction time (e.g. Grid2d::set_xy_axis_from_vec_xyz); save()
/// and load() are pass-through and never transform the stored values.
///
/// **Typical Usage:**
/// @code
/// // Create a grid with 100 bins from 0.0 to 10.0
/// Grid1d grid("example", 100, 0.0, 10.0, 0.1);
///
/// // Convert value to index
/// int idx = grid.get_index(5.5);
///
/// // Get bin boundaries
/// double lower = grid.get_lower_value(idx);
/// double upper = grid.get_upper_value(idx);
///
/// // Merge bins by factor of 2
/// int merge_factor = 2;
/// Grid1d coarse = grid.get_merged(5.0, merge_factor);
/// @endcode
//####################################################################
//####################################################################
class Grid1d {
  private:
    /// @brief name of this instance
    std::string name = "none";

    /// @brief number of bins
    int nbin = 999;

    /// @brief max, min, interval
    double min = 0.0;
    double max = 999.0;
    double interval = 1.0;
  public:
    /// @brief Structure representing an index range (closed interval)
    struct RangeIndices {
      int start = 0;              ///< Start bin after extraction (closed interval)
      int end = -1;               ///< End bin after extraction (closed interval)
    };

    //==================================================================
    /// @name constant variables
    ///@{
    
    /// @brief lower range error codes
    static constexpr int OUT_OF_RANGE_LOWER = -std::numeric_limits<int>::max();
    
    /// @brief upper range error codes
    static constexpr int OUT_OF_RANGE_UPPER = std::numeric_limits<int>::max();
    ///@} ------------------------------------------------------------------

    //==================================================================
    /// @name constructor_destructor
    ///@{
    
    /// @brief Default constructor
    Grid1d() = default;

    /// @brief Copy constructor
    Grid1d(const Grid1d &org) = default;

    /// @brief Move constructor
    Grid1d(Grid1d&& other) noexcept = default;

    /// @brief Fully-parameterized constructor
    /// @param[in] name_in Name identifier for this grid instance
    /// @param[in] nbin_in Number of bins (must be >= 1)
    /// @param[in] min_in Lower bound of the grid (exclusive at upper boundary)
    /// @param[in] max_in Upper bound of the grid (exclusive)
    /// @param[in] interval_in Bin width (must be > 0.0 and equal to (max_in - min_in) / nbin_in within tolerance)
    /// @param[in] tolerance_ratio Relative tolerance for consistency checks (default: 1.0e-4)
    /// @throws std::runtime_error If grid parameters violate invariants
    Grid1d(const std::string name_in,
           const int nbin_in,
           const double min_in, const double max_in,
           const double interval_in,
           const double tolerance_ratio = 1.0e-4);

    /// @brief destructor
    virtual ~Grid1d() = default;
    ///@} ------------------------------------------------------------------

    //==================================================================
    /// @name operators
    ///@{
    
    /// @brief Assignment operator
    Grid1d& operator=(const Grid1d& other) = default;

    /// @brief Inequality operator
    /// @note Does not compare names.
    bool operator!=(const Grid1d& other) const;

    /// @brief Equality operator (defined using inequality operator)
    bool operator==(const Grid1d& other) const {
      return !(*this != other);
    }
    ///@} ------------------------------------------------------------------
    
    //==================================================================
    /// @name check functions
    ///@{
    
    /// @brief Check if an index is within valid range
    /// @param[in] index Index to check
    /// @return true if index is in [0, nbin), false otherwise
    bool is_index_inside( const int index ) const;

    /// @brief Check if a value is within grid bounds
    /// @param[in] value Value to check
    /// @return true if value is in [min, max), false otherwise
    bool is_value_inside( const double value ) const;

    /// @brief check index and if it is outside, throw error
    /// @param index Index value
    /// @param loc Source location information (automatically captured)
    void check_inside( const int index,
      srcloc loc = srcloc::current() ) const;

    /// @brief check value and if it is outside, throw error
    /// @param value Value to check
    /// @param loc Source location information (automatically captured)
    void check_inside( const double value,
      srcloc loc = srcloc::current() ) const;

    /// @brief Validate grid consistency and invariants
    /// @param[in] tolerance_ratio Relative tolerance for interval consistency checks (default: 1.0e-4)
    /// @return true if all invariants hold, false otherwise
    /// @note Checks: nbin >= 1, min < max, interval > 0, interval == (max-min)/nbin within tolerance
    bool is_check_variables(
      const double tolerance_ratio = 1.0e-4 ) const;

    /// @brief return true if index is above lower limit. if not, return false. 
    bool is_above_lower_limit(int index_in) const { return index_in >= 0;   };
    
    /// @brief return true if index is below upper limit. if not, return false.
    bool is_below_upper_limit(int index_in) const { return index_in <  nbin; };

    ///@} ------------------------------------------------------------------

    //==================================================================
    /// @name setter_functions
    ///@{

    /// @brief Assign new grid parameters
    /// @param[in] name_in Name identifier for this grid instance
    /// @param[in] nbin_in Number of bins (must be >= 1)
    /// @param[in] min_in Lower bound of the grid
    /// @param[in] max_in Upper bound of the grid
    /// @param[in] interval_in Bin width (must be > 0.0)
    /// @param[in] tolerance_ratio Relative tolerance for consistency checks (default: 1.0e-4)
    /// @throws std::runtime_error If grid parameters violate invariants
    void assign(const std::string name_in
          , const int nbin_in
          , const double min_in, const double max_in
          , const double interval_in
          , const double tolerance_ratio = 1.0e-4 );
    
    /// @brief set name of this instance
    void set_name( const std::string &name_in ){ name = name_in; };

    /// @brief set values from class Grid1d
    void set_value( const Grid1d &org ){
      name = org.name;
      nbin = org.nbin;
      min = org.min;
      max = org.max;
      interval = org.interval;
    };
    ///@} ------------------------------------------------------------------

    //==================================================================
    ///@name getter_functions
    ///@{
    
    /// @brief get name of this instance
    std::string get_name() const { return name; };

    /// @brief get nbin
    int get_nbin() const { return nbin; };

    /// @brief get min
    double get_min() const { return min; };

    /// @brief get max
    double get_max() const { return max; };

    /// @brief get interval
    double get_interval() const { return interval; };

    /// @brief Convert a spatial value to its bin index
    /// @param[in] value Spatial coordinate value
    /// @return Bin index in [0, nbin), or OUT_OF_RANGE_LOWER if value < min, or OUT_OF_RANGE_UPPER if value >= max
    /// @note Uses small epsilon tolerance to handle floating-point rounding near boundaries
    int get_index( const double value ) const;

    /// @brief Get the lower boundary of a bin by index
    /// @param[in] index Bin index
    /// @return Lower boundary value = min + index * interval
    /// @throws std::runtime_error If index is out of range
    double get_lower_value( const int index ) const;

    /// @brief Get the lower boundary of the bin containing a value
    /// @param[in] value Spatial coordinate value
    /// @return Lower boundary value of the bin containing value
    /// @throws std::runtime_error If value is out of range
    double get_lower_value( const double value ) const;

    /// @brief Get the upper boundary of a bin by index
    /// @param[in] index Bin index
    /// @return Upper boundary value = min + (index + 1) * interval
    /// @throws std::runtime_error If index is out of range
    double get_upper_value( const int index ) const;

    /// @brief Get the upper boundary of the bin containing a value
    /// @param[in] value Spatial coordinate value
    /// @return Upper boundary value of the bin containing value
    /// @throws std::runtime_error If value is out of range
    double get_upper_value( const double value ) const;

    /// @brief Calculate the index range of bins that fall within the specified range
    /// @param[in] lower Lower bound of the query range
    /// @param[in] upper Upper bound of the query range
    /// @param[in] eps Tolerance ratio relative to interval (default: 1.0e-6)
    /// @return RangeIndices struct with start and end indices (closed interval)
    /// @throws std::runtime_error If upper < lower, eps < 0, or nbin <= 0
    /// @note Applies conservative tolerance by shrinking the range (lower + tol, upper - tol)
    RangeIndices calc_range_indices(const double lower, const double upper,
                                    const double eps = 1.0e-6) const;

    /// @brief Get the center value of a bin by index
    /// @param[in] index Bin index
    /// @return Center value = min + (index + 0.5) * interval
    /// @throws std::runtime_error If index is out of range
    double get_center_value( const int index ) const;

    /// @brief Get the center value of the bin containing a value
    /// @param[in] value Spatial coordinate value
    /// @return Center value of the bin containing value
    /// @throws std::runtime_error If value is out of range
    double get_center_value( const double value ) const;

    /// @brief Find the nearest bin boundary to a given value
    /// @param[in] value_in Input spatial coordinate value
    /// @param[in] tf_cnt If true, shift value_in by -interval*0.5 before finding nearest boundary
    /// @return The bin boundary (lower or upper) closest to the (possibly shifted) value
    /// @throws std::runtime_error If the (possibly shifted) value is out of range
    double get_nearest_value( const double value_in, const bool tf_cnt ) const;

    /// @brief Create a merged (coarsened) Grid1d by combining adjacent bins
    /// @param[in] center Center point that must lie on a bin boundary of the merged grid
    /// @param[in,out] merge_factor Number of bins to merge (modified if too large for the grid)
    /// @return New Grid1d instance with interval = this->interval * merge_factor
    /// @throws std::runtime_error If center is out of range or merge_factor < 1
    /// @details
    /// Merging proceeds symmetrically in both directions from the center.
    /// Bins that cannot be merged at either end are discarded.
    /// If merge_factor is too large given icnt and nbin, it will be automatically adjusted
    /// to the maximum feasible value in the longer direction.
    Grid1d get_merged( const double center, int& merge_factor ) const;

    /// @brief Extract a sub-grid containing only bins within the specified range
    /// @param[in] lower Lower bound value
    /// @param[in] upper Upper bound value
    /// @param[in] eps Tolerance ratio relative to interval (default: 1.0e-6)
    /// @return New Grid1d instance spanning the extracted range
    /// @note If the range would extend beyond the grid, the result uses the original min/max as boundaries.
    /// The interval is preserved from the original grid.
    Grid1d cut( const double lower, const double upper, const double eps = 1.0E-6 ) const;

    /// @brief Map a merged bin index back to the original grid's index range
    /// @param[in] g1_org Original (fine) grid before merging
    /// @param[in] g1_merged Merged (coarse) grid
    /// @param[in] index_merged Bin index in the merged grid
    /// @return Array [index_min, index_max] representing the range of original bins covered by the merged bin
    /// @throws std::runtime_error If index_merged < 0
    std::array<int,2> get_original_index_min_max(
        const Grid1d &g1_org, const Grid1d &g1_merged, const int index_merged ) const;

    /// @brief Create a refined (split) Grid1d by subdividing each bin
    /// @param[in] split_factor Number of subdivisions per bin (must be >= 1)
    /// @return New Grid1d instance with interval = this->interval / split_factor and nbin = this->nbin * split_factor
    /// @throws std::runtime_error If split_factor < 1
    Grid1d get_split( const int split_factor ) const;

    ///@} ------------------------------------------------------------------

    //==================================================================
    /// @name output_functions
    ///@{

    /// @brief Output grid information to a file stream
    /// @param[in] fout Output file stream (default: stderr)
    void out_info( FILE *fout = stderr ) const;

    /// @brief Output grid information using spdlog
    /// @param[in] log_level Logging level for the output
    void out_info(spdlog::level::level_enum log_level) const;
    ///@} ------------------------------------------------------------------


    //==================================================================
    /// @name binary_io_functions
    ///@{
    
    /// @brief Serialize grid to binary stream (pass-through).
    /// @details Writes the canonical values as-is.  Binary layout:
    ///          name + nbin + min + max + interval.
    /// @note Symmetric pair with load(): the exact stored bytes are restored,
    ///       so a save/load roundtrip is bit-identical.
    /// @param[in,out] ofs Output file stream
    /// @throws std::runtime_error If stream fails during write
    void save( std::ofstream& ofs ) const;

    /// @brief Deserialize grid from binary stream (pass-through).
    /// @details Reads name, nbin, min, max, interval from the stream exactly
    ///          as written by save(); no shift transform is applied.
    /// @note Symmetric pair with save().
    /// @param[in,out] ifs Input file stream
    /// @param[in] tolerance_ratio Tolerance for Grid1d consistency checks
    /// @throws std::runtime_error If stream fails during read or validation fails
    void load( std::ifstream& ifs, double tolerance_ratio = 1.0e-4 );
    ///@} ------------------------------------------------------------------
};

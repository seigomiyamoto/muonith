/// @file cls_Grid1dXZ.hpp
/// @brief One-dimensional grid in the XZ plane
/// @details
/// Specialized 1D grid class that associates scalar z-values with uniformly-spaced x-axis positions.
/// Derived from Grid1d, this class provides:
/// - Storage and retrieval of z-values indexed by x-position
/// - Linear interpolation in both directions (x→z and z→x)
/// - Differential computation
/// - Data import/export from/to ASCII files
///
/// Thread-safety: Not thread-safe. External synchronization required for concurrent access.
#pragma once

#include <cstdio>
#include <map>
#include <fstream>
#include <iostream>
#include <sstream> // istringstream
#include <string>
#include <cmath>
#include <functional>  //for sorting
#include <algorithm>//for sorting
#include <vector>
#include <filesystem> //for std::filesystem::path

#include <Eigen/Dense>
#include "cls_Grid1d.hpp"

//###################################################################################
//###################################################################################
/// @class Grid1dXZ
/// @brief One-dimensional grid holding scalar z-values on a uniformly-spaced x-axis
/// @details
/// Grid1dXZ extends Grid1d to store scalar z-values at each x-axis position.
/// The x-axis is uniformly spaced and managed by the base Grid1d class.
///
/// **Responsibilities:**
/// - Maintain z-value array aligned with x-axis grid points
/// - Provide bidirectional interpolation (x→z and z→x)
/// - Support data import from ASCII files with validation
/// - Export grid data to files
///
/// **Invariants:**
/// - vec_z.size() == Grid1d::get_nbin()
/// - zmin <= all z-values <= zmax (updated on modification)
///
/// **Typical usage:**
/// @code
/// // Load from file with shift applied to grid bounds
/// Grid1dXZ grid("data.txt", true);
///
/// // Query z-value at specific x-position
/// double z = grid.get_linear_interpolated_z(1.5);
///
/// // Find x-positions where z == target_value (may return multiple solutions)
/// std::vector<double> x_vals = grid.get_linear_interpolated_x(target_value);
/// @endcode
///
/// @note Thread-safety: No. Concurrent access requires external synchronization.
/// @ingroup basicGridClasses
//###################################################################################
//###################################################################################
class Grid1dXZ : public Grid1d {
  private:
    //==================================================================
    /// @name constant_values
    ///@{
    /// @brief big is the maximum value of double
    static constexpr double big = std::numeric_limits<double>::max();
    ///@} ------------------------------------------------------------------

    std::string name = "none";
    /// @brief zmin/zmax is the min/max of the z value
    double zmin = big;
    double zmax = -big;
    /// @brief vec_z is the data array of some value
    std::vector<double> vec_z;

  public:
    //==================================================================
    /// @name constructor_desctructor
    ///@{
    
    /// @brief default constructor
    Grid1dXZ() = default;
    
    /// @brief copy constructor
    Grid1dXZ(const Grid1dXZ &org) = default;

    /// @brief move constructor
    Grid1dXZ(Grid1dXZ &&other) noexcept = default;

    /// @brief destructor
    virtual ~Grid1dXZ() = default;

    /// @brief constructor from file, which calls \ref build_from_ascii_xy()
    /// @param[in] path_in Input file path containing (x, z) pairs
    /// @param[in] tf_shift If true, grid bounds are shifted by ±0.5*interval
    /// @throws std::runtime_error If file does not exist or x-values are not uniformly spaced
    Grid1dXZ(const std::filesystem::path &path_in,
             const bool tf_shift);
    ///@} ------------------------------------------------------------------
    
    /// @brief assignment operator
    Grid1dXZ& operator=(const Grid1dXZ& other) = default;
    
    /// @brief read data from file
    /// @details
    /// Reads (x, z) pairs from ASCII file, validates uniform x-spacing, and initializes the grid.
    /// The file format is space/tab-separated pairs: "x0 z0\nx1 z1\n..."
    /// @param[in] path_in Input file path
    /// @param[in] tf_shift If true, min becomes xmin-interval*0.5, max becomes xmax+interval*0.5
    /// @param[in] interval_precision Tolerance for checking uniform x-spacing (default: 1.0e-4)
    /// @throws std::runtime_error If file does not exist, cannot be read, or x-values are not uniformly spaced
    /// @note After successful load, zmin and zmax are computed from the data
    void build_from_ascii_xy(
        const std::filesystem::path &path_in
      , const bool tf_shift
      , const double interval_precision = 1.0e-4 );

    //======================================================================
    /// @name getter_funtions
    ///@{
    
    /// @brief get z value from index
    /// @param[in] index Grid bin index (must be 0 <= index < nbin)
    /// @return Z-value at the specified index
    /// @throws std::runtime_error If index is out of range
    double get_z( const int index ) const;

    /// @brief get z value from x value
    /// @param[in] value X-coordinate (must satisfy get_min() <= value < get_max())
    /// @return Z-value at the grid point corresponding to value
    /// @throws std::runtime_error If value is out of grid range
    double get_z( const double value ) const;

    /// @brief return zmin and zmax
    /// @return Tuple (zmin, zmax) computed by scanning all z-values
    /// @note Complexity: O(nbin)
    std::tuple<double,double> get_zmin_zmax() const;

    /// @brief return zmin, x_zmin, zmax, x_zmax
    /// @details Finds minimum and maximum z-values and their corresponding x-positions
    /// @return std::array<double,4> {zmin, x_zmin, zmax, x_zmax}
    /// @throws std::runtime_error If grid is empty (no bins processed)
    /// @note Complexity: O(nbin)
    std::array<double,4> get_zmin_x_zmin_zmax_x_zmax() const;

    /// @brief get name of this instance
    std::string get_name() const { return name; };

    /// @brief get zmin
    double get_zmin() const { return zmin; };

    /// @brief get zmax
    double get_zmax() const { return zmax; };

    /// @brief get z_value using linear interpolation from x=value
    /// @param[in] value X-coordinate (must satisfy get_min() <= value < get_max())
    /// @return Linearly interpolated z-value at the given x-coordinate
    /// @throws std::runtime_error If value is outside the grid range
    /// @note Performs piecewise linear interpolation between adjacent grid points
    double get_linear_interpolated_z( const double value ) const;

    /// @brief get x_value from z_value using inverse linear interpolation
    /// @details Finds all x-positions where the interpolated z equals the target z_value.
    /// Multiple solutions are possible for non-monotonic functions.
    /// @param[in] z_value Target z-value to search for
    /// @return Vector of x-coordinates where z(x) == z_value (empty if none found)
    /// @note Returns empty vector if z_value is outside [zmin, zmax]
    /// @note Complexity: O(nbin)
    std::vector<double> get_linear_interpolated_x( const double z_value ) const;

    /// @brief compute numerical derivative dz/dx
    /// @details Creates a new Grid1dXZ containing dz/dx at each point using finite differences.
    /// The last point is assigned the same derivative as the second-to-last point.
    /// @return New Grid1dXZ object containing derivative values
    /// @note Writes temporary output to "g1d_diff.tmp"
    Grid1dXZ make_differential() const;

    /// @brief create vector of (10^x, 10^z) pairs
    /// @return Vector of pairs where each pair is (pow(10, x), pow(10, z))
    /// @note Useful for logarithmic-scale data processing
    std::vector<std::pair<double,double>> make_vec_double2_pow10_xz() const;

    /// @brief create vector of (x, z) pairs
    /// @return Vector of (x, z) coordinate pairs for all grid points
    std::vector<std::pair<double,double>> make_vec_double2_xz() const;

    ///@} ------------------------------------------------------------------


    //======================================================================
    /// @name setter_functions
    ///@{

    /// @brief set zmin
    /// @param[in] zmin_in New minimum z-value
    /// @note Does not validate against actual z-values in vec_z
    void set_zmin( const double zmin_in ){ zmin = zmin_in; };

    /// @brief set zmax
    /// @param[in] zmax_in New maximum z-value
    /// @note Does not validate against actual z-values in vec_z
    void set_zmax( const double zmax_in ){ zmax = zmax_in; };

    /// @brief set z=z_in at index
    /// @param[in] index Grid bin index (must be 0 <= index < nbin)
    /// @param[in] z_in New z-value
    /// @throws std::runtime_error If index is out of range
    void set_z( const int index, const double z_in );

    /// @brief set z=z_in at x=value
    /// @param[in] value X-coordinate (must satisfy get_min() <= value < get_max())
    /// @param[in] z_in New z-value
    /// @throws std::runtime_error If value or computed index is out of range
    void set_z( const double value, const double z_in );

    /// @brief multiply all z values by factor
    /// @param[in] factor Multiplication factor
    void multiply_z( const double factor );

    /// @brief add constant to all z values
    /// @param[in] z_in Value to add to each z
    void add_z( const double z_in );

    ///@} ------------------------------------------------------------------

    /// @brief allocate memory for vec_z to match grid size
    /// @note Resizes vec_z to Grid1d::get_nbin()
    void vec_memory_allocate();

    /// @brief output (x, z) pairs to ASCII file
    /// @param[in] pathout Output file path
    /// @throws std::runtime_error If pathout is empty, file cannot be opened, or fclose fails
    /// @note Format: "x z\n" with scientific notation (%E)
    void out(const std::filesystem::path& pathout) const;

};

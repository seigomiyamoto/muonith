/// @file cls_Grid2dXYZ.hpp
/// @brief 2D grid that stores a scalar Z value per (x, y) bin
///
/// @details This file defines Grid2dXYZ, a Grid2d-derived container that associates
///          a scalar value (z) with each 2D grid bin.
///
///          Typical workflow:
///          1. Load XYZ data from ASCII (build_from_ascii_xyz) or construct from axes
///          2. Allocate storage with vec_vec_memory_allocate
///          3. Query values or interpolated values for analysis
///          4. Export to ASCII or binary if needed
///
///          Coordinate system:
///          - x/y axes follow Grid2d conventions (right-handed XY plane)
///          - z is a scalar value; units depend on dataset and are documented per method
///
///          Memory layout:
///          - vec_vec_z is column-major: vec_vec_z.at(iy).at(ix)
///          - Outer vector: y-direction (size = nbiny)
///          - Inner vector: x-direction (size = nbinx)
///
///          Thread safety:
///          - Read-only (const) methods are thread-safe after construction
///          - Mutating methods are not thread-safe
///          - build_from_ascii_xyz sorts XYZ data during loading
///
///          I/O format:
///          - ASCII: lines of "x y z" in scientific notation
///          - Binary: [count] followed by (x, y, z) triples
#pragma once

#include <map>
#include <fstream>
#include <iostream>
#include <sstream> // istringstream
#include <string>
#include <cstdio>
#include <cmath>
#include <functional>  //for sorting
#include <algorithm>//for sorting

// Skip due to potential conflict with Intel oneAPI headers.
// #include <execution> // for parallel execution of sort 

#include <vector>
#include <filesystem> // for std::filesystem::path
#include <omp.h> // for OpenMP

#include <Eigen/Dense>
#include "cls_Grid2d.hpp"
#include "cls_Grid1dXZ.hpp"

namespace fs = std::filesystem;

//##################################################################################
//##################################################################################
/// @class Grid2dXYZ
/// @brief 2D scalar grid that maps z values onto Grid2d bins
///
/// @details This class stores a scalar z value per (x, y) bin using a
///          column-major 2D vector. It supports ASCII/binary I/O and
///          bilinear interpolation for continuous queries.
///
///          Example:
///          @code
///          Grid1d x_axis("x", 100, 0.0, 1000.0, 10.0);
///          Grid1d y_axis("y", 50, 0.0, 500.0, 10.0);
///          Grid2dXYZ grid(x_axis, y_axis);
///          grid.set_name("demo");
///          grid.initialize_z(0.0);
///          double z = grid.get_bilinear_interpolated_z_value(123.4, 56.7);
///          @endcode
/// @ingroup basicGridClasses
//##################################################################################
//##################################################################################
class Grid2dXYZ : public Grid2d {
  private:
    /// @brief name of this instance
    std::string name = "none";

    /// @brief scalar grid values stored as a 2D vector
    /// @remark Index order is vec_vec_z.at(iy).at(ix) (column-major).
    std::vector<std::vector<double>> vec_vec_z = {};

    /// @brief Lower clamp bound for DL in get_bilinear_interpolated_peneflux [kg/m^2].
    /// @remark NaN means "not set"; the grid's own y-axis minimum is used instead.
    double d_DL_clamp_min = std::nan("");

    /// @brief Upper clamp bound for DL in get_bilinear_interpolated_peneflux [kg/m^2].
    /// @remark NaN means "not set"; the grid's own y-axis maximum is used instead.
    double d_DL_clamp_max = std::nan("");

    /// @brief Load grid data from a g2zbin file (called by constructor).
    /// @details g2zbin version 2 stores canonical (already shifted) axis
    ///          coordinates, so no shift interpretation is needed on load.
    /// @param[in] path_in Input file path.
    /// @param[in] tolerance_ratio Tolerance for Grid1d interval consistency
    ///            checks (default: 1.0e-4).
    /// @throws std::runtime_error If file read, magic/version validation, or
    ///         grid construction fails.
    void load_g2zbin(const fs::path& path_in,
                     double tolerance_ratio = 1.0e-4);

  public:
    //======================================================================
    /// @name constructor_destructor
    ///@{

    /// @brief default constructor
    Grid2dXYZ() = default;

    /// @brief copy constructor
    Grid2dXYZ(const Grid2dXYZ &org) = default;

    /// @brief move constructor
    Grid2dXYZ(Grid2dXYZ &&other) noexcept = default;

    /// @brief destructor
    ~Grid2dXYZ() = default;

    /// @brief Construct from x and y axes.
    /// @param g1_x_in X axis definition.
    /// @param g1_y_in Y axis definition.
    /// @note Thread-safe: No (mutating).
    Grid2dXYZ(const Grid1d &g1_x_in, const Grid1d &g1_y_in)
      : Grid2dXYZ() {
      set_x_axis(g1_x_in);
      set_y_axis(g1_y_in);
      vec_vec_memory_allocate();
    }

    /// @brief Construct from data file (ASCII xyz or g2zbin).
    /// @details Dispatches based on file extension:
    ///          - `.g2zbin` -> load_g2zbin() (direct header read; tf_shift flags unused
    ///            because v2 files store canonical, already shifted axes)
    ///          - otherwise -> build_from_ascii_xyz() (sort + axis detection via set_xy_axis_from_vec_xyz)
    /// @param path_in Path to the data file.
    /// @param tf_shift_x If true, shift x-axis by half bin width (ASCII path only).
    /// @param tf_shift_y If true, shift y-axis by half bin width (ASCII path only).
    /// @param tolerance_ratio Tolerance ratio for grid interval detection (default=1.0e-4).
    /// @throws std::runtime_error if file read or axis detection fails.
    /// @note Thread-safe: No (mutating).
    Grid2dXYZ(const fs::path &path_in,
              const bool tf_shift_x, const bool tf_shift_y,
              const double tolerance_ratio = 1.0e-4);

    ///@} ------------------------------------------------------------------


    //======================================================================
    /// @name operators
    /// @{

    /// @brief Assignment operator.
    /// @note Thread-safe: No (mutating).
    Grid2dXYZ& operator=(const Grid2dXYZ& other) = default;

    /// @brief Not-equal operator; the name is not compared.
    /// @return True if any Grid2d or z-grid values differ.
    /// @note Thread-safe: Yes (read-only).
    /// @note Time complexity: O(nbinx * nbiny).
    bool operator!=(const Grid2dXYZ& other) const;

    /// @brief Equal operator defined by not-equal.
    /// @return True if all Grid2d and z-grid values are equal.
    /// @note Thread-safe: Yes (read-only).
    bool operator==(const Grid2dXYZ& other) const {
      return !(*this != other);
    };
    ///@} ------------------------------------------------------------------
    
    //======================================================================
    /// @name getter_Grid2dXYZ
    /// @{
    
    /// @brief Get the name of this instance.
    /// @return Instance name.
    /// @note Thread-safe: Yes (read-only).
    std::string get_name() const { return name; };

    /// @brief Get a copy of the z-value grid.
    /// @return Copy of the internal z grid (column-major).
    /// @note Thread-safe: Yes (read-only).
    /// @note Time complexity: O(nbinx * nbiny).
    const std::vector<std::vector<double>> get_vec_vec_z() const { return vec_vec_z; };

    /// @brief Get z value at (ix, iy).
    /// @param ix Grid index in x-direction.
    /// @param iy Grid index in y-direction.
    /// @return z value at (ix, iy).
    /// @throws std::runtime_error if indices are out of range.
    /// @note Thread-safe: Yes (read-only).
    double get_z( const int ix, const int iy ) const;

    /// @brief Get z value at (x, y).
    /// @param x X coordinate (same unit as x-axis).
    /// @param y Y coordinate (same unit as y-axis).
    /// @return z value at (x, y) by nearest bin index.
    /// @throws std::runtime_error if (x, y) is out of range.
    /// @note Thread-safe: Yes (read-only).
    double get_z( const double x, const double y ) const;

    /// @brief Get z value using bilinear interpolation in x-y plane.
    /// @param x X coordinate (same unit as x-axis).
    /// @param y Y coordinate (same unit as y-axis).
    /// @return Interpolated z value.
    /// @throws std::runtime_error if (x, y) is out of range.
    /// @note Thread-safe: Yes (read-only).
    /// @note Time complexity: O(1).
    double get_bilinear_interpolated_z_value(const double x, const double y ) const;
    
    /// @brief Get penetrating flux value by bilinear interpolation.
    /// @param costhz Cosine of zenith angle (dimensionless).
    /// @param DL_in Overburden (density length) in kg/m^2.
    /// @return Penetrating flux (linear scale) derived from log10(z) storage.
    /// @throws std::runtime_error if costhz is outside [0, 1].
    /// @note Thread-safe: Yes (read-only).
    /// @note This is meaningful only when z stores log10(flux) with
    ///       x = costhz and y = DL.
    /// @note DL_in is clamped to the bounds set via set_DL_clamp_bounds();
    ///       when unset, the grid's own y-axis bounds are used, so the read
    ///       never falls outside this table.
    double get_bilinear_interpolated_peneflux(
      const double costhz, const double DL_in ) const;

    ///@} ------------------------------------------------------------------

    //======================================================================
    /// @name setter_Grid2dXYZ
    /// @{
    
    /// @brief Set instance name from C-string.
    /// @param name_in Name string.
    void set_name( const char *name_in ){ name = name_in; };
    /// @brief Set instance name from std::string.
    /// @param name_in Name string.
    void set_name( std::string name_in ){ name = name_in; };

    /// @brief Set explicit DL clamp bounds used by get_bilinear_interpolated_peneflux().
    /// @param DL_min_in Lower clamp bound in kg/m^2.
    /// @param DL_max_in Upper clamp bound in kg/m^2.
    /// @note When not set (default: NaN), the grid's own y-axis bounds are used.
    /// @note Thread-safe: No (mutating).
    void set_DL_clamp_bounds( const double DL_min_in, const double DL_max_in ){
      d_DL_clamp_min = DL_min_in;
      d_DL_clamp_max = DL_max_in;
    };

    /// @brief Set z value at (ix, iy).
    /// @param ix Grid index in x-direction.
    /// @param iy Grid index in y-direction.
    /// @param z_in Z value to assign.
    /// @throws std::runtime_error if indices are out of range.
    /// @note Thread-safe: No (mutating).
    void set_z( const int ix, const int iy, const double z_in ) {
      // if out of range
      check_ix_inside(ix);
      check_iy_inside(iy);
      vec_vec_z.at(iy).at(ix) = z_in;
    };

    /// @brief Set z value at (x, y).
    /// @param x X coordinate (same unit as x-axis).
    /// @param y Y coordinate (same unit as y-axis).
    /// @param z_in Z value to assign.
    /// @throws std::runtime_error if (x, y) is out of range.
    /// @note Thread-safe: No (mutating).
    void set_z( const double x, const double y, const double z_in ){
      const int ix = get_ix(x);
      check_ix_inside(ix);
      const int iy = get_iy(y);
      check_iy_inside(iy);
      set_z(ix,iy,z_in);
    };
    ///@} ------------------------------------------------------------------

    /// @brief Allocate memory for vec_vec_z and optionally initialize values.
    /// @param z_ini Initial value for all z entries.
    /// @note Thread-safe: No (mutating).
    void vec_vec_memory_allocate( const double z_ini=0.0 );

    /// @brief Read ASCII XYZ data, sort, set axes, and assign z values.
    /// @param path_in Path to ASCII file with x y z per line.
    /// @param tf_shift_x If true, shift x-axis by half bin width.
    /// @param tf_shift_y If true, shift y-axis by half bin width.
    /// @param tolerance_ratio Tolerance ratio for detecting regular grid spacing.
    /// @throws std::runtime_error if file read or axis detection fails.
    /// @note Thread-safe: No (mutating).
    void build_from_ascii_xyz( const fs::path &path_in
      , const bool tf_shift_x, const bool tf_shift_y
      , const double tolerance_ratio = 1.0e-4 );

    /// @brief Output to ASCII file (x y z format).
    /// @param pathout Output path.
    /// @throws std::runtime_error if file open or write fails.
    /// @note Sort order is y then x (iy outer loop).
    /// @note Uses bin center values for x and y.
    /// @note Thread-safe: Yes (read-only).
    void out( const fs::path& pathout ) const;

    /// @brief Output to ASCII file (x y z format), y-major order.
    /// @param pathout Output path.
    /// @param tf_xcnt If true, use x-axis center value; otherwise use lower value.
    /// @param tf_ycnt If true, use y-axis center value; otherwise use lower value.
    /// @throws std::runtime_error if file open or write fails.
    /// @note Sort order is y then x (iy outer loop).
    /// @note Thread-safe: Yes (read-only).
    void out_yx( const fs::path& pathout
      , const bool tf_xcnt, const bool tf_ycnt ) const;

    /// @brief Output to ASCII file (x y z format), x-major order.
    /// @param pathout Output path.
    /// @param tf_xcnt If true, use x-axis center value; otherwise use lower value.
    /// @param tf_ycnt If true, use y-axis center value; otherwise use lower value.
    /// @throws std::runtime_error if file open or write fails.
    /// @note Sort order is x then y (ix outer loop).
    /// @note Thread-safe: Yes (read-only).
    void out_xy( const fs::path& pathout
      , const bool tf_xcnt, const bool tf_ycnt ) const;

    /// @brief Output ASCII data with x-major order and y lower-bound values.
    /// @param pathout Output path.
    /// @throws std::runtime_error if file open or write fails.
    /// @note Sort order is x then y (ix outer loop).
    /// @note Thread-safe: Yes (read-only).
    void out_all_xy_order_ylower(
      const fs::path& pathout) const;

    /// @brief Output to ASCII file using name-based default path.
    /// @details Format: x y z, file name is "<name>.tmp".
    /// @throws std::runtime_error if file open or write fails.
    /// @note Thread-safe: Yes (read-only).
    void out() const {
      fs::path pathout = get_name() + ".tmp";
      out(pathout);
    };

    //======================================================================
    /// @name g2zbin_io
    /// @brief Compact binary I/O for uniform-grid z-raster (g2zbin format).
    /// @{

    /// @brief Save grid data in g2zbin format to a file path.
    /// @param[in] pathout Output file path (should use .g2zbin extension).
    /// @param[in] use_float64 If true, z values are stored as float64;
    ///            otherwise float32 (default).
    /// @throws std::runtime_error If file open or write fails.
    /// @note Thread-safe: Yes (read-only on this instance).
    /// @note g2zbin version 2 stores canonical (already shifted) axis
    ///       coordinates via Grid1d::save() (pass-through).
    void save_g2zbin(const fs::path& pathout, bool use_float64 = false) const;

    /// @brief Save grid data in g2zbin format to an open output stream.
    /// @param[in,out] ofs Output file stream (binary mode, must be open).
    /// @param[in] use_float64 If true, z values are stored as float64;
    ///            otherwise float32 (default).
    /// @throws std::runtime_error If stream write fails.
    /// @note Thread-safe: Yes (read-only on this instance).
    void save_g2zbin(std::ofstream& ofs, bool use_float64 = false) const;

    ///@} ------------------------------------------------------------------

    //======================================================================
    /// @name make_dFdR_table
    /// @{

    /// @brief Extract a y-slice (Grid1dXZ) at the nearest x-bin.
    /// @param x_in Query x coordinate (same unit as x-axis).
    /// @param tf_xcnt If true, use x-axis center for nearest search.
    /// @return Grid1dXZ along y with z values at nearest x-bin.
    /// @throws std::runtime_error if x_in is out of range.
    /// @note Thread-safe: Yes (read-only).
    Grid1dXZ get_Grid1dXZ_when_x_fixed( const double x_in, const bool tf_xcnt ) const;

    /// @brief Get an interpolated y-slice at specified x.
    /// @param x_in X coordinate (same unit as x-axis).
    /// @return Grid1dXZ slice along y with bilinear interpolation in x.
    /// @throws std::runtime_error if x_in is out of range.
    /// @note Thread-safe: Yes (read-only).
    Grid1dXZ get_Grid1dXZ_interp(const double x_in) const;

    /// @brief Create dF/dR table from dF/dE and range information.
    /// @details Use when this instance represents differential flux dF/dE.
    ///          F is penetrating muon flux, R is range in m.w.e. (kg/m^2).
    ///          R in g1_logR_logE is expressed in g/cm^2.
    // Grid2dXYZ make_dFdR_table(
    //   const Grid1dXZ &g1_logR_logE, const Grid2dXYZ &g2_log_peneflux_R ) const;

    /// @brief Create dF/dR table.
    /// @param g1_logR_logE z=log10(Range[kg/m^2]), x=log10(E[GeV]).
    /// @param g2_logdFdE_logE z=log10(dF/dE[m^-2 sr^-1 s^-1 GeV^-1]),
    ///                        y=log10(E[GeV]), x=cos(zenith angle).
    /// @return Grid2dXYZ z=dF/dR[kg^-1 sr^-1 s^-1], y=log10(E[GeV]),
    ///         x=cos(zenith angle).
    /// @note Thread-safe: Yes (read-only on inputs).
    static Grid2dXYZ make_dFdR_table( const Grid1dXZ  &g1_logR_logE
              , const Grid2dXYZ &g2_logdFdE_logE );

    /// @brief Create dF/dR table from penetrating muon flux table.
    /// @param g2_log_peneflux_R_costhz z=log10(F[m^-2 sr^-1 s^-1]),
    ///                                y=R[kg/m^2], x=cos(zenith angle).
    /// @return Grid2dXYZ dF/dR table in linear scale.
    /// @note Thread-safe: Yes (read-only on input).
    static Grid2dXYZ make_dFdR_table_from_peneflux(
      const Grid2dXYZ &g2_log_peneflux_R_costhz);

    ///@} ------------------------------------------------------------------

    /// @brief Convert z to log10(z).
    /// @return New Grid2dXYZ with log10-transformed values.
    /// @note Thread-safe: Yes (read-only).
    Grid2dXYZ make_log_z() const;
    
    /// @brief Convert log10(z) to linear z.
    /// @return New Grid2dXYZ with pow10-transformed values.
    /// @note Thread-safe: Yes (read-only).
    Grid2dXYZ make_pow10_z() const;

    /// @brief Convert z to log10(abs(z)).
    /// @return New Grid2dXYZ with log10(abs(z)) values.
    /// @note Thread-safe: Yes (read-only).
    Grid2dXYZ make_log_abs_z() const;

    /// @brief Initialize all z values to a constant.
    /// @param z_ini Value to assign.
    /// @note Thread-safe: No (mutating).
    void initialize_z( const double z_ini = 0.0 );

    /// @brief Compute dz/dy using forward differences.
    /// @details The last row (iy=nbiny-1) reuses the previous derivative.
    /// @return Grid2dXYZ containing dz/dy values.
    /// @note Thread-safe: Yes (read-only).
    Grid2dXYZ make_z_dzdy() const;

    /// @brief Convert all grid elements to a vector of {x, y, z}.
    /// @return Vector of size nbinx * nbiny using lower x/y values.
    /// @note Thread-safe: Yes (read-only).
    /// @note Time complexity: O(nbinx * nbiny).
    std::vector<std::array<double, 3>> make_vec_double3_xyz() const;


};

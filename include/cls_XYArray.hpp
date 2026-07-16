/// @file cls_XYArray.hpp
/// @brief Two-dimensional coordinate array handler
/// @details Defines XYArray class for managing collections of (x,y) coordinate pairs.
///
/// ## Typical Usage
/// 1. Construct from file or vector: `XYArray arr(path)` or `XYArray arr(vec_xy)`
/// 2. Data is automatically sorted by x in increasing order after loading
/// 3. Query interpolated values: `get_interp_y(x)` or `get_interp_x(y)`
/// 4. Apply transformations: `convert_x_to_pow10()`, `make_differential()`
/// 5. Export to text file: `out(path)`
///
/// ## Coordinate System
/// - No specific coordinate system assumed; x and y are abstract scalar values
/// - Units are application-dependent (not enforced by this class)
///
/// ## Thread Safety
/// - Not thread-safe. External synchronization required for concurrent access.
#pragma once

#include <algorithm>
#include <array>
#include <cmath>
#include <cstdio>
#include <filesystem>
#include <functional>
#include <limits>
#include <string>
#include <utility>
#include <vector>

#include <Eigen/Dense>
//##################################################################################
/// @class XYArray
/// @brief Container for a collection of (x,y) coordinate pairs
/// @ingroup basicTools
///
/// @details Manages an ordered collection of 2D coordinate pairs with support
/// for interpolation, differentiation, and power-of-10 transformations.
/// Data is typically stored sorted by x in increasing order.
///
/// @note Internal storage uses `std::vector<std::pair<double,double>>`.
///       Interpolation methods assume data is sorted appropriately.
///
/// ### Example
/// @code
/// // Load from file
/// XYArray arr("/path/to/data.bin");
/// // Query interpolated y value
/// std::vector<double> y_vals = arr.get_interp_y(5.0);
/// // Compute derivative
/// XYArray dydx = arr.make_differential();
/// @endcode
//##################################################################################
class XYArray {
  private:
    /// @brief Identifier name for this array (e.g., filename or label)
    std::string name = "none";

    /// @brief Minimum x value in the array
    double xmin = std::numeric_limits<double>::max();
    /// @brief Maximum x value in the array
    double xmax = -std::numeric_limits<double>::max();
    /// @brief Minimum y value in the array
    double ymin = std::numeric_limits<double>::max();
    /// @brief Maximum y value in the array
    double ymax = -std::numeric_limits<double>::max();

    /// @brief Internal storage of (x,y) coordinate pairs
    std::vector<std::pair<double,double>> vec_xy = {};

    /// @brief Sentinel value used internally during interpolation computation
    static constexpr double INITIAL_INTERPOLATION_VALUE = 931.931;
  
  public:
    //======================================================================
    /// @name Constructors and Destructor
    /// @{

    /// @brief Default constructor creating an empty XYArray
    XYArray() = default;

    /// @brief Copy constructor
    /// @param[in] org Source object to copy from
    XYArray(const XYArray &org) = default;

    /// @brief Construct from binary file
    /// @param[in] path_in Path to binary file containing (x,y) pairs
    /// @throws std::runtime_error If file cannot be read
    /// @note Data is automatically sorted by x in increasing order after loading.
    ///       Uses myapp::read_double2() for file parsing.
    XYArray(const std::filesystem::path &path_in);

    /// @brief Construct from vector of coordinate pairs
    /// @param[in] vec_double2 Vector of (x,y) pairs to copy
    /// @note Data is NOT automatically sorted; call sort_x_increasing() if needed.
    XYArray(const std::vector<std::pair<double,double>> &vec_double2);

    /// @brief Destructor
    ~XYArray() = default;
    /// @}

    //======================================================================
    /// @name Operators
    /// @{

    /// @brief Copy assignment operator
    /// @param[in] other Source object to copy from
    /// @return Reference to this object
    XYArray& operator=(const XYArray &other) = default;
    /// @}

    //======================================================================
    /// @name Getters
    /// @{

    /// @brief Get the identifier name
    /// @return Name string (e.g., filename or user-assigned label)
    std::string get_name() const { return name; }

    /// @brief Get minimum x value
    /// @return Minimum x coordinate in the array
    double get_xmin() const { return xmin; }

    /// @brief Get maximum x value
    /// @return Maximum x coordinate in the array
    double get_xmax() const { return xmax; }

    /// @brief Get minimum y value
    /// @return Minimum y coordinate in the array
    double get_ymin() const { return ymin; }

    /// @brief Get maximum y value
    /// @return Maximum y coordinate in the array
    double get_ymax() const { return ymax; }

    /// @brief Get x-range as pair
    /// @return Pair of (xmin, xmax)
    std::pair<double,double> get_xmin_xmax() const { return std::make_pair(xmin,xmax); }

    /// @brief Get y-range as pair
    /// @return Pair of (ymin, ymax)
    std::pair<double,double> get_ymin_ymax() const { return std::make_pair(ymin,ymax); }

    /// @brief Get bounding box as array
    /// @return Array of [xmin, xmax, ymin, ymax]
    std::array<double,4> get_xmin_xmax_ymin_ymax() const {
      return {xmin, xmax, ymin, ymax};
    }

    /// @brief Get number of coordinate pairs
    /// @return Size of the internal vector
    /// @note Complexity: O(1)
    size_t get_np() const { return vec_xy.size(); }

    /// @brief Get (x,y) coordinate pair at index
    /// @param[in] i Zero-based index
    /// @return Pair of (x, y) at the specified index
    /// @throws std::out_of_range If index is out of bounds (i < 0 or i >= size)
    /// @note Complexity: O(1)
    std::pair<double,double> get_xy(const int i) const { return vec_xy.at(i); }

    /// @brief Interpolate x values for a given y
    /// @param[in] y_in The y value to query
    /// @return Vector of interpolated x values (may contain multiple values)
    /// @throws std::runtime_error If y_in < ymin
    /// @throws std::runtime_error If y_in >= ymax
    /// @note Data must be sorted by y before calling. Creates internal copy and sorts.
    /// @note Complexity: O(n log n) due to internal sorting, then O(n) for interpolation
    std::vector<double> get_interp_x( const double y_in ) const;

    /// @brief Interpolate y values for a given x
    /// @param[in] x_in The x value to query
    /// @return Vector of interpolated y values (may contain multiple values)
    /// @throws std::runtime_error If x_in < xmin
    /// @throws std::runtime_error If x_in >= xmax
    /// @note Data must be sorted by x (default after construction from file).
    /// @note Complexity: O(n) linear search and interpolation
    std::vector<double> get_interp_y( const double x_in ) const;

    /// @brief Create copy with x,y converted to pow10(x),pow10(y)
    /// @return New XYArray with transformed coordinates
    XYArray get_converted_xy_to_pow10() const;
    /// @}
    
    //======================================================================
    /// @name Setters
    /// @{

    /// @brief Set the identifier name
    /// @param[in] name_in New name string
    void set_name( const std::string &name_in ){ name = name_in; }

    /// @brief Recompute min/max values for x and y from current data
    /// @note Call this after modifying coordinates via set_xy().
    /// @note Complexity: O(n)
    void set_maxmin_xy();

    /// @brief Set coordinate pair at index
    /// @param[in] i Zero-based index
    /// @param[in] x_in New x value
    /// @param[in] y_in New y value
    /// @throws std::out_of_range If index is out of bounds (i < 0 or i >= size)
    /// @note Does NOT update min/max; call set_maxmin_xy() afterward if needed.
    /// @note Complexity: O(1)
    void set_xy(const int i, const double x_in, const double y_in){
       vec_xy.at(i) = std::make_pair(x_in,y_in);
    }
    /// @}

    //======================================================================
    /// @name Sorting
    /// @{

    /// @brief Sort data by x in increasing order
    /// @note Uses mysort::sort_vec_double2_x_increasing() internally.
    /// @note Complexity: O(n log n)
    void sort_x_increasing();

    /// @brief Sort data by y in increasing order
    /// @note Uses mysort::sort_vec_double2_y_increasing() internally.
    /// @note Complexity: O(n log n)
    void sort_y_increasing();
    /// @}

    //======================================================================
    /// @name Transformations
    /// @{

    /// @brief Compute numerical derivative dy/dx
    /// @return New XYArray containing (x, dy/dx) pairs
    /// @throws std::runtime_error If consecutive x values are too close (dx < std::numeric_limits<double>::epsilon())
    /// @throws std::runtime_error If computed derivative magnitude exceeds std::numeric_limits<double>::max()
    /// @note Result has (n-1) points for input with n points.
    ///       Uses forward difference: dy/dx at x[i] = (y[i+1]-y[i])/(x[i+1]-x[i])
    /// @note Complexity: O(n)
    XYArray make_differential() const;

    /// @brief Transform x values: x -> pow(10, x)
    /// @note Updates min/max after transformation.
    void convert_x_to_pow10();

    /// @brief Transform y values: y -> pow(10, y)
    /// @note Updates min/max after transformation.
    void convert_y_to_pow10();

    /// @brief Transform both x and y: (x,y) -> (pow(10,x), pow(10,y))
    /// @note Updates min/max after transformation.
    void convert_xy_to_pow10();
    /// @}

    //======================================================================
    /// @name File Output
    /// @{

    /// @brief Write data to text file
    /// @param[in] pathout Output file path
    /// @throws std::runtime_error If path is empty
    /// @throws std::runtime_error If file cannot be opened for writing
    /// @throws std::runtime_error If fclose fails
    /// @note Output format: "%E %E\n" (scientific notation, space-separated)
    void out( const std::filesystem::path& pathout ) const;

    /// @brief Write data to default file (name + ".tmp")
    /// @throws std::runtime_error If file cannot be opened
    void out() const;
    /// @}
};

/// @file cls_Pillar.hpp
/// @brief Pillar (rectangular parallelepiped) class for 3D geometry
/// @details
/// Defines the Pillar class representing a 3D rectangular box with vertical extent [zmin, zmax)
/// and uniform density. This class is primarily used in terrain modeling where Grid2dPillar
/// places one Pillar instance at each horizontal grid point.
///
/// @note Units: All spatial coordinates (zmin, zmax) are in meters. Density units depend on application context.
/// @note Coordinate system: Right-handed z-up coordinate system (z-axis points upward).
/// @note Thread-safety: This class has no internal synchronization. Concurrent writes to the same instance are unsafe.
/// @note The z-interval is half-open: [zmin, zmax), meaning zmin is included but zmax is excluded.
#pragma once

#include <fstream>

//################################################################
//################################################################
/// @class Pillar
/// @brief Represents a rectangular parallelepiped (cuboid) placed at each grid point in Grid2dPillar.
///
/// @details
/// This class encapsulates a 3D rectangular box defined by:
/// - Vertical extent: [zmin, zmax) (half-open interval, meters)
/// - Uniform density: constant throughout the volume
///
/// **Responsibilities:**
/// - Store and provide access to vertical bounds (zmin, zmax) and density
/// - Test whether a given z-coordinate lies within the cuboid's vertical range
/// - Support value-semantic operations (copy, move, equality comparison)
///
/// **Invariants:**
/// - No automatic enforcement that zmin <= zmax (caller's responsibility)
/// - Density can be any finite double value (including zero or negative)
///
/// **Primary use case:**
/// In terrain modeling, Grid2dPillar maintains a 2D grid where each (x, y) grid point
/// has an associated Pillar defining the vertical structure and density at that location.
///
/// **Example usage:**
/// @code
/// // Create a cuboid from z=100m to z=150m with density 2.5 g/cm^3
/// Pillar c(2.5, 100.0, 150.0);
///
/// // Query vertical bounds and density
/// double z_lower = c.get_zmin();  // 100.0
/// double z_upper = c.get_zmax();  // 150.0
/// double rho = c.get_density();   // 2.5
///
/// // Test if a point is inside (remember: [zmin, zmax) is half-open)
/// bool inside = c.is_z_inside(125.0);  // true
/// bool at_top = c.is_z_inside(150.0);  // false (zmax is excluded)
/// @endcode
///
/// @ingroup terrainClasses
/// @ingroup geometryClasses
//################################################################
//################################################################
class Pillar {
  private:
    /// @brief Density of the cuboid material (units depend on application context, typically g/cm^3 or kg/m^3)
    double density = 0.0;

    /// @brief Lower bound of the vertical extent (meters, inclusive in the interval [zmin, zmax))
    double zmin = 0.0;

    /// @brief Upper bound of the vertical extent (meters, exclusive in the interval [zmin, zmax))
    double zmax = 0.0;

  public:

    //==================================================================
    /// @name constructor_destructor
    ///@{

    /// @brief Default constructor. Initializes all members to zero.
    Pillar() = default;

    /// @brief Parameterized constructor.
    /// @param[in] density_in Density value (units application-dependent)
    /// @param[in] zmin_in Lower z-bound (meters, inclusive)
    /// @param[in] zmax_in Upper z-bound (meters, exclusive)
    /// @note Does not validate that zmin_in <= zmax_in. Caller must ensure valid bounds.
    Pillar( const double density_in, const double zmin_in, const double zmax_in )
    : density(density_in), zmin(zmin_in), zmax(zmax_in) {};

    /// @brief Copy constructor (defaulted).
    Pillar( const Pillar& org ) = default;

    /// @brief Move constructor (defaulted, noexcept).
    Pillar( Pillar&& org ) noexcept = default;

    /// @brief Destructor (defaulted).
    ~Pillar() = default;

    ///@} ------------------------------------------------------------------

    //==================================================================
    /// @name operator
    ///@{

    /// @brief Copy assignment operator (defaulted).
    /// @param[in] other Source Pillar to copy from
    /// @return Reference to this object
    Pillar& operator=(const Pillar& other) = default;

    /// @brief Equality comparison operator.
    /// @param[in] other Pillar to compare with
    /// @return true if all members (density, zmin, zmax) are bitwise equal, false otherwise
    /// @note Uses exact floating-point equality (==), which may be inappropriate for computed values.
    bool operator==(const Pillar& other) const {
      return (density == other.density) && (zmin == other.zmin) && (zmax == other.zmax);
    };

    /// @brief Inequality comparison operator.
    /// @param[in] other Pillar to compare with
    /// @return true if any member differs, false if all are equal
    bool operator!=(const Pillar& other) const {
      return !(*this == other);
    };
    ///@} ------------------------------------------------------------------

    //==================================================================
    /// @name getter_functions
    ///@{

    /// @brief Get the lower z-bound.
    /// @return Lower z-coordinate (meters, inclusive in [zmin, zmax))
    double get_zmin() const { return zmin; };

    /// @brief Get the upper z-bound.
    /// @return Upper z-coordinate (meters, exclusive in [zmin, zmax))
    double get_zmax() const { return zmax; };

    /// @brief Get the density.
    /// @return Density value (units application-dependent)
    double get_density() const { return density; };
    ///@} ------------------------------------------------------------------


    //==================================================================
    /// @name setter_functions
    ///@{
    
    /// @brief Set all three member values.
    /// @param[in] zmin_in Lower z-bound (meters, inclusive)
    /// @param[in] zmax_in Upper z-bound (meters, exclusive)
    /// @param[in] density_in Density value
    /// @note Does not validate that zmin_in <= zmax_in.
    void set_values( const double zmin_in, const double zmax_in, const double density_in );

    /// @brief Set the lower z-bound.
    /// @param[in] zmin_in New lower z-bound (meters)
    void set_zmin( const double zmin_in){ zmin = zmin_in;};

    /// @brief Set the upper z-bound.
    /// @param[in] zmax_in New upper z-bound (meters)
    void set_zmax( const double zmax_in){ zmax = zmax_in;};

    /// @brief Set the density.
    /// @param[in] density_in New density value
    void set_density( const double density_in){ density = density_in; };

    /// @brief Copy all values from another Pillar.
    /// @param[in] org Source Pillar to copy from
    void set( const Pillar &org ){ set_values( org.zmin, org.zmax, org.density ); };

    /// @brief Add a delta to the current density.
    /// @param[in] delta_density Amount to add (can be negative to subtract)
    /// @note Performs density += delta_density.
    void add_density( const double delta_density){ density += delta_density;};
    ///@} ------------------------------------------------------------------
    
    //==================================================================
    /// @name bool_functions
    ///@{
    
    /// @brief Test if a z-coordinate is inside the cuboid's vertical extent.
    /// @param[in] z_in Z-coordinate to test (meters)
    /// @return true if z_in is in [zmin, zmax) (i.e., zmin <= z_in < zmax), false otherwise
    /// @note The interval is half-open: zmin is included, zmax is excluded.
    bool is_z_inside( const double z_in ) const;

    /// @brief Test if a z-coordinate is strictly below the upper bound.
    /// @param[in] z_in Z-coordinate to test (meters)
    /// @return true if z_in < zmax, false otherwise
    /// @note This does not check the lower bound (zmin).
    bool is_below_zmax( const double z_in ) const ;

    ///@} ------------------------------------------------------------------

    //==================================================================
    /// @name binary_io_functions
    ///@{

    /// @brief Save Pillar to binary stream.
    /// @param[out] ofs Output stream (must be opened in binary mode).
    void save(std::ofstream& ofs) const;

    /// @brief Load Pillar from binary stream.
    /// @param[in] ifs Input stream (must be opened in binary mode).
    void load(std::ifstream& ifs);

    ///@} ------------------------------------------------------------------

};

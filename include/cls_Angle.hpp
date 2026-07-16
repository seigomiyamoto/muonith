/// @file cls_Angle.hpp
/// @brief Angle class for handling angles in radians and degrees
/// @details
/// Provides the Angle class with automatic unit conversion between radians and degrees.
/// The class stores an angle value with an associated unit (Radian or Degree) and provides
/// conversion methods, arithmetic operators, and comparison operators.
///
/// @note Units: radians or degrees (user-specified per instance)
/// @note Thread-safety: No. Instances are not thread-safe for concurrent modification.
#pragma once

#include <string>
#include <cstdio>
#include <iosfwd>

//################################################################
//################################################################
/// @class Angle
/// @brief Represents an angle with support for Radian and Degree units.
/// @details
/// The Angle class stores an angle value along with its unit (Radian or Degree).
/// Conversion between units is performed automatically when calling rad() or deg().
/// All comparison operators normalize to radians for consistent comparisons.
///
/// Typical usage:
/// @code
/// Angle a1(45.0, Angle::Unit::Degree);
/// Angle a2(M_PI/4, Angle::Unit::Radian);
/// double rad_value = a1.rad();  // converts to radians
/// double deg_value = a2.deg();  // converts to degrees
/// @endcode
///
/// @note Thread-safety: No. Instances are not thread-safe.
/// @ingroup basicTools
//################################################################
//################################################################
class Angle {
  public:
  /// @brief Unit enumeration for Angle (Radian or Degree).
  /// @details Can be used independently of Angle instances: Angle::Unit angleUnit = Angle::Unit::Degree;
  enum class Unit { Radian, Degree };

  //==================================================================
  /// @name constructor_destructor
  ///@{

  /// @brief Construct an Angle with a specified value and unit.
  /// @param[in] value The angle value in the specified unit
  /// @param[in] unit The unit (Radian or Degree)
  /// @details Example: Angle angleInDegrees(45.0, Angle::Unit::Degree);
  Angle(const double value, const Unit unit);

  /// @brief Copy constructor (defaulted)
  Angle(const Angle &org) = default;

  /// @brief Move constructor (defaulted)
  Angle(Angle &&org) noexcept = default;

  /// @brief Destructor (defaulted)
  ~Angle() = default;
  ///@} ------------------------------------------------------------------

  //==================================================================
  /// @name getter_functions
  ///@{

  /// @brief Get the angle value in radians.
  /// @return The angle in radians (converted if stored in degrees)
  double rad() const;

  /// @brief Get the angle value in degrees.
  /// @return The angle in degrees (converted if stored in radians)
  double deg() const;

  /// @brief Get the current storage unit.
  /// @return The unit (Radian or Degree)
  Unit getUnit() const { return unit_; }
  ///@} ------------------------------------------------------------------

  /// @brief Convert a Unit enumeration value to a string.
  /// @param[in] unit The unit to convert
  /// @return "Radian", "Degree", or "Unknown"
  static std::string unitToString(Unit unit);

  //==================================================================
  /// @name setter_functions
  ///@{

  /// @brief Set the angle value in radians.
  /// @param[in] radian The angle value in radians
  void setRadian(const double radian);

  /// @brief Set the angle value in degrees.
  /// @param[in] degree The angle value in degrees
  void setDegree(const double degree);

  /// @brief Set the angle to the principal value of arctan(tangent), in the current unit.
  /// @param[in] tangent The tangent value
  /// @throws std::runtime_error If tangent is NaN
  /// @note The solution is θ = atan(tangent) + kπ. Only the principal value in (-π/2, π/2] is used.
  void setTangent(double tangent);

  /// @brief Negate the angle value (multiply by -1).
  void negate();
  ///@} ------------------------------------------------------------------

  // disp or output
  /// @brief Display the angle value to a FILE stream.
  /// @param[in] fout File pointer (default: stderr)
  /// @details Prints the angle in its current unit with the alternate unit in parentheses.
  void out(FILE *fout = stderr) const;

  //==================================================================
  /// @name operators
  ///@{

  /// @brief Assignment operator (defaulted)
  /// @param[in] rhs Right-hand side Angle
  /// @return Reference to this
  Angle& operator=(const Angle& rhs) = default;

  /// @brief Addition operator.
  /// @param[in] rhs Right-hand side Angle
  /// @return A new Angle with the sum, in the unit of this instance
  Angle operator+(const Angle& rhs) const;

  /// @brief Subtraction operator.
  /// @param[in] rhs Right-hand side Angle
  /// @return A new Angle with the difference, in the unit of this instance
  Angle operator-(const Angle& rhs) const;

  /// @brief Unary negation operator.
  /// @return A new Angle with negated value
  Angle operator-() const;

  /// @brief Scalar multiplication operator.
  /// @param[in] scalar The scalar multiplier
  /// @return A new Angle with value multiplied by scalar
  Angle operator*(const double scalar) const;

  /// @brief Scalar division operator.
  /// @param[in] scalar The scalar divisor
  /// @return A new Angle with value divided by scalar
  /// @throws std::runtime_error If scalar is zero
  Angle operator/(const double scalar) const;

  /// @brief Equality comparison operator.
  /// @param[in] rhs Right-hand side Angle
  /// @return True if angles are equal (compared in radians)
  bool operator==(const Angle& rhs) const {
    return rad() == rhs.rad();
  }

  /// @brief Inequality comparison operator.
  /// @param[in] rhs Right-hand side Angle
  /// @return True if angles are not equal (compared in radians)
  bool operator!=(const Angle& rhs) const {
    return !(*this == rhs);
  }

  /// @brief Less-than comparison operator.
  /// @param[in] rhs Right-hand side Angle
  /// @return True if this angle is less than rhs (compared in radians)
  bool operator<(const Angle& rhs) const {
    return rad() < rhs.rad();
  }

  /// @brief Greater-than comparison operator.
  /// @param[in] rhs Right-hand side Angle
  /// @return True if this angle is greater than rhs (compared in radians)
  bool operator>(const Angle& rhs) const {
    return rad() > rhs.rad();
  }

  /// @brief Less-than-or-equal comparison operator.
  /// @param[in] rhs Right-hand side Angle
  /// @return True if this angle is less than or equal to rhs (compared in radians)
  bool operator<=(const Angle& rhs) const {
    return rad() <= rhs.rad();
  }

  /// @brief Greater-than-or-equal comparison operator.
  /// @param[in] rhs Right-hand side Angle
  /// @return True if this angle is greater than or equal to rhs (compared in radians)
  bool operator>=(const Angle& rhs) const {
    return rad() >= rhs.rad();
  }
  ///@} ------------------------------------------------------------------

private:
  /// @brief Angle value in the stored unit
  double value_;

  /// @brief Storage unit (Radian or Degree)
  Unit unit_;
};

/// @brief Stream output operator for Angle.
/// @param[in,out] os Output stream
/// @param[in] angle The Angle to output
/// @return Reference to the output stream
/// @details Outputs the angle as "X.XX rad (Y.YY deg)"
std::ostream& operator<<(std::ostream& os, const Angle& angle);
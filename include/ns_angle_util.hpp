/// @file ns_angle_util.hpp
/// @brief Angle conversion and solid angle calculation utilities
/// @details
/// This module provides utilities for:
/// - Angle unit conversions (tangent, radian, degree)
/// - Solid angle calculations using various formulations
/// - Vector-to-angle and angle-to-vector transformations
/// - 3D rotation matrix generation (ZYX Euler angles)
///
/// @note Units: All angles are in radians unless explicitly stated otherwise.
/// @note Coordinate system: Right-handed, with typical conventions for azimuth/elevation.
/// @note Thread-safety: All functions are stateless and thread-safe.
#pragma once

#include <Eigen/Dense>
#include <cmath>
#include <string_view>
#include <string>
#include <stdexcept>

class Angle;

namespace angle_util {

/// @brief Rotation matrix type for 3D space
/// @note LOCAL: After rotation, the coordinate axes are attached to the object and rotate with it.
/// @note GLOBAL: After rotation, the coordinate axes remain fixed in space.
enum class Rotation3dType {
  LOCAL, GLOBAL
};

/// @brief Enumeration representing angle units
/// @details 0=Tangent (tan θ), 1=Radian, 2=Degree
enum class AngleUnit {
  Tangent = 0,
  Radian = 1,
  Degree = 2
};

/// @brief Convert AngleUnit to string representation
/// @param[in] unit The angle unit to convert
/// @return String representation ("Tangent", "Radian", or "Degree")
/// @throws std::runtime_error If the unit value is not recognized
std::string to_string(AngleUnit unit);

/// @brief Get AngleUnit from string (case-insensitive)
/// @param[in] name String representation of the angle unit
/// @return Corresponding AngleUnit value
/// @throws std::runtime_error If the string does not match any known unit
/// @note Accepts "Tangent", "Radian", "Degree" in any case (upper, lower, or mixed)
AngleUnit angle_unit_from_string(std::string_view name);

/// @brief Convert an angle value between AngleUnit representations
/// @param[in] value The angle value to convert
/// @param[in] from_unit The source unit
/// @param[in] to_unit The target unit
/// @return The converted angle value
double convert_angle_value(double value, AngleUnit from_unit, AngleUnit to_unit);

/// @brief Convert an angle value expressed in the given unit into radians
/// @param[in] value The angle value
/// @param[in] unit The unit of the input value
/// @return The angle in radians
double to_radian(double value, AngleUnit unit);

/// @brief Convert an angle value expressed in the given unit into degrees
/// @param[in] value The angle value
/// @param[in] unit The unit of the input value
/// @return The angle in degrees
double to_degree(double value, AngleUnit unit);

/// @brief Convert an angle value expressed in the given unit into tangent (tan θ)
/// @param[in] value The angle value
/// @param[in] unit The unit of the input value
/// @return The tangent of the angle
double to_tangent(double value, AngleUnit unit);


// Removed legacy conversion functions (not currently used):
// - Basic unit conversions: deg_to_rad, rad_to_deg, degree_to_tangent, tangent_to_rad, tangent_to_degree
// - Vector-to-angle conversions: uv_to_zenith_rad, uv_to_azimuth_rad (moved to cpp for internal use)
// - Transform basis directions: calc_radial, calc_lateral, calc_radial_lateral (moved to cpp for internal use)

/// @brief Convert tangent angles (tx, ty) to a 3D unit vector
/// @param[in] tan_x Tangent of the horizontal angle
/// @param[in] tan_y Tangent of the vertical angle
/// @return Normalized 3D vector
/// @note The vector is forward-pointing if tan_y >= 0, otherwise backward
Eigen::Vector3d tx_ty_vertical_to_v3_tangent(double tan_x, double tan_y);

/// @brief Convert azimuth (phi) and elevation angles to a 3D unit vector
/// @param[in] phi Azimuth angle
/// @param[in] elev Elevation angle
/// @return Normalized 3D vector
/// @throws std::runtime_error If elevation is outside [-π/2, π/2]
/// @note The vector is backward-pointing if the z-component would be negative
Eigen::Vector3d tx_ty_vertical_to_v3(Angle phi, Angle elev);

/// @brief Calculate solid angle using the (a, b, c) formulation
/// @param[in] a First parameter
/// @param[in] b Second parameter
/// @param[in] c Third parameter (typically normalized distance = 1.0)
/// @return Solid angle in steradians
/// @note This is a helper function used by the tangent-based solid angle calculation
double calc_omega_tangent(double a, double b, double c);

/// @brief Calculate solid angle over a rectangular region in tangent space
/// @param[in] txmin Minimum tangent of horizontal angle
/// @param[in] txmax Maximum tangent of horizontal angle
/// @param[in] tymin Minimum tangent of vertical angle
/// @param[in] tymax Maximum tangent of vertical angle
/// @return Solid angle in steradians
/// @throws std::runtime_error If the computed solid angle is non-positive
/// @note Uses the four-corner summation formula
double calc_omega_tangent(double txmin, double txmax, double tymin, double tymax);

/// @brief Calculate solid angle over a rectangular bin given in projective angles (radians)
/// @param[in] txmin Minimum horizontal projective angle atan(x/z) (radians), |txmin| < pi/2
/// @param[in] txmax Maximum horizontal projective angle atan(x/z) (radians), |txmax| < pi/2
/// @param[in] tymin Minimum vertical projective angle atan(y/z) (radians), |tymin| < pi/2
/// @param[in] tymax Maximum vertical projective angle atan(y/z) (radians), |tymax| < pi/2
/// @return Solid angle in steradians
/// @throws std::runtime_error If the computed solid angle is non-positive (via calc_omega_tangent)
/// @note Converts the bounds to tangent space and delegates to calc_omega_tangent
///       (exact four-corner formula). Replaces the former dphi x delev x sin(elev)
///       approximation, which underestimated near-horizon bins.
double calc_omega_radian(double txmin, double txmax, double tymin, double tymax);

/// @brief Calculate solid angle using Salerno's formula in tangent space
/// @param[in] txmin Minimum tangent of horizontal angle
/// @param[in] txmax Maximum tangent of horizontal angle
/// @param[in] tymin Minimum tangent of vertical angle
/// @param[in] tymax Maximum tangent of vertical angle
/// @return Solid angle in steradians
/// @note Alternative formulation for solid angle calculation
double calc_omega_tangent_alternative(double txmin, double txmax, double tymin, double tymax);

/// @brief Calculate the angle formed by two 3D vectors
/// @param[in] v3_0 First vector
/// @param[in] v3_1 Second vector
/// @return Angle between the vectors in radians [0, π]
/// @throws std::runtime_error If either vector has zero norm
double calc_angle_formed_by_vector3ds(const Eigen::Vector3d& v3_0, const Eigen::Vector3d& v3_1);

/// @brief Fast calculation of the angle between two 3D vectors with early rejection
/// @param[in] v3_0 First vector
/// @param[in] v3_1 Second vector
/// @param[in] cut_angle_rad Cutoff threshold for componentwise differences (radians)
/// @return Angle between vectors (radians), or -max(double) if any component difference exceeds threshold
/// @throws std::runtime_error If either vector has zero norm (when not early-rejected)
/// @note Performs early rejection by checking if any componentwise difference exceeds cut_angle_rad
double fast_calc_angle_formed_by_vector3ds(const Eigen::Vector3d& v3_0, const Eigen::Vector3d& v3_1, double cut_angle_rad);

/// @brief Generate a 3D rotation matrix using ZYX Euler angles
/// @param[in] theta_x Rotation angle about the X-axis
/// @param[in] theta_y Rotation angle about the Y-axis
/// @param[in] theta_z Rotation angle about the Z-axis
/// @param[in] rotation_type GLOBAL (fixed axes) or LOCAL (rotating axes)
/// @return 3x3 rotation matrix
/// @throws std::runtime_error If rotation_type is neither GLOBAL nor LOCAL
/// @note GLOBAL: R = Rx * Ry * Rz (axes remain fixed)
/// @note LOCAL: R = Rz * Ry * Rx (axes rotate with the object)
Eigen::Matrix3d make_rotation_3d_matrix_ZYX(const Angle& theta_x,
                                            const Angle& theta_y,
                                            const Angle& theta_z,
                                            const angle_util::Rotation3dType& rotation_type);

} // namespace angle_util

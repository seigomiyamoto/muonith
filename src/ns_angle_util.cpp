// ns_angle_util.cpp
#include "ns_angle_util.hpp"

#include <cassert>
#include <limits>
#include <Eigen/Geometry>

#include "cls_Angle.hpp"
#include "ns_myapp.hpp"
#include "ns_mymacro.hpp"
#include "ns_mylogger.hpp"
#include "spdlog_pch.hpp"

namespace angle_util {

std::string to_string(AngleUnit unit) {
  switch (unit) {
    case AngleUnit::Tangent: return "Tangent";
    case AngleUnit::Radian:  return "Radian";
    case AngleUnit::Degree:  return "Degree";
  }
  THROW_ERROR("angle_util::to_string: Unknown AngleUnit value");
}

AngleUnit angle_unit_from_string(std::string_view name) {
  if (name == "Tangent"  || name == "tangent"  || name == "TANGENT")  return AngleUnit::Tangent;
  if (name == "Radian"   || name == "radian"   || name == "RADIAN")   return AngleUnit::Radian;
  if (name == "Degree"   || name == "degree"   || name == "DEGREE")   return AngleUnit::Degree;
  THROW_ERROR("angle_util::angle_unit_from_string: unsupported value");
}

namespace {
constexpr double kPi = 3.14159265358979323846;

double to_radian_internal(double value, AngleUnit unit) {
  switch (unit) {
    case AngleUnit::Tangent: return std::atan(value);
    case AngleUnit::Radian:  return value;
    case AngleUnit::Degree:  return value * kPi / 180.0;
  }
  THROW_ERROR("angle_util::to_radian_internal: unknown angle unit");
}

double from_radian_internal(double value_rad, AngleUnit unit) {
  switch (unit) {
    case AngleUnit::Tangent: return std::tan(value_rad);
    case AngleUnit::Radian:  return value_rad;
    case AngleUnit::Degree:  return value_rad * 180.0 / kPi;
  }
  THROW_ERROR("angle_util::from_radian_internal: unknown angle unit");
}
}

double convert_angle_value(double value, AngleUnit from_unit, AngleUnit to_unit) {
  const double value_rad = to_radian_internal(value, from_unit);
  return from_radian_internal(value_rad, to_unit);
}

double to_radian(double value, AngleUnit unit) {
  return to_radian_internal(value, unit);
}

double to_degree(double value, AngleUnit unit) {
  return convert_angle_value(value, unit, AngleUnit::Degree);
}

double to_tangent(double value, AngleUnit unit) {
  return convert_angle_value(value, unit, AngleUnit::Tangent);
}

double calc_radial(double ax, double ay, double dx, double dy) {
  const double norm = std::sqrt(ax * ax + ay * ay);
  if (norm == 0.0) return 0.0;
  return (dx * ax + dy * ay) / norm;
}

double calc_radial(const Eigen::Vector2d& v2, const Eigen::Vector2d& v2_delta) {
  return calc_radial(v2.x(), v2.y(), v2_delta.x(), v2_delta.y());
}

double calc_lateral(double ax, double ay, double dx, double dy) {
  const double norm = std::sqrt(ax * ax + ay * ay);
  if (norm == 0.0) return 0.0;
  return (dx * ay - dy * ax) / norm;
}

double calc_lateral(const Eigen::Vector2d& v0, const Eigen::Vector2d& v_delta) {
  return calc_lateral(v0.x(), v0.y(), v_delta.x(), v_delta.y());
}

Eigen::Vector2d calc_radial_lateral(const Eigen::Vector2d& v2_base, const Eigen::Vector2d& v2_delta) {
  const double radial = calc_radial(v2_base, v2_delta);
  const double lateral = calc_lateral(v2_base, v2_delta);
  return {radial, lateral};
}

double uv_to_azimuth_rad(double vx, double vy, double vz) {
  if (vx == 0.0) {
    if (vy == 0.0) return 0.0;
    return kPi * 0.5;
  }
  return std::atan2(vy, vx);
}

double uv_to_azimuth_rad(const Eigen::Vector3d& v) {
  return uv_to_azimuth_rad(v.x(), v.y(), v.z());
}

Eigen::Vector3d tx_ty_vertical_to_v3_tangent(double tan_x, double tan_y) {
  if (tan_y >= 0) { // forward
    const double norm = 1.0 / std::sqrt(1.0 + tan_x * tan_x + tan_y * tan_y);
    const double vx = tan_x * norm;
    const double vy = 1.0 * norm;
    const double vz = tan_y * norm;
    return {vx, vy, vz};
  }
  const double norm = 1.0 / std::sqrt(1.0 + tan_x * tan_x + tan_y * tan_y);
  const double vx = -tan_x * norm;
  const double vy = -1.0 * norm;
  const double vz = -tan_y * norm;
  return {vx, vy, vz};
}

Eigen::Vector3d tx_ty_vertical_to_v3(Angle phi, Angle elev) {
  if (elev.rad() < -kPi * 0.5 || elev.rad() > kPi * 0.5) {
    THROW_ERROR2("angle_util::tx_ty_vertical_to_v3: elevation angle out of valid range [-pi/2, pi/2]", elev);
  }
  const double vx = std::cos(elev.rad()) * std::sin(phi.rad());
  const double vy = std::cos(elev.rad()) * std::cos(phi.rad());
  const double vz = std::sin(elev.rad());
  if (vz < 0) { // backward
    return {-vx, -vy, -vz};
  }
  return {vx, vy, vz};
}

double calc_omega_tangent(double a, double b, double c) {
  const double omega = std::asin((a * b) / (std::sqrt(b * b + c * c) * std::sqrt(c * c + a * a)));
  return omega;
}

double calc_omega_tangent(double txmin, double txmax, double tymin, double tymax) {
  double omega = 0.0;
  omega += calc_omega_tangent(txmax, tymax, 1.0);
  omega -= calc_omega_tangent(txmin, tymax, 1.0);
  omega -= calc_omega_tangent(txmax, tymin, 1.0);
  omega += calc_omega_tangent(txmin, tymin, 1.0);
  if (omega <= 0.0) {
    THROW_ERROR("angle_util::calc_omega_tangent: computed solid angle is non-positive. omega={}", omega);
  }
  return omega;
}

double calc_omega_radian(double txmin, double txmax, double tymin, double tymax) {
  // tx/ty are projective angles (tx = atan(x/z), ty = atan(y/z)), so the exact
  // solid angle is obtained by converting the bounds to tangent space and
  // delegating to the four-corner formula used by the Tangent mode.
  return calc_omega_tangent(
    std::tan(txmin), std::tan(txmax), std::tan(tymin), std::tan(tymax));
}

double calc_omega_tangent_alternative(double txmin, double txmax, double tymin, double tymax) {
  const double dtx = std::fabs(txmax - txmin);
  const double dty = std::fabs(tymax - tymin);
  const double tx = 0.5 * (txmin + txmax);
  const double ty = 0.5 * (tymin + tymax);
  const double lensq = 1 + tx * tx + ty * ty;
  return 1 / (lensq) / std::sqrt(lensq) * dtx * dty;
}

double calc_angle_formed_by_vector3ds(const Eigen::Vector3d& v3_0, const Eigen::Vector3d& v3_1) {
  const double norm0 = v3_0.norm();
  if (norm0 == 0.0) THROW_ERROR("angle_util::calc_angle_formed_by_vector3ds v3_0.norm()==0.0");
  const double norm1 = v3_1.norm();
  if (norm1 == 0.0) THROW_ERROR("angle_util::calc_angle_formed_by_vector3ds v3_1.norm()==0.0");
  const double vec_dot = v3_0.dot(v3_1);
  return std::acos(vec_dot / norm0 / norm1);
}

double fast_calc_angle_formed_by_vector3ds(const Eigen::Vector3d& v3_0, const Eigen::Vector3d& v3_1, double cut_angle_rad) {
  if (std::fabs(v3_0[0] - v3_1[0]) > cut_angle_rad) return -std::numeric_limits<double>::max();
  if (std::fabs(v3_0[1] - v3_1[1]) > cut_angle_rad) return -std::numeric_limits<double>::max();
  if (std::fabs(v3_0[2] - v3_1[2]) > cut_angle_rad) return -std::numeric_limits<double>::max();
  const double norm0 = v3_0.norm();
  const double norm1 = v3_1.norm();
  if (norm0 == 0.0) THROW_ERROR("angle_util::fast_calc_angle_formed_by_vector3ds v3_0.norm()==0.0");
  if (norm1 == 0.0) THROW_ERROR("angle_util::fast_calc_angle_formed_by_vector3ds v3_1.norm()==0.0");
  const double vec_dot = v3_0.dot(v3_1);
  return std::acos(vec_dot / norm0 / norm1);
}

Eigen::Matrix3d make_rotation_3d_matrix_ZYX(const Angle& theta_x,
                                            const Angle& theta_y,
                                            const Angle& theta_z,
                                            const angle_util::Rotation3dType& rotation_type) {
  Eigen::Matrix3d rotationMatrix;

  if (rotation_type == angle_util::Rotation3dType::GLOBAL) {
    const Eigen::AngleAxisd x_Angle(theta_x.rad(), Eigen::Vector3d::UnitX());
    const Eigen::AngleAxisd y_Angle(theta_y.rad(), Eigen::Vector3d::UnitY());
    const Eigen::AngleAxisd z_Angle(theta_z.rad(), Eigen::Vector3d::UnitZ());
    const Eigen::Quaterniond quat = x_Angle * y_Angle * z_Angle;
    rotationMatrix = quat.matrix();
  } else if (rotation_type == angle_util::Rotation3dType::LOCAL) {
    const Eigen::Vector3d euler_angles(theta_x.rad(), theta_y.rad(), theta_z.rad());
    rotationMatrix =
        Eigen::AngleAxisd(euler_angles[2], Eigen::Vector3d::UnitZ()).toRotationMatrix() *
        Eigen::AngleAxisd(euler_angles[1], Eigen::Vector3d::UnitY()).toRotationMatrix() *
        Eigen::AngleAxisd(euler_angles[0], Eigen::Vector3d::UnitX()).toRotationMatrix();
  } else {
    THROW_ERROR("angle_util::make_rotation_3d_matrix_ZYX: rotation_type must be LOCAL or GLOBAL");
  }
  return rotationMatrix;
}

} // namespace angle_util

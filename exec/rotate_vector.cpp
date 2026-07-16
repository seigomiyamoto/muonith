/// @file rotate_vector.cpp
/// @brief User-facing tool: rotate vectors by yaw/pitch/roll and show the resulting attitude.
/// @details
/// Reads a yaw/pitch/roll triple (degrees), a rotation frame (LOCAL or GLOBAL), and zero or
/// more 3D vectors. Builds a ZYX Euler rotation matrix with
/// angle_util::make_rotation_3d_matrix_ZYX and applies it to each vector, then reports the
/// rotated vector together with its azimuth and elevation. When no vector is given, the three
/// body axes (x, y, z) are rotated so that Euler-angle degeneracies such as gimbal lock
/// (pitch = +-90 deg, where yaw and roll act about the same axis) stay visible on a single run.
///
/// @note Angle mapping: yaw -> Z axis, pitch -> Y axis, roll -> X axis (ZYX Euler order).
/// @note Coordinate system: right-handed, z-up. Input angles are in degrees.
/// @note Frame: LOCAL rotates the axes with the object (R = Rz*Ry*Rx); GLOBAL keeps the axes
///       fixed in space (R = Rx*Ry*Rz). Both come from make_rotation_3d_matrix_ZYX.
/// @note Output goes through the project spdlog logger; no std::cout/std::cerr is added.
///
/// Usage:
///   rotate_vector.exe <yaw_deg> <pitch_deg> <roll_deg> <LOCAL|GLOBAL> [vx vy vz ...]

#include <cctype>
#include <cmath>
#include <numbers>
#include <string>
#include <utility>
#include <vector>

#include <Eigen/Dense>

#include "cls_Angle.hpp"
#include "ns_angle_util.hpp"
#include "ns_mylogger.hpp"
#include "ns_mymacro.hpp"

namespace {

/// @brief Parse a command-line token as a double.
/// @param[in] token The argument text.
/// @param[in] what Human-readable field name, used in the error message.
/// @return The parsed value.
/// @throws std::runtime_error If the token is not a full, in-range number.
double parse_double(const std::string& token, const char* what) {
  try {
    std::size_t pos = 0;                    // index of first unparsed char
    const double value = std::stod(token, &pos);
    if (pos != token.size()) {              // reject trailing garbage (e.g. "1.0abc")
      THROW_ERROR("rotate_vector: trailing characters in {} value. token={}", what, token);
    }
    return value;
  } catch (const std::invalid_argument&) {
    THROW_ERROR("rotate_vector: {} is not a number. token={}", what, token);
  } catch (const std::out_of_range&) {
    THROW_ERROR("rotate_vector: {} is out of range. token={}", what, token);
  }
  return 0.0;  // unreachable: THROW_ERROR always throws
}

/// @brief Parse the rotation frame token (case-insensitive).
/// @param[in] token "LOCAL" or "GLOBAL" in any case.
/// @return The matching Rotation3dType.
/// @throws std::runtime_error If the token is neither LOCAL nor GLOBAL.
angle_util::Rotation3dType parse_frame(const std::string& token) {
  std::string upper = token;
  for (char& c : upper) {
    c = static_cast<char>(std::toupper(static_cast<unsigned char>(c)));
  }
  if (upper == "LOCAL") return angle_util::Rotation3dType::LOCAL;
  if (upper == "GLOBAL") return angle_util::Rotation3dType::GLOBAL;
  THROW_ERROR("rotate_vector: frame must be LOCAL or GLOBAL. token={}", token);
  return angle_util::Rotation3dType::GLOBAL;  // unreachable: THROW_ERROR always throws
}

/// @brief Rotate one vector and log its input, output, azimuth and elevation.
/// @param[in] label Short name shown at the start of the line (e.g. "x-axis", "v1").
/// @param[in] v_in The vector before rotation.
/// @param[in] rotation_matrix The 3x3 rotation matrix to apply.
/// @note Azimuth = atan2(y, x); elevation = atan2(z, hypot(x, y)). Both are scale-invariant,
///       so a non-unit input still yields a meaningful attitude.
void log_rotation_result(const std::string& label,
                         const Eigen::Vector3d& v_in,
                         const Eigen::Matrix3d& rotation_matrix) {
  const Eigen::Vector3d v_out = rotation_matrix * v_in;
  constexpr double rad_to_deg = 180.0 / std::numbers::pi;
  const double azimuth_deg = std::atan2(v_out.y(), v_out.x()) * rad_to_deg;
  const double elevation_deg =
      std::atan2(v_out.z(), std::hypot(v_out.x(), v_out.y())) * rad_to_deg;
  LOG_INFO("{:<6} in ({: .4f}, {: .4f}, {: .4f}) -> out ({: .4f}, {: .4f}, {: .4f})"
           "  azimuth={: .2f} deg  elevation={: .2f} deg",
           label, v_in.x(), v_in.y(), v_in.z(), v_out.x(), v_out.y(), v_out.z(),
           azimuth_deg, elevation_deg);
}

}  // namespace

namespace {

/// @brief Core logic: parse arguments, build the rotation, and log each result.
/// @param[in] argc Argument count from main.
/// @param[in] argv Argument vector from main.
/// @return 0 on success, 1 on a usage error (too few arguments).
/// @throws std::runtime_error On an invalid frame, non-numeric input, or a bad vector count.
int run_rotate_vector(int argc, char** argv) {
  if (argc < 5) {
    LOG_ERROR("Usage: {} <yaw_deg> <pitch_deg> <roll_deg> <LOCAL|GLOBAL> [vx vy vz ...]", argv[0]);
    LOG_ERROR("  With no vector, the three body axes (x, y, z) are rotated.");
    return 1;
  }

  // Parse angles and frame.
  const double yaw_deg = parse_double(argv[1], "yaw_deg");
  const double pitch_deg = parse_double(argv[2], "pitch_deg");
  const double roll_deg = parse_double(argv[3], "roll_deg");
  const angle_util::Rotation3dType frame_type = parse_frame(argv[4]);

  // Map yaw/pitch/roll onto the ZYX Euler axes: yaw -> Z, pitch -> Y, roll -> X.
  const Angle theta_x(roll_deg, Angle::Unit::Degree);
  const Angle theta_y(pitch_deg, Angle::Unit::Degree);
  const Angle theta_z(yaw_deg, Angle::Unit::Degree);

  const Eigen::Matrix3d rotation_matrix =
      angle_util::make_rotation_3d_matrix_ZYX(theta_x, theta_y, theta_z, frame_type);

  // Echo the inputs and the resulting rotation matrix.
  const std::string frame_name =
      (frame_type == angle_util::Rotation3dType::LOCAL) ? "LOCAL" : "GLOBAL";
  LOG_INFO("yaw={} deg  pitch={} deg  roll={} deg  frame={}", yaw_deg, pitch_deg, roll_deg,
           frame_name);
  LOG_INFO("rotation matrix (row-major):");
  LOG_INFO("  [{: .4f}, {: .4f}, {: .4f}]", rotation_matrix(0, 0), rotation_matrix(0, 1),
           rotation_matrix(0, 2));
  LOG_INFO("  [{: .4f}, {: .4f}, {: .4f}]", rotation_matrix(1, 0), rotation_matrix(1, 1),
           rotation_matrix(1, 2));
  LOG_INFO("  [{: .4f}, {: .4f}, {: .4f}]", rotation_matrix(2, 0), rotation_matrix(2, 1),
           rotation_matrix(2, 2));

  // Build the list of vectors to rotate: user triples, or the three body axes by default.
  std::vector<std::pair<std::string, Eigen::Vector3d>> targets;
  const int n_vector_tokens = argc - 5;
  if (n_vector_tokens == 0) {
    targets.emplace_back("x-axis", Eigen::Vector3d(1.0, 0.0, 0.0));
    targets.emplace_back("y-axis", Eigen::Vector3d(0.0, 1.0, 0.0));
    targets.emplace_back("z-axis", Eigen::Vector3d(0.0, 0.0, 1.0));
    LOG_INFO("No vector given; rotating the three body axes so gimbal lock stays visible.");
  } else {
    if (n_vector_tokens % 3 != 0) {
      THROW_ERROR("rotate_vector: vector arguments must come in triples (vx vy vz). count={}",
                  n_vector_tokens);
    }
    const int n_vector = n_vector_tokens / 3;
    for (int i = 0; i < n_vector; ++i) {
      const int base = 5 + i * 3;
      const double vx = parse_double(argv[base + 0], "vx");
      const double vy = parse_double(argv[base + 1], "vy");
      const double vz = parse_double(argv[base + 2], "vz");
      targets.emplace_back("v" + std::to_string(i + 1), Eigen::Vector3d(vx, vy, vz));
    }
  }

  // Rotate each vector and report the resulting attitude.
  for (const auto& [label, v_in] : targets) {
    log_rotation_result(label, v_in, rotation_matrix);
  }

  return 0;
}

}  // namespace

int main(int argc, char** argv) {
  auto logger = mylogger::create_logger(spdlog::level::info, spdlog::level::warn,
                                        spdlog::level::trace, "logs/rotate_vector.log");
  // Convert a bad-input exception into a clean exit code instead of std::terminate.
  // The specific error was already logged at critical level by THROW_ERROR.
  try {
    return run_rotate_vector(argc, argv);
  } catch (const std::exception&) {
    return 1;
  }
}

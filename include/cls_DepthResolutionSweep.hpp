/// @file cls_DepthResolutionSweep.hpp
/// @brief Depth vs resolution sweep analysis class
/// @details
/// Performs sweep analysis to determine spatial resolution as a function of depth
/// for muon radiography density anomaly detection. Supports two modes:
/// - Cylinder mode: Uses explicit cylindrical geometry for density anomalies
/// - DDL mode: Delta Density-Length Approximation (Baumkuchen) for faster computation
///
/// Units: All spatial measurements in meters, angles in radians, density in kg/m³
/// Thread-safety: Not thread-safe; each instance should be used by a single thread
#pragma once

#include <string>
#include <vector>
#include <array>

#include <Eigen/Dense>
#include <nlohmann/json.hpp>
#include "cls_Angle.hpp"
#include "cls_FluxTable.hpp"
#include "cls_DetectorPanelArray.hpp"
#include "cls_DetectorPanelParameters.hpp"

class Grid2dPillar;
class DetectorPanelArray;

/// @brief Class for generating depth vs spatial resolution plots from surface
/// @details
/// This class performs systematic sweeps to evaluate the minimum detectable density
/// anomaly size as a function of depth below the surface. It computes statistical
/// significance of muon flux variations for various detector configurations.
///
/// Typical workflow:
/// 1. Construct with Parameters, FluxTable, and Grid2dPillar
/// 2. Call run() to execute the sweep (dispatches to mode-specific implementation)
/// 3. Results are written to ASCII files specified in Parameters
///
/// @note Not copyable or movable
/// @note Thread-safety: Not thread-safe
class DepthResolutionSweep {
public:
  /// @brief Structure for overriding angle bin configuration
  /// @details
  /// When tf_override_angle_bin is true, these values override the angle bins
  /// defined in detector JSON files. All detectors receive the same bin configuration.
  struct AngleBinOverride {
    int nbinx = 0;        ///< Number of x (horizontal angle) bins
    double txmin = 0.0;   ///< Minimum x angle value
    double txmax = 0.0;   ///< Maximum x angle value
    int nbiny = 0;        ///< Number of y (elevation angle) bins
    double tymin = 0.0;   ///< Minimum y angle value
    double tymax = 0.0;   ///< Maximum y angle value

    /// @brief Check if override values are set (valid)
    bool is_set() const { return nbinx > 0 && nbiny > 0; }
  };

  /// @brief Parameter structure for sweep configuration
  /// @details
  /// Contains all configurable parameters for the depth-resolution sweep.
  /// Parameters are typically loaded from JSON using from_json().
  ///
  /// @note Units: meters for spatial quantities, kg/m³ for density, radians for angles
  struct Parameters {
    enum class Mode { Cylinder, DDL }; ///< "cylinder" or "DDL(Delta Density-Length Approx.)"
    std::string section_name = "DEPTH_RESOLUTION_SWEEP"; ///< JSON section name
    DetectorPanel::ParameterLists det_param_lists; ///< Detector configuration for sweep

    /// @brief Flag to override angle bin settings from JSON
    bool tf_override_angle_bin = false;

    /// @brief Angle bin override values (used when tf_override_angle_bin is true)
    AngleBinOverride angle_bin_override;
    std::filesystem::path output_ascii_prefix;  ///< Output ASCII file prefix
    std::array<double, 2> obj_center{0.0, 0.0}; ///< Density anomaly center (x,y) [m]
    std::vector<double> vec_delta_density;      ///< Density difference candidates [kg/m3]
    double obj_size_upper_limit = -100.0;       ///< Density anomaly upper size limit [m]
    double obj_size_lower_limit = -100.0;       ///< Density anomaly lower size limit [m]
    double obj_size_step = -100.0;       ///< Density anomaly size step [m]
    double base_density = 2000.0;        ///< Base density [kg/m3]
    Mode mode = Mode::Cylinder; ///< Shape mode (default: Cylinder)

    /// @brief Elevation center angle step width, initial value is 0.01 rad
    double elev_center_step = 0.010;

    /// @brief Depth step width (meters)
    double depth_step = 5.0;

    /// @brief Cut panel tx in range of +-ang_between_unit_value * cut_factor
    double angle_between_cut_factor = 2.0;

    /// @brief Sweep in range of +-ang_between_unit_value * sweep_range_factor
    double sweep_range_factor = 1.0;

    /// @brief storage of statistical significance alpha values
    std::vector<double> vec_stat_alpha_value = {};

    /// @brief Signal/noise amplifier pairs [signal_amp, noise_amp]
    std::vector<std::array<double, 2>> vec_signal_noise_amplifiers = {};

    /// @brief Whether to use both-sided (two-sided) test. Default: false (one-sided).
    bool both_side = false;

    /// @brief Whether to output base PL/signal distribution. Default: false.
    bool tf_out_det_PL_signal = false;

    /// @brief Load parameters from JSON
    /// @param[in] root JSON root object
    /// @param[in] section JSON section name (default: "DEPTH_RESOLUTION_SWEEP")
    /// @return Populated Parameters structure
    /// @throws std::runtime_error If required fields are missing or invalid
    /// @note Validates all required fields and logs parameter values at DEBUG level
    static Parameters from_json(
        const nlohmann::json& root
      , const std::string& section = "DEPTH_RESOLUTION_SWEEP");

    /// @brief Utility to return obj_center as Eigen::Vector2d
    Eigen::Vector2d obj_center_vec() const {
      return Eigen::Vector2d(obj_center[0], obj_center[1]);
    }
  };

  //==================================================================
  /// @name json utility functions
  ///@{

  /// @brief Utility to retrieve double value from JSON
  /// @param j JSON object
  /// @param key  Key string
  /// @return Retrieved double value
  /// @throws std::runtime_error If key does not exist
  static double assign_double(const nlohmann::json& j, const char* key);

  /// @brief Utility to retrieve std::vector<double> value from JSON
  /// @param j JSON object
  /// @param key  Key string
  /// @return Retrieved std::vector<double> value
  /// @throws std::runtime_error If key does not exist
  static std::vector<double>
    assign_vec_double(const nlohmann::json& j, const char* key);

  /// @brief Utility to retrieve std::vector<std::array<double,2>> value from JSON
  /// @param j JSON object
  /// @param key  Key string
  /// @return Retrieved std::vector<std::array<double,2>> value
  /// @throws std::runtime_error If key does not exist
  static std::vector<std::array<double,2>>
    assign_vec_double2( const nlohmann::json& j, const char* key );

  /// @brief Utility to retrieve Angle value (radian) from JSON
  /// @param j JSON object
  /// @param key  Key string
  /// @return Retrieved Angle value
  /// @throws std::runtime_error If key does not exist
  static Angle assign_angle_rad(const nlohmann::json& j, const char* key);

  /// @brief Utility to retrieve double value from JSON pointer
  /// @param jptr JSON object pointer
  /// @param key  Key string
  /// @return Retrieved double value
  /// @throws std::runtime_error If key does not exist
  static double assign_double(const nlohmann::json* jptr, const char* key);

  /// @brief Utility to retrieve Angle value (radian) from JSON pointer
  /// @param jptr JSON object pointer
  /// @param key  Key string
  /// @return Retrieved Angle value
  /// @throws std::runtime_error If key does not exist
  static Angle assign_angle_rad(const nlohmann::json* jptr, const char* key);
  ///@} ------------------------------------------------------------------

  //==================================================================
  /// @name class functions
  ///@{

  /// @brief Constructor from parameters
  /// @param[in] params Parameter table with sweep configuration
  /// @param[in] ft_in Reference to FluxTable (must outlive this object)
  /// @param[in] g2pil_base Reference to base Grid2dPillar (copied internally)
  /// @param[in] arrdet_in Reference to DetectorPanelArray (copied internally)
  /// @throws std::runtime_error If params contains invalid configuration
  /// @note The FluxTable must remain valid for the lifetime of this object
  DepthResolutionSweep( const Parameters& params
  , const FluxTable& ft_in
  , const Grid2dPillar& g2pil_base
  , const DetectorPanelArray& arrdet_in);

  /// @brief copy constructor
  /// @param org original object
  DepthResolutionSweep(const DepthResolutionSweep& org) = delete;

  /// @brief assignment operator
  /// @param org original object
  DepthResolutionSweep& operator=(const DepthResolutionSweep& org) = delete;

  /// @brief move constructor
  /// @param org original object
  DepthResolutionSweep(DepthResolutionSweep&& org) = delete;

  /// @brief move assignment operator
  /// @param org original object
  DepthResolutionSweep& operator=(DepthResolutionSweep&& org) = delete;

  /// @brief destructor
  ~DepthResolutionSweep();

  /// @brief Calculate and set vec_vec_angbetween from arrdet_base_ and g2pil_base_ at params_.obj_center
  /// @param[in] obj_size_upper_limit Upper limit of object size [m]
  /// @details Computes angular resolution candidates for each detector based on
  ///          the object center position and maximum object size. For each detector,
  ///          generates a sequence of angular widths (ang_between) corresponding to
  ///          cylinder radii from the minimum bin width up to the maximum detectable size.
  /// @throws std::runtime_error If angle unit is unknown
  /// @note Units: obj_size_upper_limit in meters, angular values stored in radians internally
  void set_vec_vec_angbetween(const double obj_size_upper_limit);

  /// @brief Run sweep in cylinder mode
  /// @details Uses explicit vertical elliptic cylinder geometry to model density anomalies.
  ///          More accurate but computationally intensive. For each detector and angular size,
  ///          computes path length through cylindrical anomaly and evaluates statistical significance.
  /// @throws std::runtime_error If angle unit is unknown or file output fails
  /// @note Uses OpenMP for path length calculations (via pathcalc::g2pil::mp_add_PLDL and pathcalc::vcyl::mp_add_DL)
  /// @note Complexity: O(n_det * n_ang_between * n_delta_dens * n_ty_steps * n_bins)
  void run_mode_cylinder();

  /// @brief Run sweep in DDL (Delta Density-Length Approximation) mode
  /// @details Uses simplified delta density-length approximation (Baumkuchen) for faster computation.
  ///          Approximates path through density anomaly as uniform delta_DL = delta_dens * obj_size
  ///          across all bins within the sweep window. Less accurate than cylinder mode but faster.
  ///          This function sweeps over depth and obj_size combinations.
  /// @throws std::runtime_error If angle unit is unknown or file output fails
  /// @note DDL essence: The core idea is to add delta_dens * obj_size to ALL density-lengths
  ///       uniformly. This function sweeps that operation over detector, depth, and obj_size.
  /// @note Uses OpenMP for path length calculations (via pathcalc::g2pil::mp_add_PLDL)
  /// @note Complexity: O(n_det * n_ang_between * n_delta_dens * n_ty_steps * n_bins)
  void run_mode_ddl();

  /// @brief Optimised DDL sweep using precomputed full-panel PLDL + prefix sums
  /// @details
  /// Algorithmically equivalent to run_mode_ddl() but eliminates redundant
  /// computation by:
  ///   (1) Computing mp_add_PLDL on the FULL panel once per detector (not per ang_between)
  ///   (2) Using 2D prefix sums (PrefixSum2D) for O(1) rectangular range queries
  ///       in the y-sweep, replacing the O(nbinx × nbiny) loop
  ///   (3) Eliminating panel.cut() entirely (index ranges on the full panel instead)
  ///
  /// Expected speedup: ~100× for typical configurations (13 det × 160 ang_between).
  /// Output format is identical to run_mode_ddl().
  /// @throws std::runtime_error If angle unit is unknown, ty_step <= 0, or file output fails
  /// @note DDL essence: The core idea is to add delta_dens * obj_size to ALL density-lengths
  ///       uniformly. This function sweeps that operation over detector, depth, and obj_size.
  /// @note Uses OpenMP: parallel for over ang_between loop with dynamic scheduling
  /// @note Thread-safety: Inner operations use single-thread helpers (st_*) to avoid nested OpenMP
  void run_mode_ddl_v1();

  /// @brief Execute sweep using the mode specified in parameters
  /// @details Dispatches to run_mode_ddl() or run_mode_cylinder() based on params_.mode
  /// @throws std::runtime_error If mode is invalid
  void run();
  ///@} ------------------------------------------------------------------

private:
  /// @brief sweep parameters
  Parameters params_;

  /// @brief  spatial resolution candidates [rad] First index is Detid, second index is candidate number
  std::vector<std::vector<Angle>> vec_vec_angle_between_from_obj_center_;

  /// @brief reference of FluxTable
  const FluxTable& ft_;

  /// @brief reference of Grid2dPillar
  const Grid2dPillar& g2pil_base_;
  
  /// @brief base of DetectorPanelArray
  DetectorPanelArray arrdet_base_;
};

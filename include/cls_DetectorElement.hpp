/// @file cls_DetectorElement.hpp
/// @brief Detector element class representing individual angular bins in the muon detection system
///
/// @details
/// This file defines the DetectorElement class, which represents a single minimum angular
/// resolution element (angular bin) in a muon detection system. Each element stores:
/// - Geometric properties: position, direction (Ray3d), angular range (tx, ty)
/// - Physical measurements: effective area (m²), solid angle (sr), exposure time (sec)
/// - Muon flux data: path length (PL in meters), density length (DL in kg/m²),
///   penetrating muon flux (m⁻² sr⁻¹ sec⁻¹)
/// - Reconstruction results: projected density (kg/m³), signal/noise counts
/// - Efficiency parameters: lower/center/upper bounds for uncertainty sampling
///
/// @note Units:
/// - Length: meters (m)
/// - Density: kg/m³ (projected density), kg/m² (density length)
/// - Flux: m⁻² sr⁻¹ sec⁻¹
/// - Angles: configurable (Degree, Radian, or Tangent via AngleUnit enum)
///
/// @note Thread-safety:
/// This class is NOT thread-safe. External synchronization is required if instances
/// are shared across threads.
///
/// @note Coordinate system:
/// Assumes a right-handed coordinate system. Direction vectors are expected to be normalized.
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
#include <vector>
#include <mutex> // for std::mutex

#include <Eigen/Dense>
#include "cls_Ray.hpp"
#include "cls_Grid2dXYZ.hpp"
#include "cls_Grid3dVoxel.hpp"
#include "ns_angle_util.hpp"
#include "ns_type_definitions.hpp"
using namespace index_type_definitions;

//###############################################################
//###############################################################
/// @class DetectorElement
/// @brief Represents a single detector element with a small angular bin for muon tomography
/// @ingroup detectorClasses
///
/// @details
/// DetectorElement encapsulates all properties of a minimum angular resolution element
/// in a muon detection system. It serves as the fundamental unit for muon flux measurement
/// and density reconstruction.
///
/// **Responsibilities:**
/// - Store geometric configuration (position, direction, angular boundaries)
/// - Manage physical measurements (effective area, solid angle, exposure time)
/// - Track muon flux and path information through materials
/// - Calculate and store reconstructed densities with uncertainties
/// - Support efficiency sampling for Monte Carlo simulations
///
/// **Key invariants:**
/// - `angle_unit` determines interpretation of `txmin`, `txmax`, `tymin`, `tymax`
/// - `proj_density == invalid_proj_dens (-1.0)` indicates uninitialized/invalid density
/// - `vec_tf_in_PL` stores alternating in/out transitions (validated by check_alternating_tf_in)
///
/// **Typical workflow:**
/// @code
/// DetectorElement elem;
/// elem.set_uqid(0);
/// elem.set_detid(1);
/// elem.set_v3_position(Eigen::Vector3d(0.0, 0.0, 0.0));
/// elem.set_v3_direction(Eigen::Vector3d(0.0, 0.0, 1.0));
/// elem.set_txmin_txmax_tymin_tymax(-0.1, 0.1, -0.1, 0.1, AngleUnit::Tangent);
/// elem.set_effective_area_m2(1.0);
/// double omega = elem.calc_solid_angle();
/// elem.set_solid_angle(omega);
/// elem.set_exposure_time_sec(86400.0);
/// // ... perform ray tracing to accumulate PL, DL
/// elem.calc_set_proj_density(1.0, 1.0); // threshold PL=1m, DL=1kg/m²
/// @endcode
//###############################################################
//###############################################################
class DetectorElement final {
  public:
    /// @brief initial number of vec_tf_in_PL
    static constexpr int n_vec_tf_in_PL_const = 10;

    /// @brief AngleUnit definition
    using AngleUnit = angle_util::AngleUnit;

  private:
    /// @brief An unique number will be assigned for all detector elements.
    Uqid unique_index = UqidNotAssigned;
  
    /// @brief detector ID of this element
    Detid detid = DetidNotAssigned;
  
    /// @brief ID within the detector. 
    Inthis id_in_this_detector = InthisNotAssigned;
  
    /// @brief position and direction as Ray3d
    Ray3d ray3d = Ray3d();

    /// @brief txmin, txmax, tymin,tymax angle unit
    AngleUnit angle_unit = AngleUnit::Tangent;
  
    /// @brief (txmin,txmax,tymin,tymax) angle range. Unit depends on angle_unit.
    double txmin=0.0, txmax=0.0, tymin=0.0, tymax=0.0;
  
    /// @brief effective_area_m2 in this element 
    double effective_area_m2=0.0;

    /// @brief solid_angle in this element
    double solid_angle=0.0;

    /// @brief [sec], exposure_time_sec
    double exposure_time_sec=0.0;
  
    /// @brief rock length in the path (meters)
    double PL=0.0;
  
    /// @brief density length in the path (kg/m2)
    double DL=0.0;

    /// @brief penetrating_muon_flux
    /// @note SI units: m-2 sr-1 sec-1
    double penetrating_muon_flux=0.0;

    /// @brief muon signal
    double signal=0.0;

    /// @brief deterministic floor noise part (added without Poisson fluctuation)
    double noise_det=0.0;

    /// @brief Poisson-fluctuated noise part (value before the Poisson draw)
    double noise_poi=0.0;

    /// @brief projection density
    /// @note SI units: kg/m3
    double proj_density= -1.0;

    /// @brief projection density upper
    /// @note SI units: kg/m3
    double proj_density_upper= -1.0;

    /// @brief projection density lower
    /// @note SI units: kg/m3
    double proj_density_lower= -1.0;

    /// @brief lower efficiency value
    double eff_low=0.0;

    /// @brief center efficiency value
    double eff_cnt=0.0;

    /// @brief upper efficiency value
    double eff_upp=0.0;

    /// @brief store path length information with in/out flag
    /// @note tf_in_out==true is in, false is out
    std::vector<std::tuple<bool, double>> vec_tf_in_PL;

  public:
    //==================================================================
    /// @name static_constants
    ///@{

    /// @brief invalid_proj_dens
    static constexpr double invalid_proj_dens = -1.0;

    ///@} ------------------------------------------------------------------

    //==================================================================
    /// @name constructor_destructor
    ///@{
  
    /// @brief constructor
    /// @note if we want to use reserve for vec_tf_in_PL, we cannot use in-class initialization!
    DetectorElement()
    : unique_index(UqidNotAssigned)
    , detid(DetidNotAssigned)
    , id_in_this_detector(InthisNotAssigned)
    , ray3d()  // Ray3d default constructor
    , angle_unit(AngleUnit::Tangent)
    , txmin(0.0)
    , txmax(0.0)
    , tymin(0.0)
    , tymax(0.0)
    , effective_area_m2(0.0)
    , solid_angle(0.0)
    , exposure_time_sec(0.0)
    , PL(0.0)
    , DL(0.0)
    , penetrating_muon_flux(0.0)
    , signal(0.0)
    , noise_det(0.0)
    , noise_poi(0.0)
    , proj_density(-1.0)
    , proj_density_upper(-1.0)
    , proj_density_lower(-1.0)
    , eff_low(1.0)
    , eff_cnt(1.0)
    , eff_upp(1.0)
    { 
      vec_tf_in_PL.reserve(n_vec_tf_in_PL_const);
    };
  
    /// @brief copy constructor
    DetectorElement(const DetectorElement &org) = default;
  
    /// @brief move constructor
    DetectorElement(DetectorElement&&) noexcept = default;

    /// @brief destructor
    ~DetectorElement() = default;
  
    /// @brief virtual clone method
    virtual DetectorElement* clone() const { return new DetectorElement(*this); }

    ///@} ------------------------------------------------------------------

    //======================================================================
    /// @name operators
    ///@{

    /// @brief assignment operator
    DetectorElement& operator=(const DetectorElement& other) = default;

    /// @brief Inequality comparison operator
    /// @param[in] other DetectorElement to compare against
    /// @return true if any member differs, false if all members are equal
    /// @details
    /// Compares all member variables for inequality. In non-NODEBUG builds,
    /// logs which member differs using LOG_WARN before returning true.
    /// @note Follows project convention: NODEBUG fast path, non-NODEBUG diagnostic path
    bool operator!=(const DetectorElement& other) const;

    /// @brief Equality comparison operator
    /// @param[in] other DetectorElement to compare against
    /// @return true if all members are equal, false otherwise
    /// @details Implemented as negation of operator!=
    bool operator==(const DetectorElement& other) const {
      return !(*this != other);
    }

    /// @brief move assignment operator
    DetectorElement& operator=(DetectorElement&& other) noexcept = default;

    ///@} ------------------------------------------------------------------

    /// @brief Reset all members to default-constructed state
    /// @details Equivalent to *this = DetectorElement(). All indices set to NotAssigned, numeric values to zero.
    void reset(){ *this = DetectorElement(); };

    /// @brief Write detector element state to FILE stream in ASCII format
    /// @param[out] fout Output file pointer (must be open for writing)
    /// @details
    /// Outputs space-separated values: unique_index, direction (vx,vy,vz), position (x,y,z),
    /// angle_unit + angular range, effective_area_m2, solid_angle, exposure_time_sec,
    /// PL, DL, penetrating_muon_flux, signal, noise, proj_density, proj_density_lower,
    /// proj_density_upper, eff_low, eff_cnt, eff_upp, followed by newline.
    /// @note Does NOT write vec_tf_in_PL to ASCII output. Use save() for complete binary serialization.
    void out( FILE *fout ) const;

    //==================================================================
    /// @name getter_functions
    ///@{

    /// @brief make name from detid and id_in_this_detector
    std::string get_name() const;

    /// @brief get unique_index
    Uqid get_uqid() const { return unique_index; }
    
    /// @brief get detid
    Detid get_detid() const { return detid; }

    /// @brief get id_in_this_detector
    Inthis get_id_in_this_detector() const { return id_in_this_detector; }

    /// @brief get ray3d
    Ray3d get_ray3d() const { return ray3d; };

    /// @brief get ray2d
    /// @details const Ray2d ray2d = ray3d.toRay2d(); toRay2d  returns (v3_pos.head(2), v3_dir.head(2));
    Ray2d get_ray2d() const { return ray3d.toRay2d(); };

    /// @brief get v3_direction
    Eigen::Vector3d get_v3_dir() const { return ray3d.dir(); };

    /// @brief get v3_position
    Eigen::Vector3d get_v3_pos() const { return ray3d.pos(); };

    /// @brief get vx
    double get_vx() const { return ray3d.vx(); };

    /// @brief get vy
    double get_vy() const { return ray3d.vy(); };

    /// @brief get vz
    double get_vz() const { return ray3d.vz(); };

    /// @brief get x
    double get_x() const { return ray3d.x(); };

    /// @brief get y
    double get_y() const { return ray3d.y(); };

    /// @brief get z
    double get_z() const { return ray3d.z(); };

    /// @brief get 0.5*(txmin+txmax);
    double get_tx() const { return 0.5*(txmin+txmax); };

    /// @brief get 0.5*(tymin+tymax);
    double get_ty() const { return 0.5*(tymin+tymax); };

    /// @brief get txmin
    double get_txmin() const { return txmin; };

    /// @brief get txmax
    double get_txmax() const { return txmax; };

    /// @brief get tymin
    double get_tymin() const { return tymin; };

    /// @brief get tymax
    double get_tymax() const { return tymax; };

    /// @brief Get angular range with unit conversion
    /// @param[in] angle_unit_in Target angle unit for output
    /// @return Array {txmin, txmax, tymin, tymax} converted to angle_unit_in
    /// @details
    /// Converts the stored angular boundaries from the current angle_unit to angle_unit_in.
    /// Supports all conversions between Degree, Radian, and Tangent.
    /// If angle_unit == angle_unit_in, returns values as-is (no conversion).
    /// @throws std::runtime_error If angle_unit is not defined
    std::array<double,4> get_txmin_txmax_tymin_tymax(const AngleUnit& angle_unit_in) const;

    /// @brief get AngleUnit as string
    std::string get_angle_unit() const {
      return angle_util::to_string(angle_unit);
    };

    /// @brief get penetrating_muon_flux
    double get_peneflux() const {return penetrating_muon_flux; };

    /// @brief get effective_area_m2
    double get_effective_area_m2() const { return effective_area_m2; };

    /// @brief get solid_angle
    double get_solid_angle() const { return solid_angle; };

    /// @brief get exposure_time_sec
    double get_exposure_time_sec() const { return exposure_time_sec; };

    /// @brief get SOT with [m2, sr, sec]
    double get_SOT() const { return effective_area_m2 * solid_angle * exposure_time_sec; };

    /// @brief get path length in the rock (meters)
    double get_PL() const { return PL; };

    /// @brief get density length (kg/m2)
    double get_DL() const { return DL; };

    /// @brief get muon signal
    double get_signal() const { return signal; };

    /// @brief get total noise (deterministic floor + Poisson-bucket value)
    double get_noise() const { return noise_det + noise_poi; };

    /// @brief get deterministic floor noise part
    double get_noise_det() const { return noise_det; };

    /// @brief get Poisson-bucket noise part (value before the Poisson draw)
    double get_noise_poi() const { return noise_poi; };

    /// @brief get signal + total noise
    double get_signal_plus_noise() const { return signal + noise_det + noise_poi; };

    /// @brief get proj_density
    double get_proj_density() const { return proj_density; };

    /// @brief get proj_density_lower
    double get_proj_density_lower() const { return proj_density_lower; };

    /// @brief get proj_density_upper
    double get_proj_density_upper() const { return proj_density_upper; };

    /// @brief get eff_low
    double get_eff_low() const { return eff_low; };

    /// @brief get eff_cnt
    double get_eff_cnt() const { return eff_cnt; };

    /// @brief get eff_upp
    double get_eff_upp() const { return eff_upp; };

    /// @brief get immutable refrerence of vec_tf_in_PL
    const std::vector<std::tuple<bool,double>> &get_vec_tf_in_PL() const { return vec_tf_in_PL; };

    ///@} ------------------------------------------------------------------

    //==================================================================
    /// @name setter_functions
    ///@{
    
    /// @brief set uqid
    void set_uqid(const Uqid &uqid_in){ unique_index = uqid_in; };
    
    /// @brief set detid
    void set_detid(const Detid &detid_in){ detid = detid_in; };
    
    /// @brief set id_in_this_detector
    void set_id_in_this_detector(const Detid &id_in_this_detector_in){
      id_in_this_detector = id_in_this_detector_in; };
    
    /// @brief set v3_direction
    void set_v3_direction(const Eigen::Vector3d &v3_direction_in){
      ray3d.set_dir(v3_direction_in);
    };
    
    /// @brief set v3_position
    void set_v3_position(const Eigen::Vector3d &v3_position_in){
      ray3d.set_pos(v3_position_in);
    };

    /// @brief set ray3d
    void set_ray3d(const Ray3d &ray3d_in){
      ray3d = ray3d_in;
    };
    
    /// @brief Set angular range and angle unit
    /// @param[in] txmin_in Minimum x-angle (horizontal)
    /// @param[in] txmax_in Maximum x-angle (horizontal)
    /// @param[in] tymin_in Minimum y-angle (vertical)
    /// @param[in] tymax_in Maximum y-angle (vertical)
    /// @param[in] angle_unit_in Unit interpretation for the input angles
    /// @details All four angle values must be in the same unit specified by angle_unit_in
    void set_txmin_txmax_tymin_tymax(
        const double txmin_in, const double txmax_in
      , const double tymin_in, const double tymax_in
      , const AngleUnit angle_unit_in )
    {
      txmin = txmin_in;
      txmax = txmax_in;
      tymin = tymin_in;
      tymax = tymax_in;
      angle_unit = angle_unit_in;
    };

    /// @brief set effective_area_m2
    void set_effective_area_m2( const double effective_area_m2_in){
      effective_area_m2 = effective_area_m2_in;
    };

    /// @brief set solid_angle
    void set_solid_angle(const double solid_angle_in){
      solid_angle = solid_angle_in; };

    /// @brief set exposure_time_sec
    void set_exposure_time_sec(const double exposure_time_sec_in){
      exposure_time_sec = exposure_time_sec_in; };

    /// @brief set PL
    void set_PL( const double PL_in ){ PL = PL_in; };

    /// @brief PL += PL_in
    void add_PL( const double PL_in ){ PL += PL_in; };

    /// @brief set DL
    void set_DL( const double DL_in ){ DL = DL_in;  };

    /// @brief DL += DL_in
    void add_DL( const double DL_in ){ DL += DL_in; };

    /// @brief set penetrating_muon_flux
    void set_peneflux(const double penetrating_muon_flux_in){
      penetrating_muon_flux = penetrating_muon_flux_in; };

    /// @brief set signal and noise (the noise value goes into the Poisson bucket; floor is cleared)
    void set_signal_noise( const double signal_in, const double noise_in ){
      signal = signal_in;
      noise_poi = noise_in;
      noise_det = 0.0;
    };

    /// @brief set signal
    void set_signal( const double signal_in ){ signal = signal_in; };

    /// @brief set total noise (value goes into the Poisson bucket; floor is cleared)
    void set_noise( const double noise_in ){ noise_poi = noise_in; noise_det = 0.0; };

    /// @brief set deterministic floor noise part
    void set_noise_det( const double noise_det_in ){ noise_det = noise_det_in; };

    /// @brief set Poisson-bucket noise part
    void set_noise_poi( const double noise_poi_in ){ noise_poi = noise_poi_in; };

    /// @brief set proj_density
    void set_proj_density(const double proj_density_in)
      { proj_density = proj_density_in; };

    /// @brief set proj_density_lower
    void set_proj_density_lower(const double proj_density_lower_in)
      { proj_density_lower = proj_density_lower_in; };

    /// @brief set proj_density_upper
    void set_proj_density_upper(const double proj_density_upper_in)
    { proj_density_upper = proj_density_upper_in; };

    /// @brief reserve capacity for vec_tf_in_PL
    void reserve_vec_tf_in_PL( const int n_reserve )
    { vec_tf_in_PL.reserve(n_reserve); };

    /// @brief set vec_tf_in_PL
    void set_vec_tf_in_PL( const std::vector<std::tuple<bool,double>> &vec_tf_in_PL_in )
    { vec_tf_in_PL = vec_tf_in_PL_in; };

    /// @brief set eff_low
    void set_eff_low( const double eff_low_in ){ eff_low = eff_low_in; };

    /// @brief set eff_cnt
    void set_eff_cnt( const double eff_cnt_in ){ eff_cnt = eff_cnt_in; };

    /// @brief set eff_upp
    void set_eff_upp( const double eff_upp_in ){ eff_upp = eff_upp_in; };

    /// @brief insert (tf_in, length) to vec_tf_in_PL
    void insert_tf_in_PL( const bool tf_in, const double length )
    { vec_tf_in_PL.push_back(std::make_tuple(tf_in, length)); };

    /// @brief Calculate DL/PL ratio and set as projected density
    /// @param[in] PL_thres Minimum path length (meters) required for valid calculation
    /// @param[in] DL_thres Minimum density length (kg/m²) required for valid calculation
    /// @return Calculated projected density in kg/m³, or invalid_proj_dens (-1.0) if thresholds not met
    /// @details
    /// If PL < PL_thres or DL < DL_thres, sets proj_density to invalid_proj_dens (-1.0).
    /// Otherwise, computes proj_density = DL / PL and stores the result.
    /// @note Units: DL (kg/m²) / PL (m) = kg/m³
    double calc_set_proj_density( const double PL_thres, const double DL_thres );

    ///@} ------------------------------------------------------------------

    //==================================================================
    /// @name calc_functions
    ///@{

    /// @brief Calculate penetrating muon flux using bilinear interpolation
    /// @param[in] g2flux Grid2dXYZ containing pre-computed flux table (cosθ_z vs DL)
    /// @return Interpolated penetrating muon flux in m⁻² sr⁻¹ sec⁻¹
    /// @details
    /// Extracts cosθ_z from ray3d.vz() and DL from the stored density length,
    /// then performs bilinear interpolation in the provided flux grid.
    /// The calculated value is NOT stored in penetrating_muon_flux member.
    /// @note Thread-safety: Read-only operation, safe if g2flux is immutable
    double calc_peneflux( const Grid2dXYZ &g2flux) const;
 
    /// @brief Calculate solid angle of the angular bin
    /// @return Solid angle in steradians (sr)
    /// @details
    /// Dispatches to angle_util::calc_omega_radian or angle_util::calc_omega_tangent based on angle_unit:
    /// - AngleUnit::Degree: converts (txmin, txmax, tymin, tymax) to radians, then calls calc_omega_radian
    /// - AngleUnit::Radian: calls calc_omega_radian directly
    /// - AngleUnit::Tangent: calls calc_omega_tangent
    /// @throws std::runtime_error If angle_unit is not one of the defined enum values
    /// @note The result is always in steradians regardless of the input angle_unit
    double calc_solid_angle() const;

    /// @brief Calculate solid angle using alternative tangent formula
    /// @return Solid angle in steradians (sr)
    /// @details
    /// Similar to calc_solid_angle, but uses angle_util::calc_omega_tangent_alternative
    /// when angle_unit is Tangent. For Degree/Radian, behaves identically to calc_solid_angle.
    /// @throws std::runtime_error If angle_unit is not Degree, Radian, or Tangent
    /// @note Usually not used in production. Provided for validation or alternative calculations.
    double calc_solid_angle_alternative() const;

    /// @brief Calculate effective detection area accounting for angular divergence
    /// @param[in] v3_det_direction Normal vector of the detector surface
    /// @param[in] v3_det_length Physical dimensions (x, y, z) of the detector in meters
    /// @param[in] n_unit Number of detector units (multiplicative factor)
    /// @return Effective area in m²
    /// @details
    /// Computes the projected detector area visible to the muon ray, accounting for:
    /// 1. Angle between detector normal and ray direction (dot product)
    /// 2. Angular divergence (tx, ty) reducing effective x/y dimensions
    /// 3. Multiple detector units (n_unit scaling)
    ///
    /// Formula: effective_area = dot_product × (x_len - z_len×|tan(tx)|) × (y_len - z_len×|tan(ty)|) × n_unit
    ///
    /// The angle_unit determines how tx/ty are interpreted:
    /// - AngleUnit::Degree: converts to radians, then applies tan()
    /// - AngleUnit::Radian: applies tan() directly
    /// - AngleUnit::Tangent: uses tx/ty as tangent values directly
    /// @throws std::runtime_error If v3_det_direction.norm() < 1.0e-6 (zero vector) or angle_unit is undefined
    /// @note Units: v3_det_length in meters → result in m²
    double calc_effective_area(
      const Eigen::Vector3d v3_det_direction
    , const Eigen::Vector3d v3_det_length
    , const double n_unit ) const;

    /// @brief Calculate expected muon signal count
    /// @return Signal count (dimensionless, expected number of muons)
    /// @details
    /// Computes signal = penetrating_muon_flux × effective_area_m2 × solid_angle × exposure_time_sec.
    /// The result is NOT stored in the signal member variable.
    /// @note Units: (m⁻² sr⁻¹ sec⁻¹) × m² × sr × sec = dimensionless count
    double calc_signal() const;

    /// @brief Calculate efficiency-adjusted signal count
    /// @return Efficiency-sampled signal count
    /// @details Multiplies calc_signal() by calc_eff_sample() to incorporate detector efficiency variation
    double calc_signal_eff() const;

    /// @brief Sample efficiency value from asymmetric Gaussian distribution
    /// @return Sampled efficiency value (typically in range [0, 1+])
    /// @details Uses stats_util::sample_asymmetric_gaussian with eff_low, eff_cnt, eff_upp parameters
    double calc_eff_sample() const;

    /// @brief Calculate approximate volume traced by the detector element as a cone
    /// @return Volume in m³
    /// @details
    /// Sums cone volumes from vec_tf_in_PL entries using V = (1/3) × length³ × solid_angle.
    /// Entries with tf_in==true (sky→rock) subtract volume, tf_in==false (rock→sky) add volume.
    /// @throws std::runtime_error If the last tf_in entry is true (should always end outside rock)
    double calc_approx_volume() const;
    ///@} ------------------------------------------------------------------

    //==================================================================
    /// @name checker_functions
    ///@{

    /// @brief Check if detector element position is inside the specified box
    /// @param[in] xmin Minimum x coordinate (meters)
    /// @param[in] xmax Maximum x coordinate (meters)
    /// @param[in] ymin Minimum y coordinate (meters)
    /// @param[in] ymax Maximum y coordinate (meters)
    /// @param[in] zmin Minimum z coordinate (meters)
    /// @param[in] zmax Maximum z coordinate (meters)
    /// @return true if xmin ≤ x < xmax AND ymin ≤ y < ymax AND zmin ≤ z < zmax
    /// @note Uses half-open intervals [min, max) for all coordinates
    bool is_inside(
        const double xmin, const double xmax
      , const double ymin, const double ymax
      , const double zmin, const double zmax ) const;

    /// @brief Check if detector element position is inside the voxel bounds
    /// @param[in] g3vox Grid3dVoxel defining the bounding box
    /// @return true if position is inside the voxel's x/y/z bounds
    bool is_inside( const Grid3dVoxel &g3vox ) const;

    /// @brief Validate that vec_tf_in_PL alternates between true/false and ends with false
    /// @throws std::runtime_error If consecutive entries have the same tf_in value or if the last entry is true
    /// @details
    /// Ensures vec_tf_in_PL represents valid alternating in/out transitions through materials.
    /// The last entry must be false (exiting back to sky).
    void check_alternating_tf_in() const;
    ///@} ------------------------------------------------------------------

    //==================================================================
    /// @name binary_io_functions
    ///@{
    
    /// @brief save all variables to std::ofstream &ofs
    void save( std::ofstream &ofs ) const;

    /// @brief load all variables from std::ifstream &ifs
    void load( std::ifstream &ifs );

    ///@} ------------------------------------------------------------------

};
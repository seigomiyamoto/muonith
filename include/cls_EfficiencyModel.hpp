/// @file cls_EfficiencyModel.hpp
/// @brief Angular efficiency model evaluated at detector-element bin centers
/// @details
/// EfficiencyModel holds the coefficients of a JSON5 efficiency-model file
/// (the eff_table_lab/configs format) and evaluates the efficiency band
/// (eff_low, eff_cnt, eff_upp) at an arbitrary angular coordinate (tx, ty).
/// It replaces the legacy nbinx x nbiny efficiency-table text file: the
/// angular grid is owned by the consumer (DetectorPanel's Grid2dBinGroup),
/// so no grid information is stored here and no grid mismatch can occur.
///
/// Model (identical to scripts/make_eff_table.py, WITHOUT the randomize
/// jitter, which is intentionally not supported on the C++ side):
/// @code
///   eff_cnt   = base_eff * g_tx(tx) * g_ty(ty)            (clipped to [0,1])
///   sigma_eff = sigma.base * 0.5 * (h_tx(tx) + h_ty(ty))  (mean over axes)
///   half      = min(sigma_eff, eff_cnt, 1 - eff_cnt)
///   eff_low   = eff_cnt - half,  eff_upp = eff_cnt + half
/// @endcode
///
/// @note Units: tx and ty are angular coordinates in the unit the detector
///       grid uses (tangent by default); fit functions act on the raw value.
/// @note Thread-safety: eval() is thread-safe (read-only after loading).
#pragma once

#include <filesystem>
#include <string>
#include <vector>

#include <nlohmann/json.hpp>

//###########################################################################
//###########################################################################
/// @class EfficiencyModel
/// @brief Efficiency band (eff_low, eff_cnt, eff_upp) as a function of (tx, ty)
/// @ingroup detectorClasses
//###########################################################################
//###########################################################################
class EfficiencyModel final {
  public:
    //=========================================================================
    /// @name inner types
    ///@{

    /// @brief One per-axis fit function spec ({fn, params} block of the model JSON5)
    /// @details Supported fn values (same as scripts/make_eff_table.py):
    /// flat / poly / gauss / cos_power / sigmoid / step.
    struct AxisFn {
      /// @brief fit function name (flat / poly / gauss / cos_power / sigmoid / step)
      std::string fn = "flat";
      /// @brief poly: coefficients c0, c1, c2, ... of sum_k c_k * t^k
      std::vector<double> coef = {1.0};
      /// @brief gauss: mean of exp(-((t-mu)^2)/(2 sigma^2))
      double mu = 0.0;
      /// @brief gauss: width; must be non-zero
      double sigma = 1.0;
      /// @brief cos_power: exponent of cos(atan(|t|))^n
      double n = 1.0;
      /// @brief cos_power: true when cutoff_deg was given in the JSON
      bool tf_cutoff = false;
      /// @brief cos_power: zenith-angle cutoff in degrees (factor is 0 beyond it)
      double cutoff_deg = 0.0;
      /// @brief sigmoid / step: knee position on |t|
      double x0 = 1.0;
      /// @brief sigmoid / step: transition width; must be non-zero
      double w = 0.1;
      /// @brief step: level for |t| < x0
      double inner = 1.0;
      /// @brief step: level for |t| > x0
      double outer = 0.0;

      /// @brief Read one {fn, params} block from JSON
      /// @param[in] js_axis JSON object holding "fn" and "params"
      /// @param[in] where Label used in error messages (e.g. "tx", "sigma.ty")
      /// @throws std::runtime_error If fn is unknown or a required param is invalid
      void assign_parameters(const nlohmann::json& js_axis, const std::string& where);

      /// @brief Evaluate the axis factor at coordinate t
      /// @param[in] t Angular coordinate (raw value, e.g. tangent)
      /// @return Factor in [0, inf); clipping happens in EfficiencyModel::eval
      double eval(const double t) const;
    };

    /// @brief Efficiency band at one (tx, ty): matches DetectorElement's eff_low/cnt/upp
    struct Band {
      double eff_low = 1.0;
      double eff_cnt = 1.0;
      double eff_upp = 1.0;
    };
    ///@} ---------------------------------------------------------------------

    //=========================================================================
    /// @name constructor & destructor
    ///@{

    /// @brief default constructor (empty model; get_tf_loaded() returns false)
    EfficiencyModel() = default;
    ///@} ---------------------------------------------------------------------

    //=========================================================================
    /// @name loading
    ///@{

    /// @brief Assign model coefficients from an already-parsed JSON object
    /// @param[in] js Top-level JSON object of the model file
    /// @throws std::runtime_error If a "layers" block is present (k-of-n model
    ///         is not supported on the C++ side) or an axis spec is invalid
    /// @note A "randomize" block, if present, is ignored with a warning
    ///       (decided 2026-07-09: no random jitter on the C++ side).
    void assign_parameters(const nlohmann::json& js);

    /// @brief Load the model JSON5 file and assign coefficients
    /// @param[in] path_in Path to the model file (nlohmann-parsable JSON with comments)
    /// @throws std::runtime_error If the file cannot be opened or parsed
    void load(const std::filesystem::path& path_in);
    ///@} ---------------------------------------------------------------------

    //=========================================================================
    /// @name evaluation
    ///@{

    /// @brief Evaluate the efficiency band at one angular coordinate
    /// @param[in] tx Horizontal angular coordinate (bin center)
    /// @param[in] ty Vertical angular coordinate (bin center)
    /// @return Band with 0 <= eff_low <= eff_cnt <= eff_upp <= 1
    /// @note Thread-safety: Yes (read-only).
    Band eval(const double tx, const double ty) const;

    /// @brief true after assign_parameters()/load() succeeded
    bool get_tf_loaded() const { return tf_loaded_; };
    ///@} ---------------------------------------------------------------------

  private:
    /// @brief true after a model has been loaded
    bool tf_loaded_ = false;

    /// @brief central plateau efficiency at g_tx = g_ty = 1
    double base_eff_ = 1.0;

    /// @brief central-efficiency factor along tx
    AxisFn g_tx_ = AxisFn();

    /// @brief central-efficiency factor along ty
    AxisFn g_ty_ = AxisFn();

    /// @brief absolute-uncertainty base (sigma.base); 0 means no band
    double sigma_base_ = 0.0;

    /// @brief uncertainty factor along tx (sigma.tx)
    AxisFn h_tx_ = AxisFn();

    /// @brief uncertainty factor along ty (sigma.ty)
    AxisFn h_ty_ = AxisFn();
};

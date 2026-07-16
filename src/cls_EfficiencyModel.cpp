// src/cls_EfficiencyModel.cpp

#include "cls_EfficiencyModel.hpp"

#include <algorithm>
#include <cmath>

#include "ns_myapp.hpp"
#include "ns_mylogger.hpp"
#include "ns_mymacro.hpp"

// ##############################################
// ##############################################
//  struct EfficiencyModel::AxisFn
// ##############################################
// ##############################################

// Read one {fn, params} block. Semantics mirror scripts/make_eff_table.py.
void EfficiencyModel::AxisFn::assign_parameters(
    const nlohmann::json &js_axis,
    const std::string &where
) {
  fn = js_axis.value("fn", "flat");
  const nlohmann::json js_params = js_axis.value("params", nlohmann::json::object());

  if (fn == "flat") {
    // no params
  } else if (fn == "poly") {
    coef = js_params.value("coef", std::vector<double>{1.0});
  } else if (fn == "gauss") {
    mu = js_params.value("mu", 0.0);
    sigma = js_params.value("sigma", 1.0);
    if (sigma == 0.0) {
      THROW_ERROR("EfficiencyModel::AxisFn::assign_parameters: gauss fn requires non-zero sigma. where={}", where);
    }
  } else if (fn == "cos_power") {
    n = js_params.value("n", 1.0);
    tf_cutoff = js_params.contains("cutoff_deg");
    if (tf_cutoff) { cutoff_deg = js_params.at("cutoff_deg").get<double>(); }
  } else if (fn == "sigmoid" || fn == "step") {
    x0 = js_params.value("x0", 1.0);
    w = js_params.value("w", 0.1);
    if (w == 0.0) {
      THROW_ERROR("EfficiencyModel::AxisFn::assign_parameters: {} fn requires non-zero w. where={}", fn, where);
    }
    if (fn == "step") {
      inner = js_params.value("inner", 1.0);
      outer = js_params.value("outer", 0.0);
    }
  } else {
    THROW_ERROR("EfficiencyModel::AxisFn::assign_parameters: unknown fit fn '{}'. choose from flat/poly/gauss/cos_power/sigmoid/step. where={}", fn, where);
  }
}

// Evaluate the axis factor at t. Same formulas as scripts/make_eff_table.py.
double EfficiencyModel::AxisFn::eval(const double t) const
{
  if (fn == "flat") { return 1.0; }
  if (fn == "poly") {
    double out = 0.0;
    double pw = 1.0;                       // t^k, built incrementally
    for (const double c : coef) { out += c * pw; pw *= t; }
    return out;
  }
  if (fn == "gauss") {
    const double d = t - mu;
    return std::exp(-(d * d) / (2.0 * sigma * sigma));
  }
  if (fn == "cos_power") {
    const double theta = std::atan(std::abs(t));
    if (tf_cutoff && theta > cutoff_deg * M_PI / 180.0) { return 0.0; }
    return std::pow(std::cos(theta), n);
  }
  if (fn == "sigmoid") {
    return 1.0 / (1.0 + std::exp((std::abs(t) - x0) / w));
  }
  if (fn == "step") {
    return outer + (inner - outer) / (1.0 + std::exp((std::abs(t) - x0) / w));
  }
  // assign_parameters rejects unknown fn, so this is unreachable
  THROW_ERROR("EfficiencyModel::AxisFn::eval: unknown fit fn '{}'", fn);
}

// ##############################################
// ##############################################
//  class EfficiencyModel
// ##############################################
// ##############################################

// Assign coefficients from the parsed top-level model JSON object.
void EfficiencyModel::assign_parameters(const nlohmann::json &js)
{
  if (js.contains("layers")) {
    THROW_ERROR("EfficiencyModel::assign_parameters: the 'layers' (k-of-n) model is not supported on the C++ side. key=layers");
  }
  if (js.contains("randomize")) {
    // Decided 2026-07-09: no random jitter on the C++ side; the band comes from sigma only.
    LOG_WARN("EfficiencyModel::assign_parameters: 'randomize' block is ignored (not supported by design)");
  }

  base_eff_ = js.value("base_eff", 1.0);
  if (js.contains("tx")) { g_tx_.assign_parameters(js.at("tx"), "tx"); }
  if (js.contains("ty")) { g_ty_.assign_parameters(js.at("ty"), "ty"); }

  if (js.contains("sigma")) {
    const nlohmann::json &js_sigma = js.at("sigma");
    sigma_base_ = js_sigma.value("base", 0.0);
    if (js_sigma.contains("tx")) { h_tx_.assign_parameters(js_sigma.at("tx"), "sigma.tx"); }
    if (js_sigma.contains("ty")) { h_ty_.assign_parameters(js_sigma.at("ty"), "sigma.ty"); }
  }

  tf_loaded_ = true;
  LOG_INFO("EfficiencyModel::assign_parameters: base_eff={}, tx.fn={}, ty.fn={}, sigma.base={}, sigma.tx.fn={}, sigma.ty.fn={}",
           base_eff_, g_tx_.fn, g_ty_.fn, sigma_base_, h_tx_.fn, h_ty_.fn);
}

// Load the model file (JSON with comments) and assign coefficients.
void EfficiencyModel::load(const std::filesystem::path &path_in)
{
  LOG_INFO("EfficiencyModel::load: path={}", path_in.string());
  assign_parameters(myapp::load_json(path_in));
}

// Evaluate the band at one (tx, ty). Same math as make_eff_table.py build_eff
// (without randomize): center is a product, sigma the mean of its axis factors,
// and the half-width is shrunk so the band stays inside [0, 1].
EfficiencyModel::Band EfficiencyModel::eval(const double tx, const double ty) const
{
  Band band;
  const double eff_raw = base_eff_ * g_tx_.eval(tx) * g_ty_.eval(ty);
  band.eff_cnt = std::clamp(eff_raw, 0.0, 1.0);

  const double sigma_eff = std::max(0.0, sigma_base_ * 0.5 * (h_tx_.eval(tx) + h_ty_.eval(ty)));
  const double half = std::min(sigma_eff, std::min(band.eff_cnt, 1.0 - band.eff_cnt));
  band.eff_low = band.eff_cnt - half;
  band.eff_upp = band.eff_cnt + half;
  return band;
}

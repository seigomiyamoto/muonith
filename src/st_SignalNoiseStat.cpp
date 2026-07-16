// src/st_SignalNoiseStat.cpp
#include "st_SignalNoiseStat.hpp"

#include <map>
#include <fstream>
#include <iostream>
#include <sstream> // istringstream
#include <memory>
#include <string>
#include <cstdio>
#include <cmath>
#include <vector>
#include "ns_angle_util.hpp"
#include "ns_mylogger.hpp"
#include "ns_mymacro.hpp"
#include "ns_stats_util.hpp"

// @brief Output std::vector<SignalNoiseDepth> as ASCII format
void signal_noise_stat_util::out_signal_noise_depth_info_vector(
  const std::vector<SignalNoiseDepth> &vec_snd, FILE *fout )
{
  for (const auto &snd : vec_snd) {
    SignalNoiseSum sn = snd.sn;
    fprintf(fout, "%2d %7.4lf %7.4lf %7.4lf %7.4lf %E %E\n"
      ,snd.detid, sn.tx_lower(), sn.tx_upper(), sn.ty_lower(), sn.ty_upper()
      ,sn.signal, sn.noise);
  }
}

// @brief Output std::vector<SignalNoiseDepth> as ASCII format to specified path
void signal_noise_stat_util::out_signal_noise_depth_info_vector(
    const std::vector<SignalNoiseDepth> &vec_snd
  , const std::filesystem::path &filepath )
{
  FILE *fout = fopen(filepath.string().c_str(), "w");
  if (!fout) THROW_ERROR("signal_noise_stat_util::out_signal_noise_depth_info_vector: Failed to open output file. path={}", filepath.string());
  out_signal_noise_depth_info_vector(vec_snd, fout);
  fclose(fout);
}

/// @brief Evaluate statistical significance of signal differences and store results in vector
void signal_noise_stat_util::eval_signal_significance(
    const std::vector<SignalNoiseDepth> &vec_snd_base
  , const std::vector<SignalNoiseDepth> &vec_snd_modi
  , const double obj_size
  , const double delta_dens
  , const double stat_alpha
  , const double signal_amp
  , std::vector<SignalNoiseStatResult> &vec_stat_result
  , const bool both_side )
{
  if (vec_snd_base.size() != vec_snd_modi.size()) {
    THROW_ERROR("signal_noise_stat_util::eval_signal_significance: vector size mismatch");
  }
  // Reserve memory for vec_stat_result
  const int n_data_new = vec_snd_base.size();
  vec_stat_result.reserve(vec_stat_result.size() + n_data_new);

  // Data loop
  for (size_t i = 0; i < vec_snd_base.size(); ++i) {
    const auto &base = vec_snd_base.at(i);
    const auto &modi = vec_snd_modi.at(i);
    const double signal_base = base.signal() * signal_amp;
    const double signal_modi = modi.signal() * signal_amp;
    const stats_util::StatResult stat
      = stats_util::poisson_count_ztest(
        signal_base, signal_modi, stat_alpha, both_side);
    SignalNoiseStatResult stat_result;
    stat_result.snd_base = base;
    stat_result.snd_modi = modi;
    stat_result.obj_size = obj_size;
    stat_result.delta_dens = delta_dens;
    stat_result.z_value = stat.z_value;
    stat_result.p_value = stat.p_value;
    stat_result.stat_alpha = stat.alpha;
    stat_result.signal_amp = signal_amp;
    stat_result.is_significant = stat.is_significant;

    // Add data
    vec_stat_result.push_back(stat_result);
  } // end of data loop
  LOG_DEBUG_ND("current number of vec_stat_result data: {}", vec_stat_result.size());
}

void signal_noise_stat_util::out_signal_stat_result_vector(
  const std::vector<SignalNoiseStatResult>& vec_stat_result, FILE *fout)
{
  // print header
  fprintf(fout,"# det  txL     txU     tyL     tyU depth_anom_top depth_anom_btm  obsize  ddens   signal_base  signal_modi delta_signal        z_val        p_val alpha signal_amp is_signi\n");
  const int n_data = vec_stat_result.size();
  for (int i = 0; i < n_data; ++i) {
    const auto &res = vec_stat_result.at(i);
    if( res.elev_from_det_lower() < 0.0 ) continue; // skip negative elev_from_det_lower
    fprintf(fout, "%02d %7.4lf %7.4lf %7.4lf %7.4lf"
      , res.detid()
      , res.tx_lower(), res.tx_upper()
      , res.ty_lower(), res.ty_upper() );
    
    fprintf(fout, " %6.1lf %6.1lf %6.1lf %6.0lf "
      , res.elev_from_det_lower(), res.elev_from_det_upper()
      , res.obj_size, res.delta_dens);
    
    const double signal_base = res.signal_base() * res.signal_amp;
    const double signal_modi = res.signal_modi() * res.signal_amp;
    const double delta_signal = signal_modi - signal_base;

    fprintf(fout, " %E %E %E %E %E %E %E %d\n"
      , signal_base, signal_modi, delta_signal
      , res.z_value, res.p_value, res.stat_alpha
      , res.signal_amp
      , res.is_significant ? 1 : 0);
  }
}

void signal_noise_stat_util::out_signal_stat_result_vector(
    const std::vector<SignalNoiseStatResult>& vec_stat_result
  , const std::filesystem::path& outpath )
{
  FILE *fout = fopen(outpath.string().c_str(), "w");
  if (!fout) THROW_ERROR("signal_noise_stat_util::out_signal_stat_result_vector: Failed to open output file. path={}", outpath.string());
  out_signal_stat_result_vector(vec_stat_result, fout);
  fclose(fout);
}

void signal_noise_stat_util::out_signal_stat_result_vector_csv(
  const std::vector<SignalNoiseStatResult>& vec_stat_result, FILE *fout)
{
  if (!fout) fout = stdout;

  // CSV header row
  fprintf(fout,
    "det,txL,txU,tyL,tyU,depth_anom_top,depth_anom_btm,obsize,ddens,"
    "signal_base,signal_modi,delta_signal,"
    "z_val,p_val,alpha,signal_amp,is_signi\n");

  const int n_data = vec_stat_result.size();
  for (int i = 0; i < n_data; ++i) {
    const auto &res = vec_stat_result.at(i);
    if (res.elev_from_det_lower() < 0.0) continue; // skip negative elev_from_det_lower

    const double signal_base = res.signal_base() * res.signal_amp;
    const double signal_modi = res.signal_modi() * res.signal_amp;
    const double delta_signal = signal_modi - signal_base;

    // Data row (comma-separated)
    fprintf(fout,
      "%02d,%.4lf,%.4lf,%.4lf,%.4lf,%.1lf,%.1lf,%.1lf,%.0lf,"
      "%E,%E,%E,%E,%E,%E,%E,%d\n",
      res.detid(),
      res.tx_lower(), res.tx_upper(),
      res.ty_lower(), res.ty_upper(),
      res.elev_from_det_lower(), res.elev_from_det_upper(),
      res.obj_size, res.delta_dens,
      signal_base, signal_modi, delta_signal,
      res.z_value, res.p_value, res.stat_alpha,
      res.signal_amp,
      res.is_significant ? 1 : 0
    );
  }
}

void signal_noise_stat_util::out_signal_stat_result_vector_csv(
    const std::vector<SignalNoiseStatResult>& vec_stat_result
  , const std::filesystem::path& outpath )
{
  FILE *fout = fopen(outpath.string().c_str(), "w");
  if (!fout) THROW_ERROR("signal_noise_stat_util::out_signal_stat_result_vector_csv: Failed to open output file. path={}", outpath.string());
  out_signal_stat_result_vector_csv(vec_stat_result, fout);
  fclose(fout);
}


void signal_noise_stat_util::out_signal_sum_significance_vector(
    const std::vector<SignalNoiseDepth> &vec_snd_base
  , const std::vector<SignalNoiseDepth> &vec_snd_modi
  , const double stat_alpha, FILE *fout )
{
  if (vec_snd_base.size() != vec_snd_modi.size()) {
    THROW_ERROR("signal_noise_stat_util::out_signal_sum_significance_vector: vector size mismatch");
  }

  // print header
  fprintf(fout,"# detid tx_lower tx_upper ty_lower ty_upper "
                "depth_lower depth_upper "
                "signal_sum_base signal_sum_modified "
                "z_value p_value alpha is_significant\n");

  // data
  for (size_t i = 0; i < vec_snd_base.size(); ++i) {
    const auto &base = vec_snd_base.at(i);
    const auto &modi = vec_snd_modi.at(i);
    const stats_util::StatResult stat_result
     = stats_util::poisson_count_ztest(base.signal(), modi.signal(), stat_alpha);
    fprintf(fout, "%2d %7.4lf %7.4lf %7.4lf %7.4lf %7.1lf %7.1lf %E %E %E %E %E %d\n"
      , base.detid, base.tx_lower(), base.tx_upper(), base.ty_lower(), base.ty_upper()
      , base.elev_from_det_lower, base.elev_from_det_upper
      , base.signal(), modi.signal()
      , stat_result.z_value, stat_result.p_value
      , stat_result.alpha, stat_result.is_significant ? 1 : 0);
  }
}

void signal_noise_stat_util::out_signal_sum_significance_vector(
    const std::vector<SignalNoiseDepth> &vec_snd_base
  , const std::vector<SignalNoiseDepth> &vec_snd_modi
  , const double stat_alpha
  , const std::filesystem::path &filepath )
{
  FILE *fout = fopen(filepath.string().c_str(), "w");
  if (!fout) THROW_ERROR("signal_noise_stat_util::out_signal_sum_significance_vector: Failed to open output file. path={}", filepath.string());
  out_signal_sum_significance_vector(vec_snd_base, vec_snd_modi, stat_alpha, fout);
  fclose(fout);
}

// ============================================================
// sort functions
// ============================================================

// ============================================================
// sort function for std::vector<SignalNoiseStatResult>
// ============================================================
void signal_noise_stat_util::sort_signal_stat_results(
    std::vector<SignalNoiseStatResult>& vec_result
  , const std::vector<SignalNoiseStatSortCondition>& vec_condition)
{
  std::sort(vec_result.begin(), vec_result.end(),
    [&](const SignalNoiseStatResult& a, const SignalNoiseStatResult& b) {
      for (const auto& cond : vec_condition) {
        double va = 0.0, vb = 0.0;
        bool bool_a = false, bool_b = false;
        bool is_bool = false;

        switch (cond.key) {
          case SignalNoiseStatSortKey::Detid:
            va = static_cast<double>(a.detid());
            vb = static_cast<double>(b.detid());
            break;

          case SignalNoiseStatSortKey::Tx_lower:
            va = a.tx_lower(); vb = b.tx_lower(); break;
          case SignalNoiseStatSortKey::Tx_upper:
            va = a.tx_upper(); vb = b.tx_upper(); break;
          case SignalNoiseStatSortKey::Ty_lower:
            va = a.ty_lower(); vb = b.ty_lower(); break;
          case SignalNoiseStatSortKey::Ty_upper:
            va = a.ty_upper(); vb = b.ty_upper(); break;

          case SignalNoiseStatSortKey::Depth_lower:
            va = a.elev_from_det_lower(); vb = b.elev_from_det_lower(); break;
          case SignalNoiseStatSortKey::Depth_upper:
            va = a.elev_from_det_upper(); vb = b.elev_from_det_upper(); break;

          case SignalNoiseStatSortKey::SignalBase:
            va = a.signal_base(); vb = b.signal_base(); break;
          case SignalNoiseStatSortKey::SignalModi:
            va = a.signal_modi(); vb = b.signal_modi(); break;

          case SignalNoiseStatSortKey::ObjSize:
            va = a.obj_size; vb = b.obj_size; break;
          case SignalNoiseStatSortKey::DeltaDens:
            va = a.delta_dens; vb = b.delta_dens; break;

          case SignalNoiseStatSortKey::ZValue:
            va = a.z_value; vb = b.z_value; break;
          case SignalNoiseStatSortKey::PValue:
            va = a.p_value; vb = b.p_value; break;

          case SignalNoiseStatSortKey::StatAlpha:
            va = a.stat_alpha; vb = b.stat_alpha; break;

          case SignalNoiseStatSortKey::SignalAmp:
            va = a.signal_amp; vb = b.signal_amp; break;

          case SignalNoiseStatSortKey::IsSignificant:
            bool_a = a.is_significant;
            bool_b = b.is_significant;
            is_bool = true;
            break;
        }

        // Comparison logic
        if (is_bool) {
          if (bool_a != bool_b)
            return cond.ascending ? (bool_a < bool_b) : (bool_a > bool_b);
        } else {
          if (va != vb)
            return cond.ascending ? (va < vb) : (va > vb);
        }
      }
      return false; // All values equal
    });
}

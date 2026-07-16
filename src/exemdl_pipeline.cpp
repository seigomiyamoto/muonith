/// @file exemdl_pipeline.cpp
/// @brief Implementation of pipeline state management.

#include "exemdl_pipeline.hpp"
#include "exemdl_init_json_logger_seed.hpp"
#include "exemdl_build_geometry.hpp"
#include "exemdl_build_detector.hpp"
#include "exemdl_run_inversion.hpp"
#include "exemdl_calc_rec_error.hpp"
#include "cls_MatrixBuildParameters.hpp"
#include "cls_NagaInvParameters.hpp"
#include "cls_Grid2dPillarParameters.hpp"
#include "st_NagaInvLooper.hpp"
#include "ns_pathcalc.hpp"
#include "ns_calc_dNdD.hpp"
#include "ns_io_binary.hpp"
#include "ns_iodir.hpp"
#include "ns_mymacro.hpp"
#include "ns_mylogger.hpp"

namespace exemdl::pipeline {

//----------------------------------------------------------------------
// Internal helper for save(): manifest.json generation
//----------------------------------------------------------------------

namespace {

  /// @brief Column names of mat_element_table_det{NN}.bin, in column order.
  /// @details Shared between write_mat_element_table (value order) and
  /// write_manifest_json (documentation) so the two cannot drift apart.
  /// Angles are tangents; units are encoded in the names.
  constexpr std::array<const char*, 27> element_table_columns = {
    "tan_x_center", "tan_y_center",
    "tan_x_min", "tan_x_max", "tan_y_min", "tan_y_max",
    "path_length_m", "density_length_kg_m2",
    "penetrating_muon_flux", "signal", "noise_det", "noise_poi",
    "effective_area_m2", "solid_angle_sr", "exposure_time_sec",
    "eff_low", "eff_cnt", "eff_upp",
    "proj_density_kg_m3", "proj_density_lower_kg_m3", "proj_density_upper_kg_m3",
    "dir_vx", "dir_vy", "dir_vz",
    "pos_x_m", "pos_y_m", "pos_z_m"
  };

  /// @brief Write a manifest JSON describing the mat/ matrix binaries.
  /// @details
  /// Emits a self-describing companion file so that external readers
  /// (e.g. Python) can interpret mat_dNdD_grouped_{lower,center,upper}.bin
  /// without reading the C++ sources:
  /// - files: binary layout (uint64 rows/cols header, float32 raw data,
  ///   column-major order; see io_binary::write_matxf_stream).
  /// - detectors: per-detector angle_unit (tangent / degree / radian).
  /// - rows: matrix row index -> [detid, igroup, uqig, xmin, xmax, ymin, ymax]
  ///   (analysis-available groups only; the angle range is expressed in the
  ///   owning detector's angle_unit).
  /// - cols: matrix column j -> voxel uqiv = uqiv_min + j -> [uqiv, ix, iy, iz]
  ///   (uqiv numbering of g3vox_merged_input after re_assign_uqiv_by_nhit_det).
  /// @param[in] state Pipeline state (state.det and state.mat must be populated).
  /// @param[in] path Output path of the manifest JSON.
  /// @throws std::runtime_error If the row/column counts derived from the index
  ///         containers do not match the matrix dimensions, or the file cannot
  ///         be opened.
  void write_manifest_json(const State& state, const fs::path& path) {
    const Eigen::MatrixXf& mat_center = state.mat->mat_dNdD_grouped_center;
    const DetectorPanelArray& arrdet = state.det->arrdet_g3vox_input;
    const DetectorIndexContainer& dic = arrdet.get_dic();
    const Grid3dVoxel& g3vox = state.geom.g3vox_merged_input;

    // Row count check: matrix rows = number of analysis-available groups.
    const Index n_rows = static_cast<Index>(dic.get_size_indexer());
    if (n_rows != static_cast<Index>(mat_center.rows())) {
      THROW_ERROR("exemdl::pipeline::write_manifest_json: row count mismatch. indexer={}, matrix_rows={}",
                  n_rows, mat_center.rows());
    }

    // Column count check: matrix column j corresponds to uqiv = uqiv_min + j.
    const Grid3d::Uqiv uqiv_min = g3vox.get_uqiv_min();
    const Grid3d::Uqiv uqiv_max = g3vox.get_uqiv_max();
    const int n_cols = uqiv_max - uqiv_min + 1;
    if (n_cols != static_cast<int>(mat_center.cols())) {
      THROW_ERROR("exemdl::pipeline::write_manifest_json: column count mismatch. uqiv_range={}, matrix_cols={}",
                  n_cols, mat_center.cols());
    }

    nlohmann::json js;
    js["format_version"] = 1;

    // files: the three matrices share one binary layout
    const nlohmann::json js_layout = {
      {"dtype", "float32"},
      {"order", "column-major"},
      {"header", "rows and cols as uint64 x2 (native endianness), followed by raw float32 data"},
      {"rows", mat_center.rows()},
      {"cols", mat_center.cols()}
    };
    js["files"]["mat_dNdD_grouped_lower.bin"] = js_layout;
    js["files"]["mat_dNdD_grouped_center.bin"] = js_layout;
    js["files"]["mat_dNdD_grouped_upper.bin"] = js_layout;

    // detectors: angle_unit is a per-detector property shared by all its rows.
    // Position/direction and the per-detector PL/txty file references are
    // exported for external reconstruction demos (e.g. FBP back projection).
    nlohmann::json js_detectors = nlohmann::json::array();
    const int n_det = arrdet.get_n_det();
    const bool tf_pl_files
      = (static_cast<int>(state.det->vec_spmat_PL.size()) == n_det);
    for (Detid detid = 0; detid < n_det; detid++) {
      const DetectorPanel& panel = arrdet.getDetectorPanel(detid);
      const Eigen::Vector3d v3_pos = panel.get_v3_pos();
      const Eigen::Vector3d v3_dir = panel.get_v3_dir();
      nlohmann::json js_det = {
        {"detid", detid},
        {"angle_unit", angle_util::to_string(panel.get_angle_unit())},
        {"pos_m", {v3_pos.x(), v3_pos.y(), v3_pos.z()}},
        {"dir", {v3_dir.x(), v3_dir.y(), v3_dir.z()}},
        {"n_element", panel.get_n_element()}
      };
      if (tf_pl_files) {
        const std::string name_coo = fmt::format("mat_PL_det{:02d}_coo.bin", detid);
        const std::string name_txty = fmt::format("mat_txty_det{:02d}.bin", detid);
        const SpMatf& sp = state.det->vec_spmat_PL[static_cast<size_t>(detid)];
        js_det["file_PL_coo"] = name_coo;
        js_det["file_txty"] = name_txty;
        js["files"][name_coo] = {
          {"dtype", "float32"},
          {"order", "column-major"},
          {"header", "rows and cols as uint64 x2 (native endianness), followed by raw float32 data"},
          {"rows", sp.nonZeros()},
          {"cols", 3},
          {"columns", {"row (element id_in_this_detector)",
                       "col (voxel index, same ordering as mat_dNdD columns)",
                       "path_length_m"}},
          {"dense_shape", {sp.rows(), sp.cols()}}
        };
        js["files"][name_txty] = {
          {"dtype", "float32"},
          {"order", "column-major"},
          {"header", "rows and cols as uint64 x2 (native endianness), followed by raw float32 data"},
          {"rows", panel.get_n_element()},
          {"cols", 2},
          {"columns", {"tan_x (element center)", "tan_y (element center)"}},
          {"note", "Row i corresponds to COO row index i (element id_in_this_detector). Values are tangents regardless of angle_unit."}
        };
      }
      {
        // Per-element scalar table (written unconditionally with state.det).
        const std::string name_table
          = fmt::format("mat_element_table_det{:02d}.bin", detid);
        js_det["file_element_table"] = name_table;
        nlohmann::json js_cols = nlohmann::json::array();
        for (const char* col : element_table_columns) js_cols.push_back(col);
        js["files"][name_table] = {
          {"dtype", "float32"},
          {"order", "column-major"},
          {"header", "rows and cols as uint64 x2 (native endianness), followed by raw float32 data"},
          {"rows", panel.get_n_element()},
          {"cols", static_cast<int>(element_table_columns.size())},
          {"columns", js_cols},
          {"note", "Row i corresponds to element id_in_this_detector == i (same row order as mat_txty and the COO row indices). path_length_m and density_length_kg_m2 are the raytrace accumulators over the full topography, i.e. they include voxels outside the reconstruction target. Angles are tangents regardless of angle_unit."}
        };
      }
      js_detectors.push_back(js_det);
    }
    js["detectors"] = js_detectors;

    // rows: [detid, igroup, uqig, xmin, xmax, ymin, ymax] per matrix row
    nlohmann::json js_row_entries = nlohmann::json::array();
    for (Index index = 0; index < n_rows; index++) {
      const GroupInfo& gi = dic.getGroupInfo_by_Index(index);
      const std::array<double, 4> xyrange
        = arrdet.getDetectorPanel(gi.detid).get_g2bg().get_xmin_xmax_ymin_ymax(gi.igroup);
      js_row_entries.push_back({gi.detid, gi.igroup, gi.uqig,
                                xyrange[0], xyrange[1], xyrange[2], xyrange[3]});
    }
    js["rows"] = {
      {"meaning", "Matrix row i corresponds to entries[i] (analysis-available detector groups only, dense)."},
      {"columns", {"detid", "igroup", "uqig", "xmin", "xmax", "ymin", "ymax"}},
      {"angle_range_unit", "per-detector angle_unit (see detectors)"},
      {"entries", js_row_entries}
    };

    // cols: [uqiv, ix, iy, iz] per matrix column (uqiv = uqiv_min + j)
    nlohmann::json js_col_entries = nlohmann::json::array();
    for (Grid3d::Uqiv uqiv = uqiv_min; uqiv <= uqiv_max; uqiv++) {
      const Grid3d::Ixiyiz ixiyiz = g3vox.get_ixiyiz(uqiv);
      js_col_entries.push_back({uqiv, ixiyiz[0], ixiyiz[1], ixiyiz[2]});
    }
    js["cols"] = {
      {"meaning", "Matrix column j corresponds to voxel uqiv = uqiv_min + j of g3vox_merged_input."},
      {"columns", {"uqiv", "ix", "iy", "iz"}},
      {"uqiv_min", uqiv_min},
      {"uqiv_max", uqiv_max},
      {"entries", js_col_entries}
    };

    std::ofstream ofs(path);
    if (!ofs) THROW_ERROR("exemdl::pipeline::write_manifest_json: Cannot open '{}'", path.string());
    ofs << js.dump(2);

    LOG_INFO("exemdl::pipeline::write_manifest_json: Wrote '{}' (rows={}, cols={})",
             path.string(), n_rows, n_cols);
  }

  /// @brief Export vec_spmat_PL as per-detector COO triplets plus element
  ///        direction tangents for external reconstruction demos (e.g. FBP).
  /// @details
  /// For each detector detid the following files are written under mat/:
  /// - mat_PL_det{NN}_coo.bin: nnz x 3 float32 matrix. Each row holds
  ///   (row, col, path_length_m) of one nonzero of the sparse path-length
  ///   matrix: row = element id_in_this_detector (same row order as
  ///   mat_txty_det{NN}.bin), col = voxel index in the same ordering as the
  ///   mat_dNdD matrix columns, path length in meters.
  /// - mat_txty_det{NN}.bin: n_element x 2 float32 matrix. Row i holds the
  ///   element-center direction (tan_x, tan_y) of the element whose
  ///   id_in_this_detector == i, always expressed in tangent regardless of
  ///   the detector's angle_unit.
  /// @param[in] state Pipeline state (state.det and state.mat must be populated).
  /// @throws std::runtime_error If vec_spmat_PL does not hold one matrix per
  ///         detector, a per-detector column count differs from the mat_dNdD
  ///         column count, an element index is out of range, or an index
  ///         exceeds the float32 exactly-representable integer range (2^24).
  void write_mat_PL_coo_and_txty(const State& state) {
    const auto& vec_spmat_PL = state.det->vec_spmat_PL;
    const DetectorPanelArray& arrdet = state.det->arrdet_g3vox_input;
    const int n_det = arrdet.get_n_det();
    if (static_cast<int>(vec_spmat_PL.size()) != n_det) {
      THROW_ERROR("exemdl::pipeline::write_mat_PL_coo_and_txty: vec_spmat_PL size mismatch. size={}, n_det={}",
                  vec_spmat_PL.size(), n_det);
    }
    const Eigen::Index n_cols_dndd = state.mat->mat_dNdD_grouped_center.cols();
    constexpr Eigen::Index max_exact_float = static_cast<Eigen::Index>(1) << 24;

    for (Detid detid = 0; detid < n_det; detid++) {
      const SpMatf& sp = vec_spmat_PL[static_cast<size_t>(detid)];
      if (sp.cols() != n_cols_dndd) {
        THROW_ERROR("exemdl::pipeline::write_mat_PL_coo_and_txty: column count mismatch. detid={}, spmat_cols={}, mat_dNdD_cols={}",
                    detid, sp.cols(), n_cols_dndd);
      }
      if (sp.rows() >= max_exact_float || sp.cols() >= max_exact_float) {
        THROW_ERROR("exemdl::pipeline::write_mat_PL_coo_and_txty: index exceeds float32 exact range. detid={}, rows={}, cols={}",
                    detid, sp.rows(), sp.cols());
      }

      Eigen::MatrixXf coo(sp.nonZeros(), 3);
      Eigen::Index k = 0;
      for (int oc = 0; oc < sp.outerSize(); oc++) {
        for (SpMatf::InnerIterator it(sp, oc); it; ++it) {
          coo(k, 0) = static_cast<float>(it.row());
          coo(k, 1) = static_cast<float>(it.col());
          coo(k, 2) = it.value();
          k++;
        }
      }
      io_binary::out_matxf_bin(
        iodir::make_pathout(fmt::format("mat/mat_PL_det{:02d}_coo.bin", detid)), coo);

      const DetectorPanel& panel = arrdet.getDetectorPanel(detid);
      Eigen::MatrixXf txty(panel.get_n_element(), 2);
      for (int iy = 0; iy < panel.get_nbiny(); iy++) {
        for (int ix = 0; ix < panel.get_nbinx(); ix++) {
          const DetectorElement& ele = panel.getDetectorElement(ix, iy);
          const int irow = ele.get_id_in_this_detector();
          if (irow < 0 || irow >= static_cast<int>(txty.rows())) {
            THROW_ERROR("exemdl::pipeline::write_mat_PL_coo_and_txty: element index out of range. detid={}, irow={}, n_element={}",
                        detid, irow, txty.rows());
          }
          const std::array<double, 4> t
            = ele.get_txmin_txmax_tymin_tymax(angle_util::AngleUnit::Tangent);
          txty(irow, 0) = static_cast<float>(0.5 * (t[0] + t[1]));
          txty(irow, 1) = static_cast<float>(0.5 * (t[2] + t[3]));
        }
      }
      io_binary::out_matxf_bin(
        iodir::make_pathout(fmt::format("mat/mat_txty_det{:02d}.bin", detid)), txty);
    }
    LOG_INFO("exemdl::pipeline::write_mat_PL_coo_and_txty: Wrote COO path-length and txty binaries for {} detectors", n_det);
  }

  /// @brief Export one per-element scalar table per detector for external
  ///        analysis (path length, density length, counts, efficiencies).
  /// @details
  /// mat_element_table_det{NN}.bin is an n_element x 27 float32 matrix whose
  /// row i corresponds to the element with id_in_this_detector == i (the same
  /// row order as mat_txty_det{NN}.bin and the COO row indices). The column
  /// order is element_table_columns and is documented in mat/manifest.json.
  /// Unlike the COO export, path_length_m and density_length_kg_m2 come from
  /// the element accumulators (DetectorElement::get_PL/get_DL) and therefore
  /// include voxels outside the reconstruction target (full topography path).
  /// @param[in] state Pipeline state (state.det must be populated).
  /// @throws std::runtime_error If an element index is out of range.
  void write_mat_element_table(const State& state) {
    const DetectorPanelArray& arrdet = state.det->arrdet_g3vox_input;
    const int n_det = arrdet.get_n_det();
    for (Detid detid = 0; detid < n_det; detid++) {
      const DetectorPanel& panel = arrdet.getDetectorPanel(detid);
      Eigen::MatrixXf table(panel.get_n_element(),
                            static_cast<Eigen::Index>(element_table_columns.size()));
      for (int iy = 0; iy < panel.get_nbiny(); iy++) {
        for (int ix = 0; ix < panel.get_nbinx(); ix++) {
          const DetectorElement& ele = panel.getDetectorElement(ix, iy);
          const int irow = ele.get_id_in_this_detector();
          if (irow < 0 || irow >= static_cast<int>(table.rows())) {
            THROW_ERROR("exemdl::pipeline::write_mat_element_table: element index out of range. detid={}, irow={}, n_element={}",
                        detid, irow, table.rows());
          }
          const std::array<double, 4> t
            = ele.get_txmin_txmax_tymin_tymax(angle_util::AngleUnit::Tangent);
          int c = 0;
          table(irow, c++) = static_cast<float>(0.5 * (t[0] + t[1]));
          table(irow, c++) = static_cast<float>(0.5 * (t[2] + t[3]));
          table(irow, c++) = static_cast<float>(t[0]);
          table(irow, c++) = static_cast<float>(t[1]);
          table(irow, c++) = static_cast<float>(t[2]);
          table(irow, c++) = static_cast<float>(t[3]);
          table(irow, c++) = static_cast<float>(ele.get_PL());
          table(irow, c++) = static_cast<float>(ele.get_DL());
          table(irow, c++) = static_cast<float>(ele.get_peneflux());
          table(irow, c++) = static_cast<float>(ele.get_signal());
          table(irow, c++) = static_cast<float>(ele.get_noise_det());
          table(irow, c++) = static_cast<float>(ele.get_noise_poi());
          table(irow, c++) = static_cast<float>(ele.get_effective_area_m2());
          table(irow, c++) = static_cast<float>(ele.get_solid_angle());
          table(irow, c++) = static_cast<float>(ele.get_exposure_time_sec());
          table(irow, c++) = static_cast<float>(ele.get_eff_low());
          table(irow, c++) = static_cast<float>(ele.get_eff_cnt());
          table(irow, c++) = static_cast<float>(ele.get_eff_upp());
          table(irow, c++) = static_cast<float>(ele.get_proj_density());
          table(irow, c++) = static_cast<float>(ele.get_proj_density_lower());
          table(irow, c++) = static_cast<float>(ele.get_proj_density_upper());
          table(irow, c++) = static_cast<float>(ele.get_vx());
          table(irow, c++) = static_cast<float>(ele.get_vy());
          table(irow, c++) = static_cast<float>(ele.get_vz());
          table(irow, c++) = static_cast<float>(ele.get_x());
          table(irow, c++) = static_cast<float>(ele.get_y());
          table(irow, c++) = static_cast<float>(ele.get_z());
        }
      }
      io_binary::out_matxf_bin(
        iodir::make_pathout(fmt::format("mat/mat_element_table_det{:02d}.bin", detid)), table);
    }
    LOG_INFO("exemdl::pipeline::write_mat_element_table: Wrote element tables for {} detectors", n_det);
  }

} // namespace

//----------------------------------------------------------------------
// save() implementation
//----------------------------------------------------------------------

void save(const State& state, const fs::path& dir) {
  iodir::set_default_output_dir(dir);
  LOG_INFO("exemdl::pipeline::save: Saving state to '{}'", dir.string());

  // 1. state_meta.bin
  {
    auto path = iodir::make_pathout("state_meta.bin");
    std::ofstream ofs(path, std::ios::binary);
    if (!ofs) THROW_ERROR("exemdl::pipeline::save: Cannot open '{}'", path.string());

    auto info = io_binary::get_current_architecture_info();
    io_binary::write_architecture_info(ofs, info);
    io_binary::write_binary(ofs, PIPELINE_VERSION);
    io_binary::write_binary(ofs, state.completed_module);
    io_binary::write_path(ofs, state.path_json);
    io_binary::write_binary(ofs, state.det.has_value());
    io_binary::write_binary(ofs, state.mat.has_value());
    io_binary::write_binary(ofs, state.geom.has_shell_upper);
    io_binary::write_binary(ofs, state.geom.has_shell_lower);
    io_binary::write_binary(ofs, state.geom.has_shell_lateral);
    io_binary::write_binary(ofs, true); // has_shell_pl (always true for version 3)
  }

  // 2. app_params.json (human-readable)
  {
    auto path = iodir::make_pathout("app_params.json");
    std::ofstream ofs(path);
    if (!ofs) THROW_ERROR("exemdl::pipeline::save: Cannot open '{}'", path.string());
    ofs << state.js.dump(2);
  }

  // 3. geo/
  state.geom.g2pil_naive.save(iodir::make_pathout("geom/g2pil_naive.bin"));
  state.geom.g3vox_input.save(iodir::make_pathout("geom/g3vox_input.bin"));
  state.geom.g3vox_merged_input.save(iodir::make_pathout("geom/g3vox_merged_input.bin"));
  if (state.geom.has_shell_upper) {
    state.geom.g2pil_shell_upper.save(iodir::make_pathout("geom/g2pil_shell_upper.bin"));
  }
  if (state.geom.has_shell_lower) {
    state.geom.g2pil_shell_lower.save(iodir::make_pathout("geom/g2pil_shell_lower.bin"));
  }
  if (state.geom.has_shell_lateral) {
    state.geom.g2pil_shell_lateral.save(iodir::make_pathout("geom/g2pil_shell_lateral.bin"));
  }

  // 4. det/ (optional)
  if (state.det) {
    state.det->arrdet_g2pil_naive.save(iodir::make_pathout("det/arrdet_g2pil_naive.bin"));
    state.det->arrdet_g3vox_input.save(iodir::make_pathout("det/arrdet_g3vox_input.bin"));
    io_binary::write_vec_spmatf(iodir::make_pathout("det/vec_spmat_PL.bin"), state.det->vec_spmat_PL);
    state.det->shell_pl.save(dir / "det");
    {
      auto path = iodir::make_pathout("det/prior_info_all.bin");
      std::ofstream ofs(path, std::ios::binary);
      if (!ofs) THROW_ERROR("exemdl::pipeline::save: Cannot open '{}'", path.string());
      state.det->prior_info_all.save(ofs);
    }
    {
      auto path = iodir::make_pathout("det/avr_dens.bin");
      std::ofstream ofs(path, std::ios::binary);
      if (!ofs) THROW_ERROR("exemdl::pipeline::save: Cannot open '{}'", path.string());
      io_binary::write_binary(ofs, state.det->avr_dens_lower_center_upper[0]);
      io_binary::write_binary(ofs, state.det->avr_dens_lower_center_upper[1]);
      io_binary::write_binary(ofs, state.det->avr_dens_lower_center_upper[2]);
      io_binary::write_binary(ofs, state.det->density_quad[0]);
      io_binary::write_binary(ofs, state.det->density_quad[1]);
      io_binary::write_binary(ofs, state.det->density_quad[2]);
      io_binary::write_binary(ofs, state.det->density_quad[3]);
    }
  }

  // 5. mat/ (optional)
  if (state.mat) {
    io_binary::out_matxf_bin(iodir::make_pathout("mat/mat_dNdD_grouped_lower.bin"),
                             state.mat->mat_dNdD_grouped_lower);
    io_binary::out_matxf_bin(iodir::make_pathout("mat/mat_dNdD_grouped_center.bin"),
                             state.mat->mat_dNdD_grouped_center);
    io_binary::out_matxf_bin(iodir::make_pathout("mat/mat_dNdD_grouped_upper.bin"),
                             state.mat->mat_dNdD_grouped_upper);
    if (state.det) {
      write_mat_element_table(state);
      write_mat_PL_coo_and_txty(state);
      write_manifest_json(state, iodir::make_pathout("mat/manifest.json"));
    }
  }

  LOG_INFO("exemdl::pipeline::save: Done (completed_module={})", state.completed_module);
}

//----------------------------------------------------------------------
// load() implementation
//----------------------------------------------------------------------

State load(const fs::path& dir) {
  LOG_INFO("exemdl::pipeline::load: Loading state from '{}'", dir.string());
  State state;

  bool has_det = false;
  bool has_mat = false;
  bool has_shell_pl = false;

  // 1. state_meta.bin
  {
    auto path = dir / "state_meta.bin";
    std::ifstream ifs(path, std::ios::binary);
    if (!ifs) THROW_ERROR("exemdl::pipeline::load: Cannot open '{}'", path.string());

    auto file_info = io_binary::read_architecture_info(ifs);
    io_binary::check_architecture_compatibility_or_throw(file_info);

    uint32_t version = io_binary::read_binary<uint32_t>(ifs);
    if (version != PIPELINE_VERSION) {
      THROW_ERROR("exemdl::pipeline::load: Version mismatch (file={}, expected={})",
                  version, PIPELINE_VERSION);
    }

    state.completed_module = io_binary::read_binary<int>(ifs);
    state.path_json = io_binary::read_path(ifs);
    has_det = io_binary::read_binary<bool>(ifs);
    has_mat = io_binary::read_binary<bool>(ifs);
    state.geom.has_shell_upper = io_binary::read_binary<bool>(ifs);
    state.geom.has_shell_lower = io_binary::read_binary<bool>(ifs);
    state.geom.has_shell_lateral = io_binary::read_binary<bool>(ifs);
    has_shell_pl = io_binary::read_binary<bool>(ifs);
  }

  // 2. app_params.json
  {
    auto path = dir / "app_params.json";
    std::ifstream ifs(path);
    if (!ifs) THROW_ERROR("exemdl::pipeline::load: Cannot open '{}'", path.string());
    ifs >> state.js;
  }

  // 3. geo/
  // Checkpoint binaries store canonical (already shifted) coordinates and are
  // loaded pass-through; no tf_shift interpretation is needed on resume.
  state.geom.g2pil_naive.load(dir / "geom/g2pil_naive.bin");
  state.geom.g3vox_input.load(dir / "geom/g3vox_input.bin");
  state.geom.g3vox_merged_input.load(dir / "geom/g3vox_merged_input.bin");
  if (state.geom.has_shell_upper) {
    state.geom.g2pil_shell_upper.load(dir / "geom/g2pil_shell_upper.bin");
  }
  if (state.geom.has_shell_lower) {
    state.geom.g2pil_shell_lower.load(dir / "geom/g2pil_shell_lower.bin");
  }
  if (state.geom.has_shell_lateral) {
    state.geom.g2pil_shell_lateral.load(dir / "geom/g2pil_shell_lateral.bin");
  }

  // 4. det/ (optional)
  if (has_det) {
    state.det = DetState{};
    state.det->arrdet_g2pil_naive.load(dir / "det/arrdet_g2pil_naive.bin");
    state.det->arrdet_g3vox_input.load(dir / "det/arrdet_g3vox_input.bin");
    state.det->vec_spmat_PL = io_binary::read_vec_spmatf(dir / "det/vec_spmat_PL.bin");
    if (has_shell_pl) {
      state.det->shell_pl = exemdl::build_prior::ShellPL::load(dir / "det");
    }
    {
      auto path = dir / "det/prior_info_all.bin";
      std::ifstream ifs(path, std::ios::binary);
      if (!ifs) THROW_ERROR("exemdl::pipeline::load: Cannot open '{}'", path.string());
      state.det->prior_info_all.load(ifs);
    }
    {
      auto path = dir / "det/avr_dens.bin";
      std::ifstream ifs(path, std::ios::binary);
      if (!ifs) THROW_ERROR("exemdl::pipeline::load: Cannot open '{}'", path.string());
      state.det->avr_dens_lower_center_upper[0] = io_binary::read_binary<double>(ifs);
      state.det->avr_dens_lower_center_upper[1] = io_binary::read_binary<double>(ifs);
      state.det->avr_dens_lower_center_upper[2] = io_binary::read_binary<double>(ifs);
      state.det->density_quad[0] = io_binary::read_binary<double>(ifs);
      state.det->density_quad[1] = io_binary::read_binary<double>(ifs);
      state.det->density_quad[2] = io_binary::read_binary<double>(ifs);
      state.det->density_quad[3] = io_binary::read_binary<double>(ifs);
    }
  }

  // 5. mat/ (optional)
  if (has_mat) {
    state.mat = MatState{};
    state.mat->mat_dNdD_grouped_lower = io_binary::read_matxf_bin(dir / "mat/mat_dNdD_grouped_lower.bin");
    state.mat->mat_dNdD_grouped_center = io_binary::read_matxf_bin(dir / "mat/mat_dNdD_grouped_center.bin");
    state.mat->mat_dNdD_grouped_upper = io_binary::read_matxf_bin(dir / "mat/mat_dNdD_grouped_upper.bin");
  }

  LOG_INFO("exemdl::pipeline::load: Done (completed_module={})", state.completed_module);
  return state;
}

//----------------------------------------------------------------------
// init() implementation
//----------------------------------------------------------------------

State init(const fs::path& json_path, bool seed_given, unsigned seed_value) {
  LOG_INFO("exemdl::pipeline::init: Starting initialization with '{}'", json_path.string());
  State state;
  state.path_json = json_path;

  // init: JSON loading, logger setup, seed initialization, thread configuration
  exemdl::init_json_logger_seed::InitArgs init_args{ json_path, seed_given, seed_value };
  state.js = exemdl::init_json_logger_seed::init_all(init_args);

  // load_parameters: Load application parameters from JSON
  state.app_params = exemdl::load_parameters::load_all(state.js);

  // build_geometry: Build geometry (terrain, voxel grids, shells)
  exemdl::build_geometry::BuildArgs geom_args{ state.app_params };
  auto geom_result = exemdl::build_geometry::build_all(geom_args);
  state.geom.g2pil_naive        = std::move(geom_result.g2pil_naive);
  state.geom.g3vox_input        = std::move(geom_result.g3vox_input);
  state.geom.g3vox_merged_input = std::move(geom_result.g3vox_merged_input);
  state.geom.g2pil_shell_upper   = std::move(geom_result.g2pil_shell_upper);
  state.geom.g2pil_shell_lower   = std::move(geom_result.g2pil_shell_lower);
  state.geom.g2pil_shell_lateral = std::move(geom_result.g2pil_shell_lateral);
  state.geom.has_shell_upper     = geom_result.has_shell_upper;
  state.geom.has_shell_lower     = geom_result.has_shell_lower;
  state.geom.has_shell_lateral   = geom_result.has_shell_lateral;

  state.completed_module = 3;
  LOG_INFO("exemdl::pipeline::init: Completed init/load_parameters/build_geometry (Module 1-3)");
  return state;
}

//----------------------------------------------------------------------
// run_trace_path_lengths() implementation — density-independent, expensive
//----------------------------------------------------------------------

void run_trace_path_lengths(State& state) {
  if (state.completed_module < 3) {
    THROW_ERROR("exemdl::pipeline::run_trace_path_lengths: Build Geometry (Module 3) not completed. completed_module={}",
                state.completed_module);
  }
  LOG_INFO("exemdl::pipeline::run_trace_path_lengths: Starting (density-independent raytrace)");

  // Initialize DetState
  state.det.emplace();

  // 1. Create arrdet_template from detector parameters
  DetectorPanelArray arrdet_template(state.app_params.prm_det);
  arrdet_template.set_name("arrdet_template");

  // 2. Build arrdet_g2pil_naive (2D terrain-based detector)
  const nlohmann::json& js_proj_dens = state.js.at("PROJ_DENS_EVAL_GROUPED");
  state.det->arrdet_g2pil_naive = exemdl::build_detector::build_arrdet_g2pil(
      state.app_params
    , state.geom.g2pil_naive
    , arrdet_template
    , js_proj_dens
    , state.app_params.ft_prior
  );

  // 3. Build arrdet_g3vox_input, vec_spmat_PL, and vecxf_non_rec_vox_PL
  auto [arrdet_g3vox, vec_spmat_PL, vecxf_non_rec_vox_PL]
  = exemdl::build_detector::build_arrdet_g3vox(
        state.app_params
      , state.geom.g3vox_merged_input
      , state.geom.g2pil_shell_upper
      , state.geom.g2pil_shell_lower
      , state.geom.g2pil_shell_lateral
      , state.geom.has_shell_upper
      , state.geom.has_shell_lower
      , state.geom.has_shell_lateral
      , state.det->arrdet_g2pil_naive  // copy bimap from this
      , arrdet_template
    );
  state.det->arrdet_g3vox_input = std::move(arrdet_g3vox);
  state.det->vec_spmat_PL       = std::move(vec_spmat_PL);

  if (js_proj_dens.value("tf_exec", false)) {
    LOG_INFO("Computing projected density for arrdet_g3vox_input before g2bg output");
    // Efficiency-uncertainty flags: same gate as the C_N diagonal (tf_eff_cn_diag),
    // so the projected-density error band widens together with the forward count.
    const bool tf_eff =
        !state.app_params.vec_prm_nagainv.empty()
     && state.app_params.vec_prm_nagainv.front().get_tf_eff_cn_diag();
    const bool tf_eff_independent =
        tf_eff
     && state.app_params.vec_prm_nagainv.front().get_tf_eff_cn_diag_independent();
    state.det->arrdet_g3vox_input.mp_calc_set_proj_dens_grouped_all(
      js_proj_dens, state.app_params.ft_prior, tf_eff, tf_eff_independent);
  } else {
    LOG_INFO("PROJ_DENS_EVAL_GROUPED.tf_exec=false; outputting arrdet_g3vox_input g2bg without projected density");
  }
  if (state.app_params.prm_det.tf_out_g2bg_ascii) {
    LOG_INFO("Outputting arrdet_g3vox_input g2bg after projected density step");
    state.det->arrdet_g3vox_input.out_g2bg_all();
  } else {
    LOG_INFO("Skipped g2bg text output (DETECTOR_PARAMETER_LISTS.tf_out_g2bg_ascii=false)");
  }

  // 4. Build shell_pl (density-independent shell path lengths)
  state.det->shell_pl = exemdl::build_prior::build_shell_PL(
      state.app_params
    , state.det->arrdet_g3vox_input
    , state.geom.g2pil_shell_upper
    , state.geom.g2pil_shell_lower
    , state.geom.g2pil_shell_lateral
    , state.geom.has_shell_upper
    , state.geom.has_shell_lower
    , state.geom.has_shell_lateral
  );
  // attach filtered-out voxel PL captured in step 3
  state.det->shell_pl.vecxf_non_rec_vox_PL = std::move(vecxf_non_rec_vox_PL);

  // Do NOT set completed_module here; compute_prior will do that.
  LOG_INFO("exemdl::pipeline::run_trace_path_lengths: Completed (raytrace done, prior pending)");
}

//----------------------------------------------------------------------
// run_compute_prior() implementation — density-dependent, cheap
//----------------------------------------------------------------------

void run_compute_prior(State& state, const std::array<double, 4>& density_quad) {
  if (!state.det) {
    THROW_ERROR("exemdl::pipeline::run_compute_prior: DetState is not initialized. "
                "run_trace_path_lengths must be called first.");
  }
  const double voxel_density       = density_quad[0];
  const double shell_density_upper   = density_quad[1];
  const double shell_density_lower   = density_quad[2];
  const double shell_density_lateral = density_quad[3];
  LOG_INFO("exemdl::pipeline::run_compute_prior: Starting with density_quad="
           "[prior={}, upper={}, lower={}, lateral={}]",
           voxel_density, shell_density_upper, shell_density_lower, shell_density_lateral);

  // Build prior_info_all via cached shell PL (rebuild_prior_from_shell_PL).
  // shell_pl.vecxf_non_rec_vox_PL covers PL through n_hit_det-filtered voxels.
  std::array<double, 3> avr_dens = { voxel_density, voxel_density, voxel_density };

  const bool tf_prior_error = state.app_params.tf_run_inversion_prior_error;
  if (tf_prior_error) {
    // Build all three (lower/center/upper)
    LOG_INFO("tf_prior_error=true: rebuilding lower/center/upper priors from cached shell PL");
    state.det->prior_info_all.lower = exemdl::build_prior::rebuild_prior_from_shell_PL(
        "_lower", state.app_params,
        state.det->shell_pl, state.det->vec_spmat_PL,
        state.geom.g3vox_merged_input, avr_dens[0],
        shell_density_upper, shell_density_lower, shell_density_lateral,
        state.det->arrdet_g3vox_input, state.det->arrdet_g2pil_naive);

    state.det->prior_info_all.center = exemdl::build_prior::rebuild_prior_from_shell_PL(
        "_center", state.app_params,
        state.det->shell_pl, state.det->vec_spmat_PL,
        state.geom.g3vox_merged_input, avr_dens[1],
        shell_density_upper, shell_density_lower, shell_density_lateral,
        state.det->arrdet_g3vox_input, state.det->arrdet_g2pil_naive);

    state.det->prior_info_all.upper = exemdl::build_prior::rebuild_prior_from_shell_PL(
        "_upper", state.app_params,
        state.det->shell_pl, state.det->vec_spmat_PL,
        state.geom.g3vox_merged_input, avr_dens[2],
        shell_density_upper, shell_density_lower, shell_density_lateral,
        state.det->arrdet_g3vox_input, state.det->arrdet_g2pil_naive);
  } else {
    // Build center only, copy to lower/upper
    LOG_INFO("tf_prior_error=false: rebuilding center prior only from cached shell PL");
    auto center_prior = exemdl::build_prior::rebuild_prior_from_shell_PL(
        "_center", state.app_params,
        state.det->shell_pl, state.det->vec_spmat_PL,
        state.geom.g3vox_merged_input, voxel_density,
        shell_density_upper, shell_density_lower, shell_density_lateral,
        state.det->arrdet_g3vox_input, state.det->arrdet_g2pil_naive);
    state.det->prior_info_all.center = std::move(center_prior);
    state.det->prior_info_all.lower  = state.det->prior_info_all.center;  // copy
    state.det->prior_info_all.upper  = state.det->prior_info_all.center;  // copy
  }

  state.det->avr_dens_lower_center_upper = avr_dens;
  state.det->density_quad                = density_quad;

  state.completed_module = 5;
  LOG_INFO("exemdl::pipeline::run_compute_prior: Completed Compute Prior (Module 5)");
}

//----------------------------------------------------------------------
// run_build_observation_matrix() implementation - internal helper
//----------------------------------------------------------------------

namespace {
  /// @brief Internal helper to build a single observation matrix.
  Eigen::MatrixXf build_single_matrix(
      State& state,
      const std::vector<Eigen::VectorXf>& vec_vecxf_DL)
  {
    pathcalc::MatrixBuildParameters prm_mat(
        "prm_pathmatrix"
      , state.app_params.prm_path.tf_load_bin_obs_mat_dNdD
      , state.app_params.prm_path.tf_save_bin_obs_mat_dNdD
      , state.app_params.prm_path.path_bin_obs_mat_dNdD
      , state.det->arrdet_g3vox_input
      , state.geom.g3vox_merged_input
      , state.det->vec_spmat_PL
      , vec_vecxf_DL
      , state.app_params.ft_prior
      , state.app_params.prm_det.tf_apply_eff
    );

    return calc_dNdD::create_grouped_mat_dNdD_alldet_sprs(prm_mat);
  }

  /// @brief Overlay reconstruction densities onto the pre-merge grid and output cross sections.
  /// @param[in] g3vox_premerge Pre-merge full-extent grid (state.geom.g3vox_input).
  /// @param[in] g3vox_merged_input  Full reconst_volume merged grid.
  /// @param[in] allRes        Inversion results containing merged-grid reconstructions.
  /// @param[in] params        Application parameters (for zcross config).
  /// @param[in] density_quad  Shell density [prior, upper, lower, lateral].
  void output_overlay_cross_sections(
      const Grid3dVoxel& g3vox_premerge,
      const Grid3dVoxel& g3vox_merged_input,
      const exemdl::run_inversion::InversionResultsAll& allRes,
      const exemdl::load_parameters::AppParameters& params,
      const std::array<double,4>& density_quad = {})
  {
    // Build CrossSectionZParameters from the pre-merge grid
    Grid3dVoxel::CrossSectionZParameters prm;
    prm.xmin  = g3vox_premerge.get_xmin();
    prm.xmax  = g3vox_premerge.get_xmax();
    prm.xstep = g3vox_premerge.get_x_interval();
    prm.ymin  = g3vox_premerge.get_ymin();
    prm.ymax  = g3vox_premerge.get_ymax();
    prm.ystep = g3vox_premerge.get_y_interval();
    prm.zmin  = params.zcross.min;
    prm.zmax  = params.zcross.max;
    const double z_interval = g3vox_premerge.get_z_interval();
    prm.zstep = (params.zcross.zstep > 0.0) ? params.zcross.zstep : z_interval;
    prm.n_detector    = params.prm_det.get_n_det();
    prm.output_binary = params.zcross.output_binary;

    // Helper lambda: overlay and write cross sections for one result
    auto overlay_and_write = [&](const exemdl::run_inversion::InversionResults& inv,
                                 const std::string& label) {
      LOG_INFO("Overlaying merged density onto pre-merge grid ({})", label);
      Grid3dVoxel g3vox_overlay = g3vox_premerge.overlay_merged_density(
          g3vox_merged_input, inv.g3vox_rec, density_quad[0]);
      g3vox_overlay.apply_shell_density(density_quad);
      const std::string name = fmt::format("g3vox_overlay_rec_{}_zcross_all.tmp", label);
      const fs::path pathout = iodir::make_pathout(name);
      g3vox_overlay.out_cross_section_z_all(pathout, prm);
    };

    // Always output center
    overlay_and_write(allRes.center, "center");

    // Output lower/upper only when independent reconstructions exist
    if (allRes.has_lower) {
      overlay_and_write(allRes.lower, "lower");
    }
    if (allRes.has_upper) {
      overlay_and_write(allRes.upper, "upper");
    }
  }
} // anonymous namespace

//----------------------------------------------------------------------
// run_build_observation_matrix() implementation - using prior_info_all's DL vectors
//----------------------------------------------------------------------

void run_build_observation_matrix(State& state) {
  if (state.completed_module < 5) {
    THROW_ERROR("exemdl::pipeline::run_build_observation_matrix: Compute Prior (Module 5) not completed. completed_module={}",
                state.completed_module);
  }
  if (!state.det) {
    THROW_ERROR("exemdl::pipeline::run_build_observation_matrix: DetState is not initialized");
  }
  LOG_INFO("exemdl::pipeline::run_build_observation_matrix: Starting (using prior_info_all DL vectors)");

  // Initialize MatState
  state.mat.emplace();

  const bool tf_prior_error = state.app_params.tf_run_inversion_prior_error;
  if (tf_prior_error) {
    // Build all three matrices
    LOG_INFO("tf_prior_error=true: building lower/center/upper matrices");

    LOG_INFO("Building observation matrix for lower-bound prior");
    state.mat->mat_dNdD_grouped_lower = build_single_matrix(
        state, state.det->prior_info_all.lower.vec_vecxf_DL);

    LOG_INFO("Building observation matrix for center prior");
    state.mat->mat_dNdD_grouped_center = build_single_matrix(
        state, state.det->prior_info_all.center.vec_vecxf_DL);

    LOG_INFO("Building observation matrix for upper-bound prior");
    state.mat->mat_dNdD_grouped_upper = build_single_matrix(
        state, state.det->prior_info_all.upper.vec_vecxf_DL);
  } else {
    // Build center only, copy to lower/upper
    LOG_INFO("tf_prior_error=false: building center matrix only");

    LOG_INFO("Building observation matrix for center prior");
    state.mat->mat_dNdD_grouped_center = build_single_matrix(
        state, state.det->prior_info_all.center.vec_vecxf_DL);

    state.mat->mat_dNdD_grouped_lower = state.mat->mat_dNdD_grouped_center;
    state.mat->mat_dNdD_grouped_upper = state.mat->mat_dNdD_grouped_center;
  }

  state.completed_module = 6;
  LOG_INFO("exemdl::pipeline::run_build_observation_matrix: Completed Build Observation Matrix (Module 6)");
}

//----------------------------------------------------------------------
// run_build_observation_matrix() implementation - with specified DL vectors
//----------------------------------------------------------------------

void run_build_observation_matrix(State& state, const std::vector<Eigen::VectorXf>& vec_vecxf_DL) {
  if (state.completed_module < 5) {
    THROW_ERROR("exemdl::pipeline::run_build_observation_matrix: Compute Prior (Module 5) not completed. completed_module={}",
                state.completed_module);
  }
  if (!state.det) {
    THROW_ERROR("exemdl::pipeline::run_build_observation_matrix: DetState is not initialized");
  }
  LOG_INFO("exemdl::pipeline::run_build_observation_matrix: Starting (using provided DL vectors)");

  // Initialize MatState
  state.mat.emplace();

  // Build all three matrices using the same DL vectors (for iterative estimation)
  LOG_INFO("Building observation matrix (same DL for lower/center/upper)");
  state.mat->mat_dNdD_grouped_lower  = build_single_matrix(state, vec_vecxf_DL);
  state.mat->mat_dNdD_grouped_center = build_single_matrix(state, vec_vecxf_DL);
  state.mat->mat_dNdD_grouped_upper  = build_single_matrix(state, vec_vecxf_DL);

  state.completed_module = 6;
  LOG_INFO("exemdl::pipeline::run_build_observation_matrix: Completed Build Observation Matrix (Module 6)");
}

//----------------------------------------------------------------------
// calc_vec_vecxf_DL() implementation
//----------------------------------------------------------------------

std::vector<Eigen::VectorXf> calc_vec_vecxf_DL(
    const State& state, const Eigen::VectorXf& vecxf_density)
{
  if (!state.det) {
    THROW_ERROR("exemdl::pipeline::calc_vec_vecxf_DL: DetState is not initialized");
  }

  const auto& vec_spmat_PL = state.det->vec_spmat_PL;
  std::vector<Eigen::VectorXf> vec_vecxf_DL;
  vec_vecxf_DL.reserve(vec_spmat_PL.size());

  for (const auto& spmat_PL : vec_spmat_PL) {
    // DL = PL * density
    // spmat_PL: (n_ele, n_vox_avail)
    // vecxf_density: (n_vox_avail)
    // result: (n_ele)
    Eigen::VectorXf vecxf_DL = spmat_PL * vecxf_density;
    vec_vecxf_DL.push_back(std::move(vecxf_DL));
  }

  return vec_vecxf_DL;
}

//----------------------------------------------------------------------
// run_invert_density() implementation
//----------------------------------------------------------------------

exemdl::run_inversion::InversionResultsAll run_invert_density(
    const State& state, double sigma_rho, double corr_length)
{
  // 1. Precondition checks
  if (state.completed_module < 6) {
    THROW_ERROR("exemdl::pipeline::run_invert_density: Build Observation Matrix (Module 6) not completed. completed_module={}",
                state.completed_module);
  }
  if (!state.det) {
    THROW_ERROR("exemdl::pipeline::run_invert_density: DetState is not initialized");
  }
  if (!state.mat) {
    THROW_ERROR("exemdl::pipeline::run_invert_density: MatState is not initialized");
  }
  LOG_INFO("exemdl::pipeline::run_invert_density: Starting with sigma_rho={}, corr_length={}",
           sigma_rho, corr_length);

  // 2. Create a copy of AppParameters (manual clone due to NagaInvParameters deleted copy ctor)
  if (state.app_params.vec_prm_nagainv.empty()) {
    THROW_ERROR("exemdl::pipeline::run_invert_density: vec_prm_nagainv is empty");
  }

  exemdl::load_parameters::AppParameters params_copy;
  params_copy.zcross      = state.app_params.zcross;
  params_copy.prm_det     = state.app_params.prm_det;
  params_copy.prm_g2pil   = state.app_params.prm_g2pil;
  params_copy.prm_bingroup = state.app_params.prm_bingroup;
  params_copy.prm_g3vox   = state.app_params.prm_g3vox;
  params_copy.prm_g3merge = state.app_params.prm_g3merge;
  params_copy.prm_path    = state.app_params.prm_path;
  params_copy.prm_noise   = state.app_params.prm_noise;
  params_copy.ft_real     = state.app_params.ft_real;
  params_copy.ft_prior    = state.app_params.ft_prior;
  params_copy.tf_run_inversion_prior_error = state.app_params.tf_run_inversion_prior_error;

  // 3. Clone and modify NagaInvParameters with specified sigma_rho and corr_length
  NagaInvParameters prm_modified = state.app_params.vec_prm_nagainv.front().clone();
  prm_modified.json_params["sigma_rho"] = sigma_rho;
  prm_modified.json_params["corr_length"] = corr_length;
  params_copy.vec_prm_nagainv.push_back(std::move(prm_modified));

  // 4. Build calc_matrix::BuildResult from MatState
  exemdl::calc_matrix::BuildResult res_mat;
  res_mat.mat_dNdD_grouped_lower  = state.mat->mat_dNdD_grouped_lower;
  res_mat.mat_dNdD_grouped_center = state.mat->mat_dNdD_grouped_center;
  res_mat.mat_dNdD_grouped_upper  = state.mat->mat_dNdD_grouped_upper;

  // 5. Get observed muon counts (signal Poisson switch; noise always det floor + Poisson(poi))
  const bool tf_signal_poisson = params_copy.vec_prm_nagainv.front().get_tf_signal_poisson();
  const Eigen::VectorXf vecxf_nmuon_obs =
      state.det->arrdet_g3vox_input.get_vecxf_nmuon_poisson_all(tf_signal_poisson);

  // 5b. Efficiency uncertainty as analytic variance on the C_N diagonal (opt-in via flag).
  //     Built from the same detector array and AvailIndexer order as vecxf_nmuon_obs,
  //     so the resulting vector is row-aligned with the C_N diagonal.
  Eigen::VectorXf vecxf_var_eff; // empty -> disabled (backward compatible)
  if (params_copy.vec_prm_nagainv.front().get_tf_eff_cn_diag()) {
    const bool tf_eff_independent =
      params_copy.vec_prm_nagainv.front().get_tf_eff_cn_diag_independent();
    vecxf_var_eff =
      state.det->arrdet_g3vox_input.get_vecxf_var_eff_all(tf_eff_independent);
  }

  // 6. Build inversion arguments
  exemdl::run_inversion::BuildArgs inv_args{
    .index_run       = 0,
    .prior_info_all  = state.det->prior_info_all,
    .g3vox_input     = state.geom.g3vox_merged_input,
    .res_mat         = res_mat,
    .vecxf_nmuon_obs = vecxf_nmuon_obs,
    .tf_use_dens_input = true,
    .params          = params_copy,
    .density_quad    = state.det->density_quad,
    .vecxf_var_eff   = vecxf_var_eff
  };

  // 7. Execute inversion
  auto result = exemdl::run_inversion::build_all(inv_args);

  // 8. Overlay reconstruction densities onto pre-merge grid and output cross sections
  output_overlay_cross_sections(state.geom.g3vox_input, state.geom.g3vox_merged_input, result, params_copy, state.det->density_quad);

  LOG_INFO("exemdl::pipeline::run_invert_density: Completed");
  return result;
}

//----------------------------------------------------------------------
exemdl::run_inversion::InversionResultsAll run_invert_density(
    const State& state, double sigma_rho, double corr_length, double sigma_rho_diag)
{
  // 1. Precondition checks
  if (state.completed_module < 6) {
    THROW_ERROR("exemdl::pipeline::run_invert_density: Build Observation Matrix (Module 6) not completed. completed_module={}",
                state.completed_module);
  }
  if (!state.det) {
    THROW_ERROR("exemdl::pipeline::run_invert_density: DetState is not initialized");
  }
  if (!state.mat) {
    THROW_ERROR("exemdl::pipeline::run_invert_density: MatState is not initialized");
  }
  LOG_INFO("exemdl::pipeline::run_invert_density: Starting with sigma_rho={}, corr_length={}, sigma_rho_diag={}",
           sigma_rho, corr_length, sigma_rho_diag);

  // 2. Create a copy of AppParameters (manual clone due to NagaInvParameters deleted copy ctor)
  if (state.app_params.vec_prm_nagainv.empty()) {
    THROW_ERROR("exemdl::pipeline::run_invert_density: vec_prm_nagainv is empty");
  }

  exemdl::load_parameters::AppParameters params_copy;
  params_copy.zcross      = state.app_params.zcross;
  params_copy.prm_det     = state.app_params.prm_det;
  params_copy.prm_g2pil   = state.app_params.prm_g2pil;
  params_copy.prm_bingroup = state.app_params.prm_bingroup;
  params_copy.prm_g3vox   = state.app_params.prm_g3vox;
  params_copy.prm_g3merge = state.app_params.prm_g3merge;
  params_copy.prm_path    = state.app_params.prm_path;
  params_copy.prm_noise   = state.app_params.prm_noise;
  params_copy.ft_real     = state.app_params.ft_real;
  params_copy.ft_prior    = state.app_params.ft_prior;
  params_copy.tf_run_inversion_prior_error = state.app_params.tf_run_inversion_prior_error;

  // 3. Clone and modify NagaInvParameters
  NagaInvParameters prm_modified = state.app_params.vec_prm_nagainv.front().clone();
  prm_modified.json_params["sigma_rho"] = sigma_rho;
  prm_modified.json_params["corr_length"] = corr_length;
  prm_modified.json_params["sigma_rho_diag"] = sigma_rho_diag;
  params_copy.vec_prm_nagainv.push_back(std::move(prm_modified));

  // 4. Build calc_matrix::BuildResult from MatState
  exemdl::calc_matrix::BuildResult res_mat;
  res_mat.mat_dNdD_grouped_lower  = state.mat->mat_dNdD_grouped_lower;
  res_mat.mat_dNdD_grouped_center = state.mat->mat_dNdD_grouped_center;
  res_mat.mat_dNdD_grouped_upper  = state.mat->mat_dNdD_grouped_upper;

  // 5. Get observed muon counts (signal Poisson switch; noise always det floor + Poisson(poi))
  const bool tf_signal_poisson = params_copy.vec_prm_nagainv.front().get_tf_signal_poisson();
  const Eigen::VectorXf vecxf_nmuon_obs =
      state.det->arrdet_g3vox_input.get_vecxf_nmuon_poisson_all(tf_signal_poisson);

  // 5b. Efficiency uncertainty as analytic variance on the C_N diagonal (opt-in via flag).
  //     Built from the same detector array and AvailIndexer order as vecxf_nmuon_obs,
  //     so the resulting vector is row-aligned with the C_N diagonal.
  Eigen::VectorXf vecxf_var_eff; // empty -> disabled (backward compatible)
  if (params_copy.vec_prm_nagainv.front().get_tf_eff_cn_diag()) {
    const bool tf_eff_independent =
      params_copy.vec_prm_nagainv.front().get_tf_eff_cn_diag_independent();
    vecxf_var_eff =
      state.det->arrdet_g3vox_input.get_vecxf_var_eff_all(tf_eff_independent);
  }

  // 6. Build inversion arguments
  exemdl::run_inversion::BuildArgs inv_args{
    .index_run       = 0,
    .prior_info_all  = state.det->prior_info_all,
    .g3vox_input     = state.geom.g3vox_merged_input,
    .res_mat         = res_mat,
    .vecxf_nmuon_obs = vecxf_nmuon_obs,
    .tf_use_dens_input = true,
    .params          = params_copy,
    .density_quad    = state.det->density_quad,
    .vecxf_var_eff   = vecxf_var_eff
  };

  // 7. Execute inversion
  auto result = exemdl::run_inversion::build_all(inv_args);

  // 8. Overlay reconstruction densities onto pre-merge grid and output cross sections
  output_overlay_cross_sections(state.geom.g3vox_input, state.geom.g3vox_merged_input, result, params_copy, state.det->density_quad);

  LOG_INFO("exemdl::pipeline::run_invert_density: Completed");
  return result;
}

//----------------------------------------------------------------------
// save_recon_io() implementation
//----------------------------------------------------------------------

void save_recon_io(
    const State& state,
    const exemdl::run_inversion::InversionResultsAll& inv_result,
    double sigma_rho, double corr_length, double sigma_rho_diag)
{
  // 1. Precondition checks
  if (!state.det) {
    THROW_ERROR("exemdl::pipeline::save_recon_io: DetState is not initialized");
  }
  if (!state.mat) {
    THROW_ERROR("exemdl::pipeline::save_recon_io: MatState is not initialized");
  }
  if (state.app_params.vec_prm_nagainv.empty()) {
    THROW_ERROR("exemdl::pipeline::save_recon_io: vec_prm_nagainv is empty");
  }
  const auto& prm = state.app_params.vec_prm_nagainv.front();
  const Grid3dVoxel& g3vox = state.geom.g3vox_merged_input;
  const Eigen::MatrixXf& mat_center = state.mat->mat_dNdD_grouped_center;
  const Index n_rows = static_cast<Index>(mat_center.rows());
  const int n_cols = static_cast<int>(mat_center.cols());

  // 2. Observed muon counts actually used by run_invert_density. The Poisson values
  //    are precomputed reads stored in the detector panels, so this call
  //    reproduces the identical vector passed to the inversion.
  const bool tf_signal_poisson = prm.get_tf_signal_poisson();
  const Eigen::VectorXf vecxf_nmuon_obs =
      state.det->arrdet_g3vox_input.get_vecxf_nmuon_poisson_all(tf_signal_poisson);

  // 3. Row/column alignment checks (row order = AvailIndexer = matrix rows;
  //    column order = uqiv = matrix columns).
  auto check_rows = [&](const Eigen::VectorXf& v, const char* name) {
    if (v.size() != n_rows) {
      THROW_ERROR("exemdl::pipeline::save_recon_io: {} size mismatch. size={}, matrix_rows={}",
                  name, v.size(), n_rows);
    }
  };
  auto check_cols = [&](const Eigen::VectorXf& v, const char* name) {
    if (v.size() != n_cols) {
      THROW_ERROR("exemdl::pipeline::save_recon_io: {} size mismatch. size={}, matrix_cols={}",
                  name, v.size(), n_cols);
    }
  };

  nlohmann::json js;
  js["format_version"] = 1;
  js["density_unit"] = "kg/m^3";
  js["length_unit"] = "m";
  js["row_order"] = "identical to ../mat/manifest.json rows.entries (matrix row i)";
  js["col_order"] = "identical to ../mat/manifest.json cols.entries (voxel uqiv = uqiv_min + j)";

  // 4. Write one header+float32 column-major .bin and register it in the manifest.
  auto write_rec_bin = [&](const std::string& name, const Eigen::MatrixXf& m) {
    const fs::path path = iodir::make_pathout("rec/" + name);
    io_binary::out_matxf_bin(path, m);
    js["files"][name] = {
      {"dtype", "float32"},
      {"order", "column-major"},
      {"header", "rows and cols as uint64 x2 (native endianness), followed by raw float32 data"},
      {"rows", m.rows()},
      {"cols", m.cols()}
    };
  };

  check_rows(vecxf_nmuon_obs, "N_obs");
  write_rec_bin("vec_nmuon_obs.bin", vecxf_nmuon_obs);

  const auto& pia = state.det->prior_info_all;
  check_rows(pia.lower.vecxf_nmuon, "N_prior_lower");
  check_rows(pia.center.vecxf_nmuon, "N_prior_center");
  check_rows(pia.upper.vecxf_nmuon, "N_prior_upper");
  write_rec_bin("vec_nmuon_prior_lower.bin", pia.lower.vecxf_nmuon);
  write_rec_bin("vec_nmuon_prior_center.bin", pia.center.vecxf_nmuon);
  write_rec_bin("vec_nmuon_prior_upper.bin", pia.upper.vecxf_nmuon);

  const Eigen::VectorXf vecxf_rho_prior_lower  = pia.lower.g3vox.get_vecxf_density();
  const Eigen::VectorXf vecxf_rho_prior_center = pia.center.g3vox.get_vecxf_density();
  const Eigen::VectorXf vecxf_rho_prior_upper  = pia.upper.g3vox.get_vecxf_density();
  check_cols(vecxf_rho_prior_lower, "rho_prior_lower");
  check_cols(vecxf_rho_prior_center, "rho_prior_center");
  check_cols(vecxf_rho_prior_upper, "rho_prior_upper");
  write_rec_bin("vec_rho_prior_lower.bin", vecxf_rho_prior_lower);
  write_rec_bin("vec_rho_prior_center.bin", vecxf_rho_prior_center);
  write_rec_bin("vec_rho_prior_upper.bin", vecxf_rho_prior_upper);

  const Eigen::VectorXf vecxf_rho_true = g3vox.get_vecxf_density();
  check_cols(vecxf_rho_true, "rho_true");
  write_rec_bin("vec_rho_true.bin", vecxf_rho_true);

  const Eigen::VectorXf vecxf_rho_rec_center = inv_result.center.g3vox_rec.get_vecxf_density();
  const Eigen::VectorXf vecxf_dens_err_center = inv_result.center.g3vox_dens_err.get_vecxf_density();
  check_cols(vecxf_rho_rec_center, "rho_rec_center");
  check_cols(vecxf_dens_err_center, "dens_err_center");
  write_rec_bin("vec_rho_rec_center.bin", vecxf_rho_rec_center);
  write_rec_bin("vec_dens_err_center.bin", vecxf_dens_err_center);
  if (inv_result.has_lower) {
    const Eigen::VectorXf v = inv_result.lower.g3vox_rec.get_vecxf_density();
    check_cols(v, "rho_rec_lower");
    write_rec_bin("vec_rho_rec_lower.bin", v);
  }
  if (inv_result.has_upper) {
    const Eigen::VectorXf v = inv_result.upper.g3vox_rec.get_vecxf_density();
    check_cols(v, "rho_rec_upper");
    write_rec_bin("vec_rho_rec_upper.bin", v);
  }

  // 5. Voxel center coordinates (n_voxel x 3, meters), same column order.
  const Grid3d::Uqiv uqiv_min = g3vox.get_uqiv_min();
  Eigen::MatrixXf mat_voxel_xyz(n_cols, 3);
  for (int j = 0; j < n_cols; j++) {
    const Grid3d::Ixiyiz ixiyiz = g3vox.get_ixiyiz(uqiv_min + j);
    const AABB3d aabb = std::get<1>(g3vox.get_density_AABB3d(ixiyiz));
    const Eigen::Vector3d v3_center = aabb.center();
    mat_voxel_xyz(j, 0) = static_cast<float>(v3_center.x());
    mat_voxel_xyz(j, 1) = static_cast<float>(v3_center.y());
    mat_voxel_xyz(j, 2) = static_cast<float>(v3_center.z());
  }
  write_rec_bin("mat_voxel_xyz.bin", mat_voxel_xyz);

  // 6. Solver settings and result summary. sigma_rho / corr_length /
  //    sigma_rho_diag are the values actually passed to run_invert_density;
  //    the remaining settings come from the base NagaInvParameters.
  js["solver"] = {
    {"formula", "rho' = rho_0 + (A^T C_N^-1 A + C_rho^-1)^-1 A^T C_N^-1 (N_obs - N_prior)"},
    {"cov_muon", "diagonal; C_N(i,i) = |N_obs(i)|, replaced by nmuon_under_thres when below nmuon_thres"},
    {"cov_dens", "C_rho(j,k) = sigma_rho^2 * exp(-d_jk / corr_length) (isotropic case); diagonal = sigma_rho_diag^2"},
    {"sigma_rho", sigma_rho},
    {"corr_length", corr_length},
    {"sigma_rho_diag", sigma_rho_diag},
    {"sigma_rho_diag_note", "negative means Module 7 fell back to sigma_rho"},
    {"nmuon_thres", prm.get_nmuon_thres()},
    {"nmuon_under_thres", prm.get_nmuon_under_thres()},
    {"tf_signal_poisson", tf_signal_poisson},
    {"tf_eff_cn_diag", prm.get_tf_eff_cn_diag()},
    {"tf_aniso", prm.is_anisotropic()},
    {"aniso_cov_type", prm.get_aniso_cov_type()},
    {"corr_length_xy", prm.get_corr_length_xy()},
    {"corr_length_z", prm.get_corr_length_z()}
  };
  js["results"] = {
    {"chi2_muon", inv_result.center.chi2_muon},
    {"ndf_muon", inv_result.center.ndf_muon},
    {"has_lower", inv_result.has_lower},
    {"has_upper", inv_result.has_upper}
  };

  const fs::path manifest_path = iodir::make_pathout("rec/manifest.json");
  std::ofstream ofs(manifest_path);
  if (!ofs) {
    THROW_ERROR("exemdl::pipeline::save_recon_io: Cannot open '{}'", manifest_path.string());
  }
  ofs << js.dump(1) << "\n";
  LOG_INFO("exemdl::pipeline::save_recon_io: Done ({} files + manifest under '{}')",
           js["files"].size(), manifest_path.parent_path().string());
}

//----------------------------------------------------------------------
void run_analyze_errors(
    const State& state,
    const exemdl::run_inversion::InversionResultsAll& inv_result,
    const std::string& output_prefix)
{
  if (!inv_result.center.pNagainv) {
    THROW_ERROR("exemdl::pipeline::run_analyze_errors: pNagainv is null");
  }
  if (!state.det) {
    THROW_ERROR("exemdl::pipeline::run_analyze_errors: DetState is not initialized");
  }

  const std::string name = output_prefix.empty()
    ? "nagalooper00"
    : fmt::format("nagalooper_{}", output_prefix);

  // Clone NagaInvParameters to keep alive for looper (takes const ref)
  NagaInvParameters prm_clone = state.app_params.vec_prm_nagainv.front().clone();

  NagaInvLooper looper{
      name,
      *inv_result.center.pNagainv,
      prm_clone,
      state.geom.g3vox_merged_input,
      inv_result.center.prm_zcross,
      state.det->arrdet_g3vox_input.get_dic()
  };

  const Eigen::VectorXf vecxf_dens_prior_cnt =
    state.det->prior_info_all.center.g3vox.get_vecxf_density();
  const Eigen::VectorXf& vecxf_dens_rec =
    inv_result.center.reconst_res.vecxf_dens_rec;

  LOG_INFO("exemdl::pipeline::run_analyze_errors: exec_disable_mode_all");
  auto vec_reconst_res_disabled =
    looper.exec_disable_mode_all(vecxf_dens_prior_cnt, vecxf_dens_rec);

  if (vec_reconst_res_disabled.empty()) {
    THROW_ERROR("exemdl::pipeline::run_analyze_errors: No result from exec_disable_mode_all");
  }

  // Use PriorInfoAll overload to avoid copying DetectorPanelArray (deleted copy ctor)
  exemdl::calc_rec_error::write_disable_and_errors(
    looper, vec_reconst_res_disabled, inv_result,
    state.det->prior_info_all,
    state.geom.g3vox_merged_input, output_prefix);

  LOG_INFO("exemdl::pipeline::run_analyze_errors: Completed");
}

} // namespace exemdl::pipeline

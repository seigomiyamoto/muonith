// cls_Grid2dBinGroupParameters.cpp

#include "cls_Grid2dBinGroup.hpp"
#include "cls_Grid2dBinGroupParameters.hpp"
#include "ns_myapp.hpp"
#include "ns_mymacro.hpp"
#include "ns_mylogger.hpp"
#include "spdlog_pch.hpp"
#include "ns_mylogger.hpp"

//###############################################
// class Grid2dBinGroup::Parameters
//###############################################

// Inequality operator overload
bool Grid2dBinGroup::Parameters::operator!=(const Grid2dBinGroup::Parameters& other) const
{
  #ifdef NODEBUG
    if( signal_init != other.signal_init ) return true;
    if( noise_init != other.noise_init ) return true;
    if( is_avail_init != other.is_avail_init ) return true;
    if( PL_thres != other.PL_thres ) return true;
    if( DL_thres != other.DL_thres ) return true;
    if( signal_under_thres != other.signal_under_thres ) return true;
    if( noise_under_thres != other.noise_under_thres ) return true;
    if( is_avail_under_thres != other.is_avail_under_thres ) return true;
    if( tf_run_1st_grouping != other.tf_run_1st_grouping ) return true;
    if( tf_run_auto_grouping != other.tf_run_auto_grouping ) return true;
    if( igroup_start != other.igroup_start ) return true;
    if( nx_div_init != other.nx_div_init ) return true;
    if( ny_div_init != other.ny_div_init ) return true;
    if( signal_noise_group_trig != other.signal_noise_group_trig ) return true;
    if( ixlen_min != other.ixlen_min ) return true;
    if( iylen_min != other.iylen_min ) return true;
    if( tf_prefer_split_x != other.tf_prefer_split_x ) return true;
    if( nloop_limit != other.nloop_limit ) return true;
    if( n_detector_grouping_manual != other.n_detector_grouping_manual ) return true;
    if( vec_tf_read_bin_group_list != other.vec_tf_read_bin_group_list ) return true;
    if( vec_file_path_bin_group_list != other.vec_file_path_bin_group_list ) return true;
  #else
    if( signal_init != other.signal_init ) { LOG_WARN("signal_init differs"); return true; }
    if( noise_init != other.noise_init ) { LOG_WARN("noise_init differs"); return true; }
    if( is_avail_init != other.is_avail_init ) { LOG_WARN("is_avail_init differs"); return true; }
    if( PL_thres != other.PL_thres ) { LOG_WARN("PL_thres differs"); return true; }
    if( DL_thres != other.DL_thres ) { LOG_WARN("DL_thres differs"); return true; }
    if( signal_under_thres != other.signal_under_thres ) { LOG_WARN("signal_under_thres differs"); return true; }
    if( noise_under_thres != other.noise_under_thres ) { LOG_WARN("noise_under_thres differs"); return true; }
    if( is_avail_under_thres != other.is_avail_under_thres ) { LOG_WARN("is_avail_under_thres differs"); return true; }
    if( tf_run_1st_grouping != other.tf_run_1st_grouping ) { LOG_WARN("tf_run_1st_grouping differs"); return true; }
    if( tf_run_auto_grouping != other.tf_run_auto_grouping ) { LOG_WARN("tf_run_auto_grouping differs"); return true; }
    if( igroup_start != other.igroup_start ) { LOG_WARN("igroup_start differs"); return true; }
    if( nx_div_init != other.nx_div_init ) { LOG_WARN("nx_div_init differs"); return true; }
    if( ny_div_init != other.ny_div_init ) { LOG_WARN("ny_div_init differs"); return true; }
    if( signal_noise_group_trig != other.signal_noise_group_trig ) { LOG_WARN("signal_noise_group_trig differs"); return true; }
    if( ixlen_min != other.ixlen_min ) { LOG_WARN("ixlen_min differs"); return true; }
    if( iylen_min != other.iylen_min ) { LOG_WARN("iylen_min differs"); return true; }
    if( tf_prefer_split_x != other.tf_prefer_split_x ) { LOG_WARN("tf_prefer_split_x differs"); return true; }
    if( nloop_limit != other.nloop_limit ) { LOG_WARN("nloop_limit differs"); return true; }
    if( n_detector_grouping_manual != other.n_detector_grouping_manual ) { LOG_WARN("n_detector_grouping_manual differs"); return true; }
    if( vec_tf_read_bin_group_list != other.vec_tf_read_bin_group_list ) { LOG_WARN("vec_tf_read_bin_group_list differs"); return true; }
    if( vec_file_path_bin_group_list != other.vec_file_path_bin_group_list ) { LOG_WARN("vec_file_path_bin_group_list differs"); return true; }
  #endif
  return false;
}

void Grid2dBinGroup::Parameters::assign_parameters(const nlohmann::json& js, const std::string &section_name)
{
  LOG_INFO("...");
  std::string key;

  // name
  key = TOSTRING( name );
  if(js.contains(section_name) && js.at(section_name).contains(key))
    name = js.at(section_name).at(key).get<std::string>();

  // signal_init
  key = TOSTRING( signal_init );
  if(js.contains(section_name) && js.at(section_name).contains(key))
    signal_init = js.at(section_name).at(key).get<double>();

  // noise_init
  key = TOSTRING( noise_init );
  if(js.contains(section_name) && js.at(section_name).contains(key))
    noise_init = js.at(section_name).at(key).get<double>();

  // is_avail_init
  key = TOSTRING( is_avail_init );
  if(js.contains(section_name) && js.at(section_name).contains(key))
    is_avail_init = js.at(section_name).at(key).get<bool>();

  // tf_run_1st_grouping
  key = TOSTRING( tf_run_1st_grouping );
  if(js.contains(section_name) && js.at(section_name).contains(key))
    tf_run_1st_grouping = js.at(section_name).at(key).get<bool>();

  // tf_run_auto_grouping
  key = TOSTRING( tf_run_auto_grouping );
  if(js.contains(section_name) && js.at(section_name).contains(key))
    tf_run_auto_grouping = js.at(section_name).at(key).get<bool>();

  // PL_thres
  key = TOSTRING( PL_thres );
  if(js.contains(section_name) && js.at(section_name).contains(key))
    PL_thres = js.at(section_name).at(key).get<double>();
  assert( PL_thres >= 0.0 );

  // DL_thres
  key = TOSTRING( DL_thres );
  if(js.contains(section_name) && js.at(section_name).contains(key))
    DL_thres = js.at(section_name).at(key).get<double>();
  assert( DL_thres >= 0.0 );

  // signal_under_thres
  key = TOSTRING( signal_under_thres );
  if(js.contains(section_name) && js.at(section_name).contains(key))
    signal_under_thres = js.at(section_name).at(key).get<double>();

  // noise_under_thres
  key = TOSTRING( noise_under_thres );
  if(js.contains(section_name) && js.at(section_name).contains(key))
    noise_under_thres = js.at(section_name).at(key).get<double>();

  // is_avail_under_thres
  key = TOSTRING( is_avail_under_thres );
  if(js.contains(section_name) && js.at(section_name).contains(key))
    is_avail_under_thres = js.at(section_name).at(key).get<bool>();

  // igroup_start
  key = TOSTRING( igroup_start );
  if(js.contains(section_name) && js.at(section_name).contains(key))
    igroup_start = js.at(section_name).at(key).get<int>();
  assert( igroup_start >= 0 );

  // nx_div_init
  key = TOSTRING( nx_div_init );
  if(js.contains(section_name) && js.at(section_name).contains(key))
    nx_div_init = js.at(section_name).at(key).get<int>();
  assert( nx_div_init >= 1 );

  // ny_div_init
  key = TOSTRING( ny_div_init );
  if(js.contains(section_name) && js.at(section_name).contains(key))
    ny_div_init = js.at(section_name).at(key).get<int>();
  assert( ny_div_init >= 1 );

  // signal_noise_group_trig
  key = TOSTRING( signal_noise_group_trig );
  if(js.contains(section_name) && js.at(section_name).contains(key))
    signal_noise_group_trig = js.at(section_name).at(key).get<double>();

  // ixlen_min
  key = TOSTRING( ixlen_min );
  if(js.contains(section_name) && js.at(section_name).contains(key))
    ixlen_min = js.at(section_name).at(key).get<int>();
  assert( ixlen_min >= 1 );

  // iylen_min
  key = TOSTRING( iylen_min );
  if(js.contains(section_name) && js.at(section_name).contains(key))
    iylen_min = js.at(section_name).at(key).get<int>();
  assert( iylen_min >= 1 );

  // tf_prefer_split_x
  key = TOSTRING( tf_prefer_split_x );
  if(js.contains(section_name) && js.at(section_name).contains(key))
    tf_prefer_split_x = js.at(section_name).at(key).get<bool>();

  // nloop_limit
  key = TOSTRING( nloop_limit );
  if(js.contains(section_name) && js.at(section_name).contains(key))
    nloop_limit = js.at(section_name).at(key).get<int>();
  assert( nloop_limit >= 1 );

  // If both tf_run_1st_grouping and tf_run_auto_grouping are false
  if( !tf_run_1st_grouping && !tf_run_auto_grouping ){
    // Read n_detector_grouping_manual
    key = TOSTRING( n_detector_grouping_manual );
    if(js.contains(section_name) && js.at(section_name).contains(key))
      n_detector_grouping_manual = js.at(section_name).at(key).get<int>();
    assert( n_detector_grouping_manual >= 0 );

    if( n_detector_grouping_manual > 0 ){
      // Read vec_tf_read_bin_group_list directly from JSON array
      key = "vec_tf_read_bin_group_list";
      if(js.contains(section_name) && js.at(section_name).contains(key))
        vec_tf_read_bin_group_list = js.at(section_name).at(key).get<std::vector<bool>>();

      // Read vec_file_path_bin_group_list directly from JSON array
      key = "vec_file_path_bin_group_list";
      if(js.contains(section_name) && js.at(section_name).contains(key)){
        std::vector<std::string> temp = js.at(section_name).at(key).get<std::vector<std::string>>();
        vec_file_path_bin_group_list.clear();
        for(const auto &str : temp)
          vec_file_path_bin_group_list.push_back(std::filesystem::path(str));
      }
    }
  }
  LOG_INFO("...OK");
}

void Grid2dBinGroup::Parameters::save( std::ofstream& ofs ) const
{
  io_binary::write_string( ofs, name ); // save name
  io_binary::write_binary( ofs, signal_init ); // save signal_init
  io_binary::write_binary( ofs, noise_init ); // save noise_init
  io_binary::write_bool( ofs, is_avail_init ); // save is_avail_init
  io_binary::write_binary( ofs, PL_thres ); // save PL_thres
  io_binary::write_binary( ofs, DL_thres ); // save DL_thres
  io_binary::write_binary( ofs, signal_under_thres ); // save signal_under_thres
  io_binary::write_binary( ofs, noise_under_thres ); // save noise_under_thres
  io_binary::write_bool( ofs, is_avail_under_thres ); // save is_avail_under_thres
  io_binary::write_bool( ofs, tf_run_1st_grouping ); // save tf_run_1st_grouping
  io_binary::write_bool( ofs, tf_run_auto_grouping ); // save tf_run_auto_grouping
  io_binary::write_binary( ofs, igroup_start ); // save igroup_start
  io_binary::write_binary( ofs, nx_div_init ); // save nx_div_init
  io_binary::write_binary( ofs, ny_div_init ); // save ny_div_init
  io_binary::write_binary( ofs, signal_noise_group_trig ); // save signal_noise_group_trig
  io_binary::write_binary( ofs, ixlen_min ); // save ixlen_min
  io_binary::write_binary( ofs, iylen_min ); // save iylen_min
  io_binary::write_bool( ofs, tf_prefer_split_x ); // save tf_prefer_split_x
  io_binary::write_binary( ofs, nloop_limit ); // save nloop_limit
  io_binary::write_binary( ofs, n_detector_grouping_manual ); // save n_detector_grouping_manual
  io_binary::write_vec_bool( ofs, vec_tf_read_bin_group_list ); // save vec_tf_read_bin_group_list
  io_binary::write_vec( ofs, vec_file_path_bin_group_list ); // save vec_file_path_bin_group_list
  if( ofs.fail() ) THROW_ERROR("Grid2dBinGroup::Parameters::save: Output stream failed.");
}

void Grid2dBinGroup::Parameters::load( std::ifstream& ifs )
{
  name = io_binary::read_string( ifs ); // load name
  signal_init = io_binary::read_binary<double>( ifs ); // load signal_init
  noise_init = io_binary::read_binary<double>( ifs ); // load noise_init
  is_avail_init = io_binary::read_bool( ifs ); // load is_avail_init
  PL_thres = io_binary::read_binary<double>( ifs ); // load PL_thres
  DL_thres = io_binary::read_binary<double>( ifs ); // load DL_thres
  signal_under_thres = io_binary::read_binary<double>( ifs ); // load signal_under_thres
  noise_under_thres = io_binary::read_binary<double>( ifs ); // load noise_under_thres
  is_avail_under_thres = io_binary::read_bool( ifs ); // load is_avail_under_thres
  tf_run_1st_grouping = io_binary::read_bool( ifs ); // load tf_run_1st_grouping
  tf_run_auto_grouping = io_binary::read_bool( ifs ); // load tf_run_auto_grouping
  igroup_start = io_binary::read_binary<int>( ifs ); // load igroup_start
  nx_div_init = io_binary::read_binary<int>( ifs ); // load nx_div_init
  ny_div_init = io_binary::read_binary<int>( ifs ); // load ny_div_init
  signal_noise_group_trig = io_binary::read_binary<double>( ifs ); // load signal_noise_group_trig
  ixlen_min = io_binary::read_binary<int>( ifs ); // load ixlen_min
  iylen_min = io_binary::read_binary<int>( ifs ); // load iylen_min
  tf_prefer_split_x = io_binary::read_bool( ifs ); // load tf_prefer_split_x
  nloop_limit = io_binary::read_binary<int>( ifs ); // load nloop_limit
  n_detector_grouping_manual = io_binary::read_binary<int>( ifs ); // load n_detector_grouping_manual
  vec_tf_read_bin_group_list = io_binary::read_vec_bool( ifs ); // load vec_tf_read_bin_group_list
  vec_file_path_bin_group_list = io_binary::read_vec<std::filesystem::path>( ifs ); // load vec_file_path_bin_group_list
  if( ifs.fail() ) THROW_ERROR("Grid2dBinGroup::Parameters::load: Input stream failed.");
}

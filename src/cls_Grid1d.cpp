// cls_Grid1d.cpp
#include "cls_Grid1d.hpp"
#include "ns_myapp.hpp"
#include "ns_mymacro.hpp"
#include "ns_mylogger.hpp"
#include "spdlog_pch.hpp"
//####################################################
//####################################################
// class Grid1d
//####################################################
//####################################################

// constructor with values 1
Grid1d::Grid1d( const std::string name_in
  , const int nbin_in, const double min_in
  , const double max_in, const double interval_in
  , const double tolerance_ratio )
  : name(name_in), nbin(nbin_in), min(min_in), max(max_in)
  , interval(interval_in)
{
  // Rounding max and min to cleaner values - removed as effectiveness is questionable (2023-04-11).
  // set_round_min_max();
  if( false==is_check_variables(tolerance_ratio) )
   THROW_ERROR3(
    "Grid1d::Grid1d false==is_check_variables(tolerance_ratio)"
    , name, tolerance_ratio);
}

/// @brief Inequality operator
/// @note Does not compare names.
bool Grid1d::operator!=(const Grid1d& other) const
{
  #ifdef NODEBUG
    if (nbin != other.nbin) return true;
    if (min != other.min) return true;
    if (max != other.max) return true;
    if (interval != other.interval) return true;
  #else      
    if (nbin != other.nbin) { LOG_WARN("Grid1d: nbin differs"); return true; }
    if (min != other.min) { LOG_WARN("Grid1d: min differs"); return true; }
    if (max != other.max) { LOG_WARN("Grid1d: max differs"); return true; }
    if (interval != other.interval) { LOG_WARN("Grid1d: interval differs"); return true; }
  #endif
  return false;
}

//==================================================================
// check functions
//==================================================================

/// @brief return true if index is inside
bool Grid1d::is_index_inside( const int index ) const {
  if( index < 0 ) return false;
  if( index > get_nbin()-1 ) return false;
  return true;
}

/// @brief return true if value is inside
bool Grid1d::is_value_inside( const double value ) const
{
  if( value < min ) return false;
  if( value >= max ) return false;
  return true;
}

/// @brief check index and if it is outside, throw error
void Grid1d::check_inside( const int index, srcloc loc ) const
{
  if(is_index_inside(index)==true) return;
  fprintf(stderr,"index is out of range.\n file=%s, function=%s, line=%d\n index=%d\n"
    ,loc.file_name(), loc.function_name(), loc.line(), index);
  LOG_ERROR("index is out of range. file={}, function={}, line={}, index={}"
    ,loc.file_name(), loc.function_name(), loc.line(), index);
  THROW_ERROR("Grid1d::check_inside: index is out of range. index={}", index);
}

/// @brief check value and if it is outside, throw error
void Grid1d::check_inside( const double value, srcloc loc ) const
{
  if(is_value_inside(value)==true) return;
  fprintf(stderr,"value is out of range.\n file=%s, function=%s, line=%d\n value=%E\n"
    ,loc.file_name(), loc.function_name(), loc.line(), value);
  LOG_ERROR("value is out of range. file={}, function={}, line={}, value={}"
    ,loc.file_name(), loc.function_name(), loc.line(), value);
  THROW_ERROR("Grid1d::check_inside: value is out of range. value={}", value);
}

void Grid1d::assign(const std::string name_in, const int nbin_in
  , const double min_in, const double max_in, const double interval_in
  , const double tolerance_ratio )
{
  name = name_in;
  nbin = nbin_in;
  min = min_in;
  max = max_in;
  interval = interval_in;
  // Rounding max and min to cleaner values - removed as effectiveness is questionable (2023-04-11).
  // set_round_min_max();

  if( false==is_check_variables(tolerance_ratio) ){
   THROW_ERROR3(
    "Grid1d::assign false==is_check_variables(tolerance_ratio)"
    , name, tolerance_ratio);
  }
}

bool Grid1d::is_check_variables(const double tolerance_ratio) const
{
  if( nbin < 1 ) {
    LOG_WARN_ND("nbin < 1");
    return false;
  }

  if( min >= max ) {
    LOG_WARN_ND("min >= max");
    return false;
  }

  if( interval <= 0.0 ){
    LOG_WARN_ND("interval <= 0.0");
    return false;
  }
  const double delta = max-min;
  const double value_fmod = std::remainder(delta,interval);
  const double value_fmod_thres = fabs(interval)*tolerance_ratio;
  if( fabs(value_fmod) > value_fmod_thres ){
    LOG_TRACE_ND("nbin={}, max={}, min={}, interval={}, delta={}, value_fmod={}"
      ,nbin,max,min,interval,delta,value_fmod);
    fprintf(stderr,"Grid1d::is_check_variables : std::remainder(max-min,interval) exceeds tolerance\n");
    LOG_WARN_ND("std::remainder(max-min,interval) exceeds tolerance");
  }

  const double interval_tmp = delta/((double)nbin);
  if( fabs(interval-interval_tmp) > fabs(interval)*tolerance_ratio ){
    LOG_WARN_ND("FALSE!! {} , {:E} != {:E}",name,interval,interval_tmp);
    return false;
  }
  return true;
}

double Grid1d::get_lower_value( const int index ) const
{
  check_inside(index);
  const double lower = min + ((double)index + 0.0) *interval;
  return lower;
}
double Grid1d::get_lower_value( const double value ) const
{
  const int index = get_index(value);
  check_inside(index);
  return get_lower_value(index);
}

double Grid1d::get_upper_value( const int index ) const
{
  check_inside(index);
  const double upper = min + ((double)index + 1.0)*interval;
  return upper;
}

double Grid1d::get_upper_value( const double value ) const
{
  const int index = get_index(value);
  check_inside(index);
  return get_upper_value(index);
}

Grid1d::RangeIndices Grid1d::calc_range_indices(
  const double lower, const double upper, const double eps) const
{
  if (upper < lower) THROW_ERROR("Grid1d::calc_range_indices : upper < lower");
  if (eps < 0.0) THROW_ERROR("Grid1d::calc_range_indices : eps < 0.0");
  if (get_nbin() <= 0) THROW_ERROR("Grid1d::calc_range_indices : nbin <= 0");

  const double interval_current = get_interval();
  if( interval_current <= 0.0 )
    THROW_ERROR("Grid1d::calc_range_indices : interval_current <= 0.0");
  const double tol = interval_current * eps;
  // Apply tolerance to shrink the range slightly for conservative bin selection
  const double lower_bound = lower + tol;
  const double upper_bound = upper - tol;

  RangeIndices range;
  int index_begin = get_index(lower_bound);
  if (index_begin == OUT_OF_RANGE_LOWER) {
    index_begin = -1;
  } else if (index_begin == OUT_OF_RANGE_UPPER) {
    THROW_ERROR("Grid1d::calc_range_indices : lower bound is above the upper limit");
  }

  int index_end = get_index(upper_bound);
  if (index_end == OUT_OF_RANGE_UPPER) {
    index_end = get_nbin();
  } else if (index_end == OUT_OF_RANGE_LOWER) {
    THROW_ERROR("Grid1d::calc_range_indices : upper bound is below the lower limit");
  }

  range.start = index_begin;
  range.end   = index_end;
  return range;
}

double Grid1d::get_center_value( const int index ) const
{
  check_inside(index);
  const double center = min + ((double)index + 0.5)*interval;
  return center;
}

double Grid1d::get_center_value( const double value ) const
{
  const int index = get_index(value);
  check_inside(index);
  return get_center_value(index);
}

double Grid1d::get_nearest_value( const double value_in, const bool tf_cnt ) const
{
  double value2 = value_in;
  if( tf_cnt == true ) value2 -= interval*0.5;

  const int index = get_index(value2);
  check_inside(index);

  const double lower = get_lower_value(index);
  const double upper = get_upper_value(index);

  // Return the boundary value (lower or upper) that is closer to value2
  if (std::abs(lower - value2) < std::abs(upper - value2)) {
    return lower;
  } else {
    return upper;
  }
}

int Grid1d::get_index( const double value ) const
{
  // define small value
  const double safe_factor = 1.0;
  const double eps = std::numeric_limits<double>::epsilon() * safe_factor;
  const double small_len = interval * eps;
  // if value is out of range, return OUT_OF_RANGE_LOWER or OUT_OF_RANGE_UPPER
  if( value < min - small_len ){
    LOG_TRACE_ND("name={}, value = {} < min-small_len = {}",name,value,min-small_len);
    fprintf(stderr,"Grid1d::get_index | name=%s, value = %E < min = %E\r",name.c_str(),value,min-small_len);  
    return OUT_OF_RANGE_LOWER;
  }
  if( value >= max + small_len ){
    LOG_TRACE_ND("name={}, value = {} >= max+small_len = {}",name,value,max+small_len);
    fprintf(stderr,"Grid1d::get_index | name=%s, value = %E >= max = %E\r",name.c_str(),value,max+small_len);
    return OUT_OF_RANGE_UPPER;
  }

  // calculate index number from value
  // the meaning of small is to avoid round error
  const int index = (int)floor( ( value - min + small_len )/interval );

  // boundary check: index must be in [0, nbin)
  if( index >= nbin ){
    return OUT_OF_RANGE_UPPER;
  }

  return index;
}

Grid1d Grid1d::get_merged( const double center, int& merge_factor ) const
{
  // if merge_factor is 1, return the copy of this instance
  if( merge_factor == 1 ){
    Grid1d new_axis(*this);
    return new_axis;
  }

  // center should be inside the range
  if( center < get_min() || center >= get_max() ){
    LOG_ERROR_ND("center should be in the raneg of {} <= center <= {}",get_min(),get_max());
    THROW_ERROR_NAME2("center should be in the raneg of " 
    + std::to_string(min) + " <= center <= " + std::to_string(max), center );
  }

  // merge_factor should be >= 2
  if( merge_factor < 1 ){
    LOG_ERROR_ND("merge_factor should be >= 1");
    THROW_ERROR_NAME2("merge_factor should be >= 1", merge_factor);
  }

  // Identify the index that corresponds to the center value
  int icnt = get_index(center);
  if (icnt == OUT_OF_RANGE_LOWER || icnt == OUT_OF_RANGE_UPPER){
    LOG_ERROR_ND("center is out of range");
    THROW_ERROR_NAME2("center is out of range", center);
  }

  // nbin_lower is the number of bins that can be merged on the lower side of the center
  // nbin_upper is the number of bins that can be merged on the upper side of the center
  const int nbin_lower = (icnt) / merge_factor;
  const int nbin_upper = (nbin - 1 - icnt) / merge_factor;

  // If neither side can be merged, adjust merge_factor to allow at least one side to merge.
  if( nbin_lower==0 && nbin_upper==0 ){
    // If upper direction is longer, set merge_factor = nbin - 1 - icnt
    if( nbin - 1 - icnt > icnt ){
      merge_factor = nbin - 1 - icnt;
      LOG_WARN_ND("instance name = {} : nbin_lower==0 && nbin_upper==0.\n nbin - 1 - icnt > icnt.\n merge_factor is changed to nbin - 1 - icnt = {}"
        ,get_name(),merge_factor);
    }
    else{
      // If lower direction is longer, set merge_factor = icnt
      merge_factor = icnt;
      LOG_WARN_ND("instance name = {} : nbin_lower==0 && nbin_upper==0.\n nbin - 1 - icnt <= icnt.\n merge_factor is changed to icnt = {}"
        ,get_name(),merge_factor);
    }
  }

  // Calculate new values
  const std::string name_new = name + "_merged_" + std::to_string(merge_factor);
  const int nbin_new = nbin_lower + nbin_upper;
  const int imin = icnt - nbin_lower * merge_factor;
  const int imax = icnt + nbin_upper * merge_factor -1;
  const double min_new = get_lower_value(imin);
  const double max_new = get_upper_value(imax);
  const double interval_new = get_interval() * static_cast<double>(merge_factor);

  Grid1d new_axis(name_new, nbin_new, min_new, max_new, interval_new);
  return new_axis;
}

Grid1d Grid1d::cut( const double lower, const double upper, const double eps ) const
{
  const RangeIndices range = calc_range_indices(lower, upper, eps);
  const double interval_curr = get_interval();
  const double min_current = get_min();
  const double max_current = get_max();

  double min_new = min_current;
  double max_new = max_current;

  if (range.start >= 0) {
    min_new = get_lower_value(range.start);
  }
  if (range.end <= get_nbin() - 1) {
    max_new = get_upper_value(range.end);
  }

  const int nbin_new = range.end - range.start + 1;
  const std::string name_new = get_name() + "_cut";

  fprintf(stderr,"Grid1d::cut : index_begin=%d, index_end=%d, nbin_new=%d, min_new=%E, max_new=%E\r"
    ,range.start,range.end,nbin_new,min_new,max_new);
  Grid1d new_axis(name_new, nbin_new, min_new, max_new, interval_curr);
  return new_axis;
}


std::array<int,2> Grid1d::get_original_index_min_max(
    const Grid1d &g1_org, const Grid1d &g1_merged, const int index_merged ) const
{
  if(index_merged < 0) THROW_ERROR("index_merged < 0");

  const double lower_value_merged = g1_merged.get_lower_value(index_merged);
  const double upper_value_merged = g1_merged.get_upper_value(index_merged);
  const int index_min = g1_org.get_index(lower_value_merged);
  const int index_max = g1_org.get_index(upper_value_merged)-1;

  fprintf(stderr,"Grid1d::get_original_index_min_max index_min=%d, index_max=%d\r",index_min,index_max);

  return std::array<int,2>{index_min, index_max};
}

Grid1d Grid1d::get_split( const int split_factor ) const
{
  if( split_factor < 1 ) THROW_ERROR_NAME2(" split_factor should be >=1 ",split_factor);

  if( split_factor == 1 ){
    Grid1d new_axis(*this);
    new_axis.set_name( name + "_split_" +std::to_string(split_factor) );
    return new_axis;
  }

  const std::string name_new( name + "_split_" +std::to_string(split_factor) );
  const int nbin_new = get_nbin()*split_factor;
  const double min_new = get_min();
  const double max_new = get_max();
  const double interval_new = get_interval()/(double)split_factor;

  Grid1d new_axis(name_new,nbin_new,min_new,max_new,interval_new);
  return new_axis;
}

void Grid1d::out_info( FILE *fout ) const
{
  fprintf(fout," bins_info name=%s: %d %E %E %E\n"
  ,name.c_str(),nbin,min,max,interval);
}

void Grid1d::out_info(spdlog::level::level_enum log_level) const
{
  LOG_INFO("bins_info name={} : {} {:E} {:E} {:E}", name, nbin, min, max, interval );
}

void Grid1d::save( std::ofstream& ofs ) const
{
  // Pass-through: canonical min/max are written as-is (no shift transform).
  io_binary::write_string(ofs, name);
  io_binary::write_binary(ofs, nbin);
  io_binary::write_binary(ofs, min);
  io_binary::write_binary(ofs, max);
  io_binary::write_binary(ofs, interval);

  if( ofs.fail() ) THROW_ERROR("Grid1d::save: std::ofstream& ofs failed.");
}

void Grid1d::load( std::ifstream& ifs, double tolerance_ratio )
{
  // Pass-through: values are read exactly as written by save().
  std::string name_in = io_binary::read_string(ifs);
  const int nbin_in = io_binary::read_binary<int>(ifs);
  const double min_in = io_binary::read_binary<double>(ifs);
  const double max_in = io_binary::read_binary<double>(ifs);
  const double interval_in = io_binary::read_binary<double>(ifs);

  if( ifs.fail() ) THROW_ERROR("Grid1d::load: std::ifstream& ifs failed.");

  assign(name_in, nbin_in, min_in, max_in, interval_in, tolerance_ratio);
}

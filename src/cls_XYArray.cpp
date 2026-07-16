/// @file cls_XYArray.cpp
/// @brief Implementation of XYArray class for (x,y) coordinate pair collection
#include "cls_XYArray.hpp"
#include "ns_myapp.hpp"
#include "ns_mymacro.hpp"
#include "ns_mylogger.hpp"
#include "spdlog_pch.hpp"

//#########################################################################
//#########################################################################
// class XYArray
//#########################################################################
//#########################################################################

// constructor from file
XYArray::XYArray( const std::filesystem::path &path_in )
{
  set_name(path_in.filename().string());
  LOG_INFO("reading {} ...",path_in.string());

  vec_xy = myapp::read_double2(path_in);
  LOG_INFO("reading done.");

  LOG_INFO("sorting ...");
  sort_x_increasing();
  set_maxmin_xy();
  LOG_INFO("sorting done.");
}

// constructor from vector
XYArray::XYArray( const std::vector<std::pair<double,double>> &vec_double2 )
{
  set_name("made from vec_double2");
  const size_t np = vec_double2.size();
  vec_xy.resize(np);
  for(size_t i=0; i<np; i++) vec_xy.at(i) = vec_double2.at(i);
  set_maxmin_xy();
}

// Recompute min/max values from current data
void XYArray::set_maxmin_xy()
{
  constexpr double big = std::numeric_limits<double>::max();
  double xmin_tmp =  big;
  double xmax_tmp = -big;
  double ymin_tmp =  big;
  double ymax_tmp = -big;
  const size_t np = get_np();
  for(size_t i=0; i<np; i++){
    const auto&[x,y] = this->get_xy(static_cast<int>(i));
    xmin_tmp = std::min(xmin_tmp,x);
    xmax_tmp = std::max(xmax_tmp,x);
    ymin_tmp = std::min(ymin_tmp,y);
    ymax_tmp = std::max(ymax_tmp,y);
  }
  xmin = xmin_tmp;
  xmax = xmax_tmp;
  ymin = ymin_tmp;
  ymax = ymax_tmp;
}

//======================================================================
// Sorting

void XYArray::sort_x_increasing()
{
  mysort::sort_vec_double2_x_increasing(vec_xy);
}

void XYArray::sort_y_increasing()
{
  mysort::sort_vec_double2_y_increasing(vec_xy);
}


// Interpolate y values for given x
std::vector<double> XYArray::get_interp_y( const double x_in ) const
{
  if( x_in <  xmin )
    THROW_ERROR_NAME("XYArray::get_interp_y: x_in < xmin. x_in={}, xmin={}", x_in, xmin);
  if( x_in >= xmax )
    THROW_ERROR_NAME("XYArray::get_interp_y: x_in >= xmax. x_in={}, xmax={}", x_in, xmax);

  std::vector<double> vec_y_interp;
  const size_t np = this->get_np();
  for(size_t i=0; i<np-1; i++){
    const auto&[x0, y0] = this->get_xy(static_cast<int>(i));
    const auto&[x1, y1] = this->get_xy(static_cast<int>(i+1));
    if( x_in <  x0 ) continue;
    if( x_in >= x1 ) continue;
    // x0 <= x_in < x1
    const double y_interp = myapp::linear_interpolation(x0,y0,x1,y1,x_in);
    vec_y_interp.push_back(y_interp);
  }
  return vec_y_interp;
}

// Interpolate x values for given y
std::vector<double> XYArray::get_interp_x( const double y_in ) const
{
  if( y_in <  ymin )
    THROW_ERROR_NAME("XYArray::get_interp_x: y_in < ymin. y_in={}, ymin={}", y_in, ymin);
  if( y_in >= ymax )
    THROW_ERROR_NAME("XYArray::get_interp_x: y_in >= ymax. y_in={}, ymax={}", y_in, ymax);

  std::vector<double> vec_x_interp;

  // Copy and sort by y for interpolation along y-axis
  XYArray xy_arr_yx(*this);
  xy_arr_yx.sort_y_increasing();

  const size_t np = this->get_np();
  for(size_t i=0; i<np-1; i++){
    const auto&[x0, y0] = xy_arr_yx.get_xy(static_cast<int>(i));
    const auto&[x1, y1] = xy_arr_yx.get_xy(static_cast<int>(i+1));
    if( y_in <  y0 ) continue;
    if( y_in >= y1 ) continue;
    // y0 <= y_in < y1
    const double x_interp = myapp::linear_interpolation(y0,x0,y1,x1,y_in);
    vec_x_interp.push_back(x_interp);
  }
  return vec_x_interp;
}

// Compute numerical derivative dy/dx
XYArray XYArray::make_differential() const
{
  const size_t np = this->get_np();
  XYArray xy_arr_diff;
  xy_arr_diff.vec_xy.resize(np-1);
  constexpr double small = std::numeric_limits<double>::epsilon();
  constexpr double large_value = std::numeric_limits<double>::max();
  for(size_t i=0; i<np-1; i++){
    const auto&[x0, y0] = this->get_xy(static_cast<int>(i));
    const auto&[x1, y1] = this->get_xy(static_cast<int>(i+1));
    const double dx = x1-x0;
    if( std::fabs(dx) < small )
      THROW_ERROR_NAME("XYArray::make_differential: dx too small. dx={}, epsilon={}", dx, small);
    const double dy = y1-y0;
    const double dydx = dy/dx;
    // Avoid direct infinity comparison (unreliable on some platforms); compare to large_value instead
    if (dydx >= large_value || dydx <= -large_value) {
      THROW_ERROR_NAME("XYArray::make_differential: dy/dx magnitude exceeds limit. dydx={}", dydx);
    }
    xy_arr_diff.vec_xy.at(i) = std::make_pair(x0,dydx);
  }
  xy_arr_diff.set_maxmin_xy();
  return xy_arr_diff;
}

// Transform x values: x -> pow(10, x)
void XYArray::convert_x_to_pow10()
{
  const size_t np = get_np();
  for(size_t i=0; i<np; i++){
    const auto&[x,y] = get_xy(static_cast<int>(i));
    const double pow10_x = std::pow(10.0, x);
    set_xy(static_cast<int>(i), pow10_x, y);
  }
  set_maxmin_xy();
}

// Transform y values: y -> pow(10, y)
void XYArray::convert_y_to_pow10()
{
  const size_t np = get_np();
  for(size_t i=0; i<np; i++){
    const auto&[x,y] = get_xy(static_cast<int>(i));
    const double pow10_y = std::pow(10.0, y);
    set_xy(static_cast<int>(i), x, pow10_y);
  }
  set_maxmin_xy();
}

// Transform both x and y: (x,y) -> (pow(10,x), pow(10,y))
void XYArray::convert_xy_to_pow10()
{
  const size_t np = get_np();
  for(size_t i=0; i<np; i++){
    const auto&[x,y] = get_xy(static_cast<int>(i));
    const double pow10_x = std::pow(10.0, x);
    const double pow10_y = std::pow(10.0, y);
    set_xy(static_cast<int>(i), pow10_x, pow10_y);
  }
  set_maxmin_xy();
}

// Create copy with x,y converted to pow10(x),pow10(y)
XYArray XYArray::get_converted_xy_to_pow10() const
{
  XYArray xyarr(*this);
  const size_t np = get_np();
  for(size_t i=0; i<np; i++){
    const auto&[x,y] = get_xy(static_cast<int>(i));
    const double pow10_x = std::pow(10.0, x);
    const double pow10_y = std::pow(10.0, y);
    xyarr.set_xy(static_cast<int>(i), pow10_x, pow10_y);
  }
  xyarr.set_maxmin_xy();
  return xyarr;
}

// Write data to text file
void XYArray::out( const std::filesystem::path& pathout ) const
{
  LOG_INFO("to pathout={}", pathout.string());
  if(pathout.empty())
    THROW_ERROR_NAME("XYArray::out: output path is empty");
  FILE *fout = fopen(pathout.c_str(), "wt");
  if(fout == nullptr)
    THROW_ERROR_NAME("XYArray::out: cannot open file. path={}", pathout.string());
  const size_t np = get_np();
  for(size_t i=0; i<np; i++){
    const auto&[x,y] = get_xy(static_cast<int>(i));
    fprintf(fout, "%E %E\n", x, y);
  }
  if(fclose(fout) == EOF)
    THROW_ERROR_NAME("XYArray::out: fclose failed. path={}", pathout.string());
}

// Write data to default file (name + ".tmp")
void XYArray::out() const
{
  const std::filesystem::path pathout = get_name() + ".tmp";
  out(pathout);
}
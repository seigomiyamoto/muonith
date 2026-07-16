/// @file cls_RectAngularBinGroup.cpp
/// @brief Implementation of RectAngularBinGroup class

#include "cls_RectAngularBinGroup.hpp"
#include "ns_io_binary.hpp"
#include "ns_mymacro.hpp"
#include "ns_mylogger.hpp"

#include <cmath>
#include <fstream>
#include <limits>

//==========================================================================
// class RectAngularBinGroup
//==========================================================================

//--------------------------------------------------------------------------
// Constructor from file
//--------------------------------------------------------------------------
RectAngularBinGroup::RectAngularBinGroup(const std::filesystem::path& path_in)
{
  LOG_INFO("Loading RectAngularBinGroup from {}", path_in.string());
  std::ifstream ifs = io_binary::open_ifstream(path_in);

  double xmin_val = std::numeric_limits<double>::max();
  double xmax_val = std::numeric_limits<double>::lowest();
  double ymin_val = std::numeric_limits<double>::max();
  double ymax_val = std::numeric_limits<double>::lowest();

  // File format: space-separated 4-column text
  // Column order: xlow xup ylow yup
  double xlow, xup, ylow, yup;
  while (ifs >> xlow >> xup >> ylow >> yup) {
    // Note: AABB2d constructor takes (xmin, ymin, xmax, ymax) order
    AABB2d rect(xlow, ylow, xup, yup);
    xmin_val = std::min(xmin_val, xlow);
    xmax_val = std::max(xmax_val, xup);
    ymin_val = std::min(ymin_val, ylow);
    ymax_val = std::max(ymax_val, yup);
    vec_rect_.push_back(rect);
  }

  // Set bounding box of base class
  set_xmin(xmin_val);
  set_xmax(xmax_val);
  set_ymin(ymin_val);
  set_ymax(ymax_val);

  ifs.close();

  LOG_INFO("Loaded RectAngularBinGroup from {} (nbin={})", path_in.string(), vec_rect_.size());
}

//--------------------------------------------------------------------------
// Inequality operator
//--------------------------------------------------------------------------
bool RectAngularBinGroup::operator!=(const RectAngularBinGroup& other) const
{
#ifdef NODEBUG
  if (AABB2d::operator!=(other)) return true;
  if (vec_rect_ != other.vec_rect_) return true;
#else
  if (AABB2d::operator!=(other)) { LOG_WARN("RectAngularBinGroup: AABB2d differs");  return true;  }
  if (vec_rect_ != other.vec_rect_) { LOG_WARN("RectAngularBinGroup: vec_rect_ differs"); return true; }
#endif
  return false;
}

//--------------------------------------------------------------------------
// get_index
//--------------------------------------------------------------------------
int RectAngularBinGroup::get_index(double x, double y) const
{
  const int n = get_nbin();
  for (int i = 0; i < n; ++i) {
    const AABB2d& rect = vec_rect_.at(i);
    if (rect.is_inside(x, y)) {
      return i;
    }
  }
  return no_index;
}

//--------------------------------------------------------------------------
// check_no_overlap
//--------------------------------------------------------------------------
void RectAngularBinGroup::check_no_overlap() const
{
  const int n = static_cast<int>(vec_rect_.size());
  for (int i = 0; i < n; ++i) {
    const AABB2d& rect_i = vec_rect_.at(i);
    for (int j = i + 1; j < n; ++j) {
      const AABB2d& rect_j = vec_rect_.at(j);
      if (rect_i.is_overlap(rect_j)) {
        LOG_ERROR("overlap detected");
        LOG_DEBUG("i={}, j={}, bin_i=({},{})..({},{}) bin_j=({},{})..({},{})",
                  i, j,
                  rect_i.xmin(), rect_i.ymin(), rect_i.xmax(), rect_i.ymax(),
                  rect_j.xmin(), rect_j.ymin(), rect_j.xmax(), rect_j.ymax());
        THROW_ERROR("RectAngularBinGroup::check_no_overlap: overlap detected. i={}, j={}", i, j);
      }
    }
  }
}

//--------------------------------------------------------------------------
// is_tessellated
//--------------------------------------------------------------------------
bool RectAngularBinGroup::is_tessellated() const
{
  // 1) Reject empty bin collection
  if (vec_rect_.empty()) {
    THROW_ERROR("RectAngularBinGroup::is_tessellated: vec_rect_ is empty");
  }

  // 2) Compute bounding box area
  const double w = xmax() - xmin();
  const double h = ymax() - ymin();
  if (w <= 0.0) {
    THROW_ERROR("RectAngularBinGroup::is_tessellated: bounding box width is zero or negative");
  }
  if (h <= 0.0) {
    THROW_ERROR("RectAngularBinGroup::is_tessellated: bounding box height is zero or negative");
  }
  const double bounding_area = w * h;

  // 3) Check all bins are within the bounding box
  for (const auto& rect : vec_rect_) {
    if (rect.xmin() < xmin() || rect.xmax() > xmax() ||
        rect.ymin() < ymin() || rect.ymax() > ymax()) {
      THROW_ERROR("RectAngularBinGroup::is_tessellated: bin is outside bounding box");
    }
  }

  // 4) Check for overlaps
  check_no_overlap();

  // 5) Sum bin areas
  double sum_area = 0.0;
  for (const auto& rect : vec_rect_) {
    const double rw = rect.xmax() - rect.xmin();
    const double rh = rect.ymax() - rect.ymin();
    if (rw <= 0.0) {
      THROW_ERROR("RectAngularBinGroup::is_tessellated: bin width is zero or negative");
    }
    if (rh <= 0.0) {
      THROW_ERROR("RectAngularBinGroup::is_tessellated: bin height is zero or negative");
    }
    sum_area += rw * rh;
  }

  // 6) Compare with bounding box area (with floating-point tolerance)
  constexpr double EPS = 1.0e-6;
  const double tol = EPS * bounding_area;

  if (std::fabs(sum_area - bounding_area) > tol) {
    // Area mismatch implies gaps exist
    return false;
  }

  return true;
}

//--------------------------------------------------------------------------
// check_no_void
//--------------------------------------------------------------------------
void RectAngularBinGroup::check_no_void() const
{
  if (!is_tessellated()) {
    LOG_ERROR("void (gap) detected in tessellation");
    THROW_ERROR("RectAngularBinGroup::check_no_void: void (gap) detected in tessellation");
  }
}

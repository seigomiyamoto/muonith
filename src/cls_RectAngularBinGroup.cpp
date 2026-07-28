/// @file cls_RectAngularBinGroup.cpp
/// @brief Implementation of RectAngularBinGroup class

#include "cls_RectAngularBinGroup.hpp"
#include "ns_io_binary.hpp"
#include "ns_mymacro.hpp"
#include "ns_mylogger.hpp"

#include <cmath>
#include <fstream>
#include <limits>
#include <sstream>
#include <string>

//==========================================================================
// class RectAngularBinGroup
//==========================================================================

namespace {

// Append to vec_out the part of rect_a lying outside rect_b, as up to 4
// axis-aligned pieces (left strip, right strip, middle-bottom, middle-top).
void subtract_rect(const AABB2d& rect_a, const AABB2d& rect_b, std::vector<AABB2d>& vec_out)
{
  if (!rect_a.is_overlap(rect_b)) {
    vec_out.push_back(rect_a);
    return;
  }
  if (rect_b.xmin() > rect_a.xmin()) { // left strip
    vec_out.push_back(AABB2d(rect_a.xmin(), rect_b.xmin(), rect_a.ymin(), rect_a.ymax()));
  }
  if (rect_b.xmax() < rect_a.xmax()) { // right strip
    vec_out.push_back(AABB2d(rect_b.xmax(), rect_a.xmax(), rect_a.ymin(), rect_a.ymax()));
  }
  const double xlo = std::max(rect_a.xmin(), rect_b.xmin());
  const double xhi = std::min(rect_a.xmax(), rect_b.xmax());
  if (rect_b.ymin() > rect_a.ymin()) { // middle-bottom strip
    vec_out.push_back(AABB2d(xlo, xhi, rect_a.ymin(), rect_b.ymin()));
  }
  if (rect_b.ymax() < rect_a.ymax()) { // middle-top strip
    vec_out.push_back(AABB2d(xlo, xhi, rect_b.ymax(), rect_a.ymax()));
  }
}

// Rebuild rectangles into mutually disjoint ones covering the same union (OR).
// Each input rectangle keeps only the part not already covered by earlier
// ones; a rectangle fully covered before contributes nothing.
std::vector<AABB2d> normalize_disjoint(const std::vector<AABB2d>& vec_in)
{
  std::vector<AABB2d> vec_out;
  vec_out.reserve(vec_in.size());
  for (const auto& rect : vec_in) {
    std::vector<AABB2d> vec_pieces;
    vec_pieces.push_back(rect);
    const std::size_t n_kept = vec_out.size(); // pieces appended below stay untouched
    for (std::size_t i = 0; i < n_kept && !vec_pieces.empty(); ++i) {
      std::vector<AABB2d> vec_next;
      for (const auto& piece : vec_pieces) subtract_rect(piece, vec_out.at(i), vec_next);
      vec_pieces.swap(vec_next);
    }
    vec_out.insert(vec_out.end(), vec_pieces.begin(), vec_pieces.end());
  }
  return vec_out;
}

} // namespace

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

  // File format: one rectangle per line, "xlow xup ylow yup", optionally
  // preceded by exclude_keyword. Read line by line so that an unreadable line
  // stops the run instead of silently dropping the rest of the file.
  std::string line;
  int lineno = 0;
  while (std::getline(ifs, line)) {
    lineno++;

    std::istringstream iss(line);
    std::string tok_head;
    if (!(iss >> tok_head)) continue;  // blank or whitespace-only line
    if (tok_head.front() == '#') continue; // comment line

    // "exclude" marks a region left out on purpose; anything else must be a number
    double xlow = 0.0;
    const bool tf_exclude = (tok_head == exclude_keyword);
    if (tf_exclude) {
      if (!(iss >> xlow)) {
        THROW_ERROR("RectAngularBinGroup::RectAngularBinGroup: line {} of {} declares '{}' but has no numbers"
          , lineno, path_in.string(), exclude_keyword);
      }
    }
    else {
      std::istringstream iss_head(tok_head);
      if (!(iss_head >> xlow) || !iss_head.eof()) {
        THROW_ERROR("RectAngularBinGroup::RectAngularBinGroup: unreadable token at line {} of {}. token='{}'"
          , lineno, path_in.string(), tok_head);
      }
    }

    double xup = 0.0, ylow = 0.0, yup = 0.0;
    if (!(iss >> xup >> ylow >> yup)) {
      THROW_ERROR("RectAngularBinGroup::RectAngularBinGroup: line {} of {} does not have 4 numbers"
        , lineno, path_in.string());
    }
    std::string tok_extra;
    if (iss >> tok_extra) {
      THROW_ERROR("RectAngularBinGroup::RectAngularBinGroup: extra token at line {} of {}. token='{}'"
        , lineno, path_in.string(), tok_extra);
    }

    // Note: AABB2d constructor takes (xmin, xmax, ymin, ymax) order
    AABB2d rect(xlow, xup, ylow, yup);
    xmin_val = std::min(xmin_val, xlow);
    xmax_val = std::max(xmax_val, xup);
    ymin_val = std::min(ymin_val, ylow);
    ymax_val = std::max(ymax_val, yup);
    if (tf_exclude) vec_exclude_.push_back(rect);
    else            vec_rect_.push_back(rect);
  }

  // Set bounding box of base class
  set_xmin(xmin_val);
  set_xmax(xmax_val);
  set_ymin(ymin_val);
  set_ymax(ymax_val);

  ifs.close();

  // Overlapping excluded regions mean "exclude both" (OR, union). Rebuild
  // them into disjoint rectangles with the same union so that downstream
  // checks keep their no-overlap premise; untouched when nothing overlaps.
  bool tf_exclude_overlap = false;
  for (std::size_t i = 0; i + 1 < vec_exclude_.size() && !tf_exclude_overlap; ++i) {
    for (std::size_t j = i + 1; j < vec_exclude_.size(); ++j) {
      if (vec_exclude_.at(i).is_overlap(vec_exclude_.at(j))) { tf_exclude_overlap = true; break; }
    }
  }
  if (tf_exclude_overlap) {
    const std::size_t n_before = vec_exclude_.size();
    vec_exclude_ = normalize_disjoint(vec_exclude_);
    LOG_INFO("RectAngularBinGroup: overlapping excluded regions treated as union (OR): rebuilt {} region(s) into {} disjoint rectangle(s) in {}"
      , n_before, vec_exclude_.size(), path_in.string());
  }

  // A group rectangle overlapping a declared excluded region is dropped whole
  // instead of stopping the run; get_index() reports its bins as excluded.
  // Group indices stay dense 0..nbin-1 because they are the vector positions.
  if (!vec_exclude_.empty() && !vec_rect_.empty()) {
    std::vector<AABB2d> vec_kept;
    vec_kept.reserve(vec_rect_.size());
    for (const auto& rect : vec_rect_) {
      bool tf_drop = false;
      for (const auto& rect_ex : vec_exclude_) {
        if (rect.is_overlap(rect_ex)) { tf_drop = true; break; }
      }
      if (tf_drop) vec_dropped_.push_back(rect);
      else         vec_kept.push_back(rect);
    }
    vec_rect_.swap(vec_kept);
    if (!vec_dropped_.empty()) {
      LOG_INFO("RectAngularBinGroup: dropped {} group rectangle(s) overlapping {} excluded region(s) in {}. The removed area can be wider than the declared exclusion."
        , vec_dropped_.size(), vec_exclude_.size(), path_in.string());
    }
  }

  LOG_INFO("Loaded RectAngularBinGroup from {} (nbin={}, nbin_exclude={}, nbin_dropped={})"
    , path_in.string(), vec_rect_.size(), vec_exclude_.size(), vec_dropped_.size());
}

//--------------------------------------------------------------------------
// Inequality operator
//--------------------------------------------------------------------------
bool RectAngularBinGroup::operator!=(const RectAngularBinGroup& other) const
{
#ifdef NODEBUG
  if (AABB2d::operator!=(other)) return true;
  if (vec_rect_ != other.vec_rect_) return true;
  if (vec_exclude_ != other.vec_exclude_) return true;
  if (vec_dropped_ != other.vec_dropped_) return true;
#else
  if (AABB2d::operator!=(other)) { LOG_WARN("RectAngularBinGroup: AABB2d differs");  return true;  }
  if (vec_rect_ != other.vec_rect_) { LOG_WARN("RectAngularBinGroup: vec_rect_ differs"); return true; }
  if (vec_exclude_ != other.vec_exclude_) { LOG_WARN("RectAngularBinGroup: vec_exclude_ differs"); return true; }
  if (vec_dropped_ != other.vec_dropped_) { LOG_WARN("RectAngularBinGroup: vec_dropped_ differs"); return true; }
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
  // Not in any group bin: tell a declared exclusion apart from a plain gap
  const int n_exclude = get_nbin_exclude();
  for (int i = 0; i < n_exclude; ++i) {
    const AABB2d& rect = vec_exclude_.at(i);
    if (rect.is_inside(x, y)) {
      return excluded_index;
    }
  }
  // A dropped rectangle counts as excluded as well
  const int n_dropped = get_nbin_dropped();
  for (int i = 0; i < n_dropped; ++i) {
    const AABB2d& rect = vec_dropped_.at(i);
    if (rect.is_inside(x, y)) {
      return excluded_index;
    }
  }
  return no_index;
}

//--------------------------------------------------------------------------
// check_no_overlap
//--------------------------------------------------------------------------
void RectAngularBinGroup::check_no_overlap() const
{
  // Group bins, excluded regions, and dropped rectangles are walked as one
  // list. The only overlap tolerated is dropped x exclude: that overlap is
  // the reason the rectangle was dropped. Every other pair is rejected.
  const int n_group = get_nbin();
  const int n_exclude = get_nbin_exclude();
  std::vector<const AABB2d*> vec_all;
  vec_all.reserve(vec_rect_.size() + vec_exclude_.size() + vec_dropped_.size());
  for (const auto& rect : vec_rect_)    vec_all.push_back(&rect);
  for (const auto& rect : vec_exclude_) vec_all.push_back(&rect);
  for (const auto& rect : vec_dropped_) vec_all.push_back(&rect);

  // kind: 0 = group bin, 1 = exclude, 2 = dropped (by segment of vec_all)
  static constexpr const char* kind_names[3] = {"bin", "exclude", "dropped"};
  const int n = static_cast<int>(vec_all.size());
  for (int i = 0; i < n; ++i) {
    const AABB2d& rect_i = *vec_all.at(i);
    const int kind_i = (i < n_group) ? 0 : (i < n_group + n_exclude) ? 1 : 2;
    for (int j = i + 1; j < n; ++j) {
      const int kind_j = (j < n_group) ? 0 : (j < n_group + n_exclude) ? 1 : 2;
      if (kind_i + kind_j == 3) continue; // dropped x exclude: tolerated
      const AABB2d& rect_j = *vec_all.at(j);
      if (rect_i.is_overlap(rect_j)) {
        const int base_i = (kind_i == 0) ? 0 : (kind_i == 1) ? n_group : n_group + n_exclude;
        const int base_j = (kind_j == 0) ? 0 : (kind_j == 1) ? n_group : n_group + n_exclude;
        LOG_ERROR("overlap detected");
        LOG_DEBUG("i={}, j={}, bin_i=({},{})..({},{}) bin_j=({},{})..({},{})",
                  i, j,
                  rect_i.xmin(), rect_i.ymin(), rect_i.xmax(), rect_i.ymax(),
                  rect_j.xmin(), rect_j.ymin(), rect_j.xmax(), rect_j.ymax());
        THROW_ERROR("RectAngularBinGroup::check_no_overlap: overlap detected. i={} ({}), j={} ({})"
          , i - base_i, kind_names[kind_i], j - base_j, kind_names[kind_j]);
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
  for (const auto& rect : vec_exclude_) {
    if (rect.xmin() < xmin() || rect.xmax() > xmax() ||
        rect.ymin() < ymin() || rect.ymax() > ymax()) {
      THROW_ERROR("RectAngularBinGroup::is_tessellated: excluded region is outside bounding box");
    }
  }

  // 4) Check for overlaps
  check_no_overlap();

  // 5) Sum bin areas. Declared exclusions count as tiled, so a file that
  //    leaves a band out on purpose is not reported as having a gap.
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
  for (const auto& rect : vec_exclude_) {
    const double rw = rect.xmax() - rect.xmin();
    const double rh = rect.ymax() - rect.ymin();
    if (rw <= 0.0) {
      THROW_ERROR("RectAngularBinGroup::is_tessellated: excluded region width is zero or negative");
    }
    if (rh <= 0.0) {
      THROW_ERROR("RectAngularBinGroup::is_tessellated: excluded region height is zero or negative");
    }
    sum_area += rw * rh;
  }
  // A dropped rectangle contributes only its area outside the excluded
  // regions; the overlapped part is already counted by the loop above.
  // Excluded regions never overlap each other, so the subtraction is exact.
  for (const auto& rect : vec_dropped_) {
    const double rw = rect.xmax() - rect.xmin();
    const double rh = rect.ymax() - rect.ymin();
    if (rw <= 0.0) {
      THROW_ERROR("RectAngularBinGroup::is_tessellated: dropped rectangle width is zero or negative");
    }
    if (rh <= 0.0) {
      THROW_ERROR("RectAngularBinGroup::is_tessellated: dropped rectangle height is zero or negative");
    }
    double area_add = rw * rh;
    for (const auto& rect_ex : vec_exclude_) {
      const double ow = std::min(rect.xmax(), rect_ex.xmax()) - std::max(rect.xmin(), rect_ex.xmin());
      const double oh = std::min(rect.ymax(), rect_ex.ymax()) - std::max(rect.ymin(), rect_ex.ymin());
      if (ow > 0.0 && oh > 0.0) area_add -= ow * oh; // subtract the overlapped part
    }
    sum_area += area_add;
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

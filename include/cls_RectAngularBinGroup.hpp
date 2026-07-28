/// @file cls_RectAngularBinGroup.hpp
/// @brief Rectangular angular bin group container
/// @details Defines RectAngularBinGroup class for managing collections of
///          rectangular angular bins (AABB2d) in 2D angular space.
///
/// ## Typical Workflow
/// 1. Load bin definitions from a text file via the file constructor
/// 2. Query which bin contains a given (x, y) coordinate via get_index()
/// 3. Optionally validate via check_no_overlap() and check_no_void()
///
/// ## File Format
/// Space-separated 4-column text file:
///   xlow xup ylow yup
/// Each line defines one rectangular bin.
/// A line starting with the keyword "exclude" declares a region that is
/// deliberately left out of the analysis:
///   exclude xlow xup ylow yup
/// A group-bin rectangle overlapping an excluded region is dropped whole at
/// load time (with a log notice) instead of stopping the run, so the removed
/// area can be wider than the declared exclusion.
/// Excluded regions may overlap each other: overlapping declarations mean the
/// union (OR) of the regions, and are rebuilt at load time into disjoint
/// rectangles covering the same union.
/// Blank lines and lines starting with '#' are skipped. Any other leading
/// token stops the run; the file is never truncated silently.
///
/// ## Coordinate System
/// - x: azimuthal angle (or first angular axis)
/// - y: polar angle (or second angular axis)
/// - Units must match the detector panel's angle_unit (default: Tangent, dimensionless).
///
/// ## Thread Safety
/// - Thread-safe for read operations (const methods).
/// - Not thread-safe for write operations without external synchronization.
#pragma once

#include <filesystem>
#include <vector>

#include "cls_AABB.hpp"

//==========================================================================
/// @class RectAngularBinGroup
/// @brief Container class for a collection of rectangular angular bins.
/// @details Manages a vector of AABB2d objects representing rectangular bins
///          in 2D angular space. Inherits from AABB2d to represent the overall
///          bounding box of all contained bins.
///
/// ## Responsibilities
/// - Load rectangular bin definitions from a text file.
/// - Query which bin contains a given 2D point.
/// - Validate bin geometry (no overlap, complete tessellation).
///
/// ## Invariants
/// - The bounding box (inherited AABB2d) encloses all contained bins,
///   both the group bins and the excluded regions.
/// - Each bin is an AABB2d with non-negative width and height.
/// - Group bins and excluded regions never overlap each other after loading:
///   a group bin overlapping an excluded region is moved to the dropped list
///   by the file constructor.
///
/// ## Usage Example
/// @code
/// // Load from file
/// RectAngularBinGroup group("/path/to/bins.txt");
///
/// // Query bin index for a point
/// int idx = group.get_index(1.2, 0.5);
/// if (idx != RectAngularBinGroup::no_index) {
///     const AABB2d& bin = group.get_rect(idx);
/// }
///
/// // Validate geometry
/// group.check_no_overlap();
/// group.check_no_void();
/// @endcode
//==========================================================================
class RectAngularBinGroup : public AABB2d {
private:
  /// @brief Name of this instance (for debugging/logging purposes)
  std::string name_ = "RectAngularBinGroup";

  /// @brief Collection of rectangular bins
  std::vector<AABB2d> vec_rect_ = {};

  /// @brief Collection of regions deliberately excluded from the analysis
  /// @details Filled from lines beginning with exclude_keyword. Bins whose
  ///          center falls in one of these regions are reported as
  ///          excluded_index by get_index() instead of no_index.
  std::vector<AABB2d> vec_exclude_ = {};

  /// @brief Group rectangles dropped for overlapping an excluded region
  /// @details Filled by the file constructor: a group rectangle overlapping
  ///          any entry of vec_exclude_ is moved here instead of stopping the
  ///          run. Points inside are reported as excluded_index by
  ///          get_index(), and the area outside the excluded regions counts
  ///          as tiled in is_tessellated().
  std::vector<AABB2d> vec_dropped_ = {};

public:
  //======================================================================
  /// @name Constructors and Destructor
  ///@{

  /// @brief Default constructor
  /// @details Creates an empty bin group with no bins.
  RectAngularBinGroup() = default;

  /// @brief Copy constructor
  RectAngularBinGroup(const RectAngularBinGroup& org) = default;

  /// @brief Construct from a text file
  /// @param[in] path_in Path to the bin definition file.
  /// @throws std::runtime_error If the file cannot be opened, or if any line
  ///         is neither blank, a '#' comment, four numbers, nor
  ///         exclude_keyword followed by four numbers.
  /// @note File format: space-separated columns (xlow xup ylow yup) per line,
  ///       optionally preceded by exclude_keyword to declare an excluded region.
  RectAngularBinGroup(const std::filesystem::path& path_in);

  /// @brief Destructor
  ~RectAngularBinGroup() = default;

  ///@} ------------------------------------------------------------------

  //======================================================================
  /// @name Operators
  ///@{

  /// @brief Copy assignment operator
  RectAngularBinGroup& operator=(const RectAngularBinGroup& other) = default;

  /// @brief Move assignment operator
  RectAngularBinGroup& operator=(RectAngularBinGroup&& other) noexcept = default;

  /// @brief Inequality operator
  /// @details Compares base AABB2d, vec_rect_, vec_exclude_, and vec_dropped_.
  ///          The name_ member is not compared.
  /// @param[in] other The other RectAngularBinGroup to compare with.
  /// @return true if not equal, false otherwise.
  bool operator!=(const RectAngularBinGroup& other) const;

  /// @brief Equality operator
  /// @param[in] other The other RectAngularBinGroup to compare with.
  /// @return true if equal, false otherwise.
  bool operator==(const RectAngularBinGroup& other) const {
    return !(*this != other);
  }

  ///@} ------------------------------------------------------------------

  //======================================================================
  /// @name Getter Functions
  ///@{

  /// @brief Get the number of rectangular bins
  /// @return Number of bins in the group.
  int get_nbin() const { return static_cast<int>(vec_rect_.size()); }

  /// @brief Get immutable reference to a bin by index
  /// @param[in] index_in Zero-based index of the bin.
  /// @return Const reference to the AABB2d bin.
  /// @throws std::out_of_range If index_in is out of bounds.
  const AABB2d& get_rect(int index_in) const {
    return vec_rect_.at(index_in);
  }

  /// @brief Get the number of excluded regions
  /// @return Number of regions declared with exclude_keyword.
  int get_nbin_exclude() const { return static_cast<int>(vec_exclude_.size()); }

  /// @brief Get immutable reference to an excluded region by index
  /// @param[in] index_in Zero-based index of the excluded region.
  /// @return Const reference to the AABB2d region.
  /// @throws std::out_of_range If index_in is out of bounds.
  const AABB2d& get_rect_exclude(int index_in) const {
    return vec_exclude_.at(index_in);
  }

  /// @brief Get the number of dropped group rectangles
  /// @return Number of group rectangles dropped for overlapping an excluded region.
  int get_nbin_dropped() const { return static_cast<int>(vec_dropped_.size()); }

  /// @brief Sentinel value indicating no bin found
  static constexpr int no_index = -1;

  /// @brief Sentinel value indicating the point lies in a declared excluded region
  /// @details Distinct from no_index so that a deliberate exclusion can be told
  ///          apart from a gap left by a mistake in the bin-list file.
  static constexpr int excluded_index = -2;

  /// @brief Leading token that declares an excluded region in the bin-list file
  /// @note Matched exactly; "EXCLUDE" and other spellings are rejected.
  static constexpr const char* exclude_keyword = "exclude";

  /// @brief Find the bin index containing a given point
  /// @param[in] x The x-coordinate (first angular axis).
  /// @param[in] y The y-coordinate (second angular axis).
  /// @return Index of the bin containing (x, y), excluded_index if the point
  ///         lies in a declared excluded region or in a dropped rectangle,
  ///         or no_index if neither.
  /// @note Complexity: O(n + m + k) where n = get_nbin(), m = get_nbin_exclude(),
  ///       k = get_nbin_dropped().
  int get_index(double x, double y) const;

  /// @brief Get minimum x-coordinate of the bounding box
  double xmin() const { return AABB2d::xmin(); }

  /// @brief Get minimum y-coordinate of the bounding box
  double ymin() const { return AABB2d::ymin(); }

  /// @brief Get maximum x-coordinate of the bounding box
  double xmax() const { return AABB2d::xmax(); }

  /// @brief Get maximum y-coordinate of the bounding box
  double ymax() const { return AABB2d::ymax(); }

  ///@} ------------------------------------------------------------------

  //======================================================================
  /// @name Modifier Functions
  ///@{

  /// @brief Get mutable reference to a bin by index
  /// @param[in] index_in Zero-based index of the bin.
  /// @return Mutable reference to the AABB2d bin.
  /// @throws std::out_of_range If index_in is out of bounds.
  AABB2d& call_rect(int index_in) {
    return vec_rect_.at(index_in);
  }

  /// @brief Set minimum x-coordinate of the bounding box
  /// @param[in] xmin_in New minimum x value.
  void set_xmin(double xmin_in) {
    AABB2d::set_xmin(xmin_in);
  }

  /// @brief Set minimum y-coordinate of the bounding box
  /// @param[in] ymin_in New minimum y value.
  void set_ymin(double ymin_in) {
    AABB2d::set_ymin(ymin_in);
  }

  /// @brief Set maximum x-coordinate of the bounding box
  /// @param[in] xmax_in New maximum x value.
  void set_xmax(double xmax_in) {
    AABB2d::set_xmax(xmax_in);
  }

  /// @brief Set maximum y-coordinate of the bounding box
  /// @param[in] ymax_in New maximum y value.
  void set_ymax(double ymax_in) {
    AABB2d::set_ymax(ymax_in);
  }

  ///@} ------------------------------------------------------------------

  //======================================================================
  /// @name Validation Functions
  ///@{

  /// @brief Check that no bins overlap
  /// @details Group bins, excluded regions, and dropped rectangles are checked
  ///          together. The only overlap tolerated is a dropped rectangle
  ///          against an excluded region: that overlap is the reason the
  ///          rectangle was dropped. Every other pair is rejected.
  /// @throws std::runtime_error If any two rectangles overlap (except the
  ///         tolerated dropped-vs-excluded pair).
  /// @note Complexity: O((n+m+k)^2) where n = get_nbin(), m = get_nbin_exclude(),
  ///       k = get_nbin_dropped().
  void check_no_overlap() const;

  /// @brief Check if bins completely tessellate the bounding box
  /// @details Excluded regions count towards the tiled area, so a file that
  ///          declares the uncovered band with exclude_keyword still passes.
  ///          A dropped rectangle contributes only the part of its area that
  ///          lies outside the excluded regions, so the overlap is not counted
  ///          twice.
  /// @return true if bins, excluded regions, and dropped rectangles perfectly
  ///         tile the bounding box.
  /// @throws std::runtime_error If vec_rect_ is empty or dimensions are invalid.
  /// @note Complexity: O((n+m+k)^2) due to overlap check.
  bool is_tessellated() const;

  /// @brief Verify bins completely tile the bounding box (throws if not)
  /// @throws std::runtime_error If tessellation check fails (gaps detected).
  void check_no_void() const;

  ///@} ------------------------------------------------------------------
};

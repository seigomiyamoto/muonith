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
///
/// ## Coordinate System
/// - x: azimuthal angle (or first angular axis)
/// - y: polar angle (or second angular axis)
/// - Units depend on the input file; typically radians or degrees.
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
/// - The bounding box (inherited AABB2d) encloses all contained bins.
/// - Each bin is an AABB2d with non-negative width and height.
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
  /// @throws std::runtime_error If the file cannot be opened.
  /// @note File format: space-separated columns (xlow xup ylow yup) per line.
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
  /// @details Compares base AABB2d and vec_rect_. The name_ member is not compared.
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

  /// @brief Sentinel value indicating no bin found
  static constexpr int no_index = -1;

  /// @brief Find the bin index containing a given point
  /// @param[in] x The x-coordinate (first angular axis).
  /// @param[in] y The y-coordinate (second angular axis).
  /// @return Index of the bin containing (x, y), or no_index if not found.
  /// @note Complexity: O(n) where n = get_nbin().
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
  /// @throws std::runtime_error If any two bins overlap.
  /// @note Complexity: O(n^2) where n = get_nbin().
  void check_no_overlap() const;

  /// @brief Check if bins completely tessellate the bounding box
  /// @return true if bins perfectly tile the bounding box without gaps.
  /// @throws std::runtime_error If vec_rect_ is empty or dimensions are invalid.
  /// @note Complexity: O(n^2) due to overlap check.
  bool is_tessellated() const;

  /// @brief Verify bins completely tile the bounding box (throws if not)
  /// @throws std::runtime_error If tessellation check fails (gaps detected).
  void check_no_void() const;

  ///@} ------------------------------------------------------------------
};

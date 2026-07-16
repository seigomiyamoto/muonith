// cls_AABB.cpp
#include "cls_AABB.hpp"
#include "ns_mymacro.hpp"
#include "spdlog_pch.hpp"
#include "ns_mylogger.hpp"
//###########################################
// class AABB2d
//###########################################

bool AABB2d::operator!=(const AABB2d& other) const
{
  if (v2_min != other.v2_min) return true;
  if (v2_max != other.v2_max) return true;
  return false;
}

double AABB2d::min( const int i ) const
{
  if(i<0) THROW_ERROR("AABB2d::min: index out of range. i={}", i);
  if(i>1) THROW_ERROR("AABB2d::min: index out of range. i={}", i);
  return v2_min(i);
};

double AABB2d::max( const int i ) const
{
  if(i<0) THROW_ERROR("AABB2d::max: index out of range. i={}", i);
  if(i>1) THROW_ERROR("AABB2d::max: index out of range. i={}", i);
  return v2_max(i);
}

bool AABB2d::is_inside(const Eigen::Vector2d &v2_pos ) const
{
  // Check if the coordinate is inside the rectangle
  if (v2_pos.x() <  this->xmin()) return false;
  if (v2_pos.y() <  this->ymin()) return false;
  if (v2_pos.x() >= this->xmax()) return false;
  if (v2_pos.y() >= this->ymax()) return false;
  return true;
}

bool AABB2d::is_inside(const Eigen::Vector2d &v2_pos, const Angle &ang_theta) const
{
  // Compute the centroid of the rectangle
  Eigen::Vector2d v2_center = this->center();

  // Create inverse rotation matrix
  Eigen::Rotation2Dd inv_rot(-ang_theta.rad());

  // Apply inverse rotation around v2_center
  Eigen::Vector2d v2_rotated = inv_rot * (v2_pos - v2_center) + v2_center;

  // Check if the inversely rotated coordinate is inside the rectangle
  return this->is_inside(v2_rotated);
}

AABB2DVertices AABB2d::get_vertices() const
{
  const Eigen::Vector2d v2_xmym(v2_min.x(),v2_min.y());
  const Eigen::Vector2d v2_xpym(v2_max.x(),v2_min.y());
  const Eigen::Vector2d v2_xmyp(v2_min.x(),v2_max.y());
  const Eigen::Vector2d v2_xpyp(v2_max.x(),v2_max.y());
  return {v2_xmym, v2_xpym
        , v2_xmyp, v2_xpyp};
};


bool AABB2d::is_overlap(const AABB2d &other) const
{
  if ( other.v2_max.x() <= v2_min.x() ) return false;
  if ( v2_max.x() <= other.v2_min.x() ) return false;
  if ( other.v2_max.y() <= v2_min.y() ) return false;
  if ( v2_max.y() <= other.v2_min.y() ) return false;
  return true;
}


AABB2d::Adjacency AABB2d::is_adjacent(const AABB2d &other) const
{

  if (is_overlap(other)) {
    return Adjacency::Overlap; // overlapping case
  }
  if (xmax() == other.xmin() && ymin() < other.ymax() && ymax() > other.ymin()) {
    return Adjacency::Xp; // other rectangle is adjacent on the right
  }
  if (xmin() == other.xmax() && ymin() < other.ymax() && ymax() > other.ymin()) {
    return Adjacency::Xm;  // other rectangle is adjacent on the left
  }
  if (ymax() == other.ymin() && xmin() < other.xmax() && xmax() > other.xmin()) {
    return Adjacency::Yp;   // other rectangle is adjacent above
  }
  if (ymin() == other.ymax() && xmin() < other.xmax() && xmax() > other.xmin()) {
    return Adjacency::Ym;// other rectangle is adjacent below
  }
  return Adjacency::None; // not adjacent
}

void AABB2d::checkMinMax() const
{
  if (v2_min.x() > v2_max.x() || v2_min.y() > v2_max.y()) {
    THROW_ERROR("AABB2d::checkMinMax: min must be smaller than max. min=({},{}), max=({},{})",
                v2_min.x(), v2_min.y(), v2_max.x(), v2_max.y());
  }
};

//###########################################
// class AABB3d
//###########################################
void AABB3d::checkMinMax() const
{
  if ( v3_min.x() > v3_max.x() ){
    THROW_ERROR("AABB3d::checkMinMax: min must be smaller than max. v3_min.x()={}, v3_max.x()={}", v3_min.x(), v3_max.x());
  }
  if ( v3_min.y() > v3_max.y() ){
    THROW_ERROR("AABB3d::checkMinMax: min must be smaller than max. v3_min.y()={}, v3_max.y()={}", v3_min.y(), v3_max.y());
  }
  if ( v3_min.z() > v3_max.z() ){
    THROW_ERROR("AABB3d::checkMinMax: min must be smaller than max. v3_min.z()={}, v3_max.z()={}", v3_min.z(), v3_max.z());
  }
};

AABB3d::Adjacency AABB3d::is_adjacent(const AABB3d &other) const
{
  if (is_overlap(other)) {
    return Adjacency::Overlap; // overlapping case
  }
  if (xmax() == other.xmin() && ymin() < other.ymax() && ymax() > other.ymin()) {
    return Adjacency::Xp; // other cuboid is adjacent on the right
  }
  if (xmin() == other.xmax() && ymin() < other.ymax() && ymax() > other.ymin()) {
    return Adjacency::Xm;  // other cuboid is adjacent on the left
  }
  if (zmax() == other.zmin() && xmin() < other.xmax() && xmax() > other.xmin()) {
    return Adjacency::Zp;   // other cuboid is adjacent above
  }
  if (zmin() == other.zmax() && xmin() < other.xmax() && xmax() > other.xmin()) {
    return Adjacency::Zm;// other cuboid is adjacent below
  }
  if (ymax() == other.ymin() && xmin() < other.xmax() && xmax() > other.xmin()) {
    return Adjacency::Yp;   // other cuboid is adjacent behind
  }
  if (ymin() == other.ymax() && xmin() < other.xmax() && xmax() > other.xmin()) {
    return Adjacency::Ym;// other cuboid is adjacent in front
  }
  return Adjacency::None; // not adjacent
}

void AABB3d::disp(FILE *fout) const
{
  fprintf(fout,"%s | x:%.1lf-%.1lf, y:%.1lf-%.1lf, z:%.1lf-%.1lf\n"
    ,__FUNCTION__
    ,this->xmin(),this->xmax()
    ,this->ymin(),this->ymax()
    ,this->zmin(),this->zmax()
  );
};


bool AABB3d::operator!=(const AABB3d& other) const
{
  if (v3_min != other.v3_min) return true;
  if (v3_max != other.v3_max) return true;
  return false;
}

double AABB3d::min( const int i ) const
{
  if(i<0) THROW_ERROR("AABB3d::min: index out of range. i={}", i);
  if(i>2) THROW_ERROR("AABB3d::min: index out of range. i={}", i);
  return v3_min(i);
}

double AABB3d::max( const int i ) const
{
  if(i<0) THROW_ERROR("AABB3d::max: index out of range. i={}", i);
  if(i>2) THROW_ERROR("AABB3d::max: index out of range. i={}", i);
  return v3_max(i);
}


bool AABB3d::is_inside(const Eigen::Vector3d &v3_pos ) const
{
  // Check if the coordinate is inside the cuboid
  if (v3_pos.x() <  this->xmin()) return false;
  if (v3_pos.y() <  this->ymin()) return false;
  if (v3_pos.z() <  this->zmin()) return false;
  if (v3_pos.x() >= this->xmax()) return false;
  if (v3_pos.y() >= this->ymax()) return false;
  if (v3_pos.z() >= this->zmax()) return false;
  return true;
}

/// @brief check overlap exit or not
bool AABB3d::is_overlap(const AABB3d &other) const
{
  if ( other.v3_max.x() <= v3_min.x() ) return false;
  if ( v3_max.x() <= other.v3_min.x() ) return false;
  if ( other.v3_max.y() <= v3_min.y() ) return false;
  if ( v3_max.y() <= other.v3_min.y() ) return false;
  if ( other.v3_max.z() <= v3_min.z() ) return false;
  if ( v3_max.z() <= other.v3_min.z() ) return false;
  return true;
}

AABB3DVertices AABB3d::get_vertices() const
{
  // z lower
  const Eigen::Vector3d v3_xmymzm(v3_min.x(),v3_min.y(),v3_min.z()); // x-,y-,z-
  const Eigen::Vector3d v3_xpymzm(v3_max.x(),v3_min.y(),v3_min.z()); // x+,y-,z-
  const Eigen::Vector3d v3_xmypzm(v3_min.x(),v3_max.y(),v3_min.z()); // x-,y+,z-
  const Eigen::Vector3d v3_xpypzm(v3_max.x(),v3_max.y(),v3_min.z()); // x+,y+,z-
  // z upper
  const Eigen::Vector3d v3_xmymzp(v3_min.x(),v3_min.y(),v3_max.z()); // x-,y-,z+
  const Eigen::Vector3d v3_xpymzp(v3_max.x(),v3_min.y(),v3_max.z()); // x+,y-,z+
  const Eigen::Vector3d v3_xmypzp(v3_min.x(),v3_max.y(),v3_max.z()); // x-,y+,z+
  const Eigen::Vector3d v3_xpypzp(v3_max.x(),v3_max.y(),v3_max.z()); // x+,y+,z+

  return {v3_xmymzm, v3_xpymzm, v3_xmypzm, v3_xpypzm
        , v3_xmymzp, v3_xpymzp, v3_xmypzp, v3_xpypzp};
};

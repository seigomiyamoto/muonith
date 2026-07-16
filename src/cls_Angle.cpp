// cls_Angle.cpp
#include "cls_Angle.hpp"
#include "spdlog_pch.hpp"
#include "ns_mylogger.hpp"
#include "ns_mymacro.hpp"
#include <cmath>
#include <iostream>

//######################################
// class Angle
//######################################

Angle::Angle(const double value, const Unit unit)
  : value_(value), unit_(unit) {}

double Angle::rad() const
{
  return (unit_ == Unit::Radian) ? value_ : value_ * M_PI / 180.0;
}

double Angle::deg() const
{
  return (unit_ == Unit::Degree) ? value_ : value_ * 180.0 / M_PI;
}


void Angle::out( FILE *fout ) const
{
  if (getUnit() == Angle::Unit::Degree) {
    fprintf(fout, "Angle = %lf %s ", deg(), unitToString(getUnit()).c_str());
    fprintf(fout, "( = %lf %s)\n", rad(), unitToString(Angle::Unit::Radian).c_str());
  }
  else if (getUnit() == Angle::Unit::Radian) {
    fprintf(fout, "Angle = %lf %s ", rad(), unitToString(getUnit()).c_str());
    fprintf(fout, "( = %lf %s)\n", deg(), unitToString(Angle::Unit::Degree).c_str());
  }
  else {
    fprintf(fout, "Current unit = unknown.\n");
  }
}

std::string Angle::unitToString(Unit unit)
{
  switch (unit) {
    case Unit::Radian : return "Radian";
    case Unit::Degree : return "Degree";
    default : return "Unknown";
  }
}

Angle Angle::operator+(const Angle& rhs) const
{
  double newValue = (unit_ == Unit::Radian) ? rad() + rhs.rad() : deg() + rhs.deg();
  return Angle(newValue, unit_);
}

Angle Angle::operator-(const Angle& rhs) const
{
  double newValue = (unit_ == Unit::Radian) ? rad() - rhs.rad() : deg() - rhs.deg();
  return Angle(newValue, unit_);
}

Angle Angle::operator-() const
{
  return Angle(-value_, unit_);
}

Angle Angle::operator*(const double scalar) const
{
  return Angle(value_ * scalar, unit_);
}

Angle Angle::operator/(const double scalar) const
{
  if (scalar == 0) {
    THROW_ERROR("Angle::operator/: Division by zero is not allowed.");
  }
  return Angle(value_ / scalar, unit_);
}

std::ostream& operator<<(std::ostream& os, const Angle& angle)
{
  os << angle.rad() << " rad (" << angle.deg() << " deg)";
  return os;
}

void Angle::setRadian(const double radian)
{
  value_ = radian;
  unit_ = Unit::Radian;
}

void Angle::setDegree(const double degree)
{
  value_ = degree;
  unit_ = Unit::Degree;
}

void Angle::setTangent(double tangent)
{
  if (std::isnan(tangent)) {
    THROW_ERROR("Angle::setTangent: tangent is NaN.");
  }
  const double ang_rad = std::atan(tangent); // principal value in (-π/2, π/2]
  if (unit_ == Unit::Radian) {
    value_ = ang_rad;
  } else {
    value_ = ang_rad * 180.0 / M_PI;
  }
}

void Angle::negate()
{
  value_ = -value_;
}

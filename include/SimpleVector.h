//
// Created by paolo on 04/11/2023.
//

#ifndef DEEMBROJECT_INCLUDE_SIMPLEVECTOR_H
#define DEEMBROJECT_INCLUDE_SIMPLEVECTOR_H

#include <cmath>

template<typename T>
class SimpleVector {
  // coordinate cartesiane
  T x_;
  T y_;
  T z_;

  SimpleVector(const T x, const T y, const T z) : x_{x}, y_{y}, z_{z}{};

 public:

  static  SimpleVector createEmpty() {
    return SimpleVector(0, 0, 0);
  }

  static SimpleVector createCartesian(const T x, const T y, const T z) {
        return SimpleVector(x, y, z);
  }

  static SimpleVector createPolar(double phi, double theta, double r) {
    double rp = r * std::cos(theta);// raggio piano

    T x = rp * std::cos(phi);
    T y = rp * std::sin(phi);
    T z = r * std::sin(theta);

    return SimpleVector(x, y, z);
  }

  SimpleVector operator+(const SimpleVector& other) const {
      return createCartesian(x_ + other.x_, y_ + other.y_, z_ + other.z_);
  }

  SimpleVector operator-(const SimpleVector& other) const {
      return createCartesian(x_ - other.x_, y_ - other.y_, z_ - other.z_);
  }

  auto operator*(const SimpleVector& other) const {
      return x_ * other.x_ + y_ * other.y_ + z_ * other.z_;
  }

  template<typename S>
  SimpleVector operator*(const S scalar) const {
      return createCartesian(x_ * scalar, y_ * scalar, z_ * scalar);
  }

  template<typename S>
  SimpleVector operator/(const S scalar) const {
      return createCartesian(x_ / scalar, y_ / scalar, z_ / scalar);
  }

  T x() const {
      return x_;
  }

  T y() const {
      return y_;
  }

  T z() const {
      return z_;
  }

  [[nodiscard]] double magnitude() const {
    return std::sqrt(x_ * x_ + y_ * y_ + z_ * z_);
  }

  [[nodiscard]] double phi() const {
    double phi;
    if (x_ == 0) {
      if (y_ < 0) {
        phi = -M_PI / 2;
      } else {
        phi = M_PI / 2;
      }
    } else {
      phi = std::atan(y_ / x_);
      if (x_ < 0) {
        phi += M_PI / 2;
      }
    }
    return phi;
  }

  [[nodiscard]] double theta() const {
    double rp = std::sqrt(x_ * x_ + y_ * y_);// raggio piano

    double theta;
    if (rp == 0) {
      if (z_ < 0) {
        theta = -M_PI / 2;
      } else {
        theta = M_PI / 2;
      }
    } else {
      theta = std::atan(z_ / rp);
    }
    return theta;
  }
};

#endif//DEEMBROJECT_INCLUDE_SIMPLEVECTOR_H

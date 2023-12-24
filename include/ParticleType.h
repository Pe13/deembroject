//
// Created by paolo on 02/11/2023.
//

#ifndef DEEMBROJECT_INCLUDE_PARTICLETYPE_H_
#define DEEMBROJECT_INCLUDE_PARTICLETYPE_H_

#include <map>
#include <string>

enum Type {
  undefined = -1,
  pp = 0,
  pm = 1,
  Kp = 2,
  Km = 3,
  Pp = 4,
  Pm = 5,
  Ks = 6
};

class ParticleType {

  const std::string name_;
  const bool isResonance_;
  const double mass_;
  const int charge_;
  const double width_{};

 public:
  ParticleType();
  ParticleType(const std::string &name, double mass, int charge);
  ParticleType(const std::string &name, double mass, int charge, double width);

  ParticleType &operator=(ParticleType const &particle_type);

  static int getNParticleTypes();
  static const std::map<int, std::string> &getParticleNames();

  [[nodiscard]] bool isResonance() const;
  [[nodiscard]] const std::string &getName() const;
  [[nodiscard]] double getMass() const;
  [[nodiscard]] int getCharge() const;
  [[nodiscard]] double getWidth() const;
  void print() const;
};

#endif// DEEMBROJECT_INCLUDE_PARTICLETYPE_H_

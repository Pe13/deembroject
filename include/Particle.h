//
// Created by pol on 02/11/2023.
//

#ifndef DEEMBROJECT_INCLUDE_PARTICLE_H_
#define DEEMBROJECT_INCLUDE_PARTICLE_H_

#include "ParticleType.h"

#include "SimpleVector.h"
#include <array>

class Particle {

  static std::array<ParticleType, 10> particleTypes_;// lista di tutti i tipi ordinati
  static int nParticleTypes_;                        // numero di tipi validi

  Type typeIndex_;

  SimpleVector<double> p_;

  void boost(const SimpleVector<double> &other);

 public:
  Particle(Type type, double px = 0, double py = 0, double pz = 0);
  Particle();

  [[nodiscard]] Type getType() const;
  bool setType(Type type);
  [[nodiscard]] const SimpleVector<double> &getP() const;
  void setP(const SimpleVector<double> &p);
  [[nodiscard]] double getMass() const;
  [[nodiscard]] int getCharge() const;
  [[nodiscard]] double getEnergy() const;
  [[nodiscard]] double invMass(const Particle &other) const;

  static int addParticleType(const std::string &name, double mass, int charge);
  static int addParticleType(const std::string &name, double mass, int charge, double width);

  static int getNParticleTypes();

  int decay2body(Particle &dau1, Particle &dau2) const;

  void print() const;
  static void printParticles();
};

#endif// DEEMBROJECT_INCLUDE_PARTICLE_H_

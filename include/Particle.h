//
// Created by pol on 02/11/2023.
//

#ifndef DEEMBROJECT_INCLUDE_PARTICLE_H_
#define DEEMBROJECT_INCLUDE_PARTICLE_H_

#include "ParticleType.h"

#include <array>

class Particle {

  static std::array<ParticleType, 10> particleTypes_;// lista di tutti i tipi ordinati
  static int nParticleTypes_;                        // numero di tipi validi

  Type typeIndex_;
  // componenti dell'impulso
  double px_;
  double py_;
  double pz_;

  void boost(double bx, double by, double bz);

 public:
  Particle(Type type, double px = 0, double py = 0, double pz = 0);
  Particle();

  [[nodiscard]] Type getType() const;
  [[nodiscard]] bool setType(Type type);
  [[nodiscard]] std::array<double, 3> getP() const;
  void setP(double px, double py, double pz);
  [[nodiscard]] double getMass() const;
  [[nodiscard]] double getEnergy() const;
  [[nodiscard]] double InvMass(const Particle &other) const;

  static int addParticleType(const std::string &name, double mass, int charge);
  static int addParticleType(const std::string &name, double mass, int charge, double width);

  int decay2body(Particle &dau1, Particle &dau2) const;

  void print() const;
  static void printParticles();
};

#endif// DEEMBROJECT_INCLUDE_PARTICLE_H_

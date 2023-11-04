//
// Created by paolo on 02/11/2023.
//

#include "Particle.h"

#include <TLorentzVector.h>

#include <algorithm>
#include <cmath>
#include <cstdlib>
#include <iostream>
#include <stdexcept>

std::array<ParticleType, 10> Particle::particleTypes_ = {
    ParticleType(),                       // non inizializzati
    ParticleType(),
    ParticleType(),
    ParticleType(),
    ParticleType(),
    ParticleType(),
    ParticleType(),
    ParticleType(),
    ParticleType(),
    ParticleType(),
};

int Particle::nParticleTypes_ = 0

void Particle::boost(double bx, double by, double bz) {

  double energy = getEnergy();

  //boost this Lorentz vector
  double b2 = bx * bx + by * by + bz * bz;
  double gamma = 1.0 / sqrt(1.0 - b2);
  double bp = bx * px_ + by * py_ + bz * pz_;
  double gamma2 = b2 > 0 ? (gamma - 1.0) / b2 : 0.0;

  px_ += gamma2 * bp * bx + gamma * bx * energy;
  py_ += gamma2 * bp * by + gamma * by * energy;
  pz_ += gamma2 * bp * bz + gamma * bz * energy;
}

Particle::Particle(Type type, double px, double py, double pz) : typeIndex_{type}, px_{px}, py_{py}, pz_{pz} {
  if (typeIndex_ >= nParticleTypes_ || type < 0) {
    throw std::runtime_error("Il tipo di particella richiesta non esiste");
  }
}

Particle::Particle(): typeIndex_{undefined}, px_{0}, py_{0}, pz_{0} {}

Type Particle::getType() const { return typeIndex_; }

bool Particle::setType(Type type) {
  if (type >= nParticleTypes_ || type < 0) {
    return false;
  }
  typeIndex_ = type;
  return true;
}

std::array<double, 3> Particle::getP() const { return {px_, py_, pz_}; }

void Particle::setP(const double px, const double py, const double pz) {
  px_ = px;
  py_ = py;
  pz_ = pz;
}

double Particle::getMass() const { return particleTypes_[typeIndex_].getMass(); }

double Particle::getEnergy() const {
  return std::sqrt(getMass() * getMass() + px_ * px_ + py_ * py_ + pz_ * pz_);
}

double Particle::InvMass(const Particle &other) const {
  auto thisImpulse = TLorentzVector(px_, py_, pz_, 0);
  auto otherImpulse = TLorentzVector(other.px_, other.py_, other.pz_, 0);

  return std::sqrt((getEnergy() + other.getEnergy()) * (getEnergy() + other.getEnergy())
                   + (thisImpulse + otherImpulse) * (thisImpulse + otherImpulse));
}

int Particle::addParticleType(const std::string &name, double mass, int charge) {
  if (nParticleTypes_ >= particleTypes_.size()) {
    return 0;
  }
//  particleTypes_[nParticleTypes_] = ParticleType(name, mass, charge);
  particleTypes_[nParticleTypes_] = ParticleType(name, mass, charge);
  nParticleTypes_++;
  return nParticleTypes_;
}
int Particle::addParticleType(const std::string &name, double mass, int charge, double width) {
  if (nParticleTypes_ >= particleTypes_.size()) {
    return 0;
  }
  particleTypes_[nParticleTypes_] = ParticleType(name, mass, charge, width);
  nParticleTypes_++;
  return nParticleTypes_;
}

int Particle::decay2body(Particle &dau1, Particle &dau2) const {
  if (getMass() == 0.0) {
    printf("Decayment cannot be preformed if mass is zero\n");
    return 1;
  }

  double massMot = getMass();
  double massDau1 = dau1.getMass();
  double massDau2 = dau2.getMass();

  if (particleTypes_[typeIndex_].isResonance()) { // add width effect

    // gaussian random numbers

    float x1, x2, w, y1;

    double invnum = 1. / RAND_MAX;
    do {
      x1 = 2.0 * rand() * invnum - 1.0;
      x2 = 2.0 * rand() * invnum - 1.0;
      w = x1 * x1 + x2 * x2;
    } while (w >= 1.0);

    w = std::sqrt((-2.0 * log(w)) / w);
    y1 = x1 * w;

    massMot += particleTypes_[typeIndex_].getWidth() * y1;
  }

  if (massMot < massDau1 + massDau2) {
    printf("Decayment cannot be preformed because mass is too low in this channel\n");
    return 2;
  }

  double pout = sqrt((massMot * massMot - (massDau1 + massDau2) * (massDau1 + massDau2)) * (massMot * massMot - (massDau1 - massDau2) * (massDau1 - massDau2))) / massMot * 0.5;

  double norm = 2 * M_PI / RAND_MAX;

  double phi = rand() * norm;
  double theta = rand() * norm * 0.5 - M_PI / 2.;
  dau1.setP(pout * sin(theta) * cos(phi), pout * sin(theta) * sin(phi), pout * cos(theta));
  dau2.setP(-pout * sin(theta) * cos(phi), -pout * sin(theta) * sin(phi), -pout * cos(theta));

  double energy = sqrt(px_ * px_ + py_ * py_ + pz_ * pz_ + massMot * massMot);

  double bx = px_ / energy;
  double by = py_ / energy;
  double bz = pz_ / energy;

  dau1.boost(bx, by, bz);
  dau2.boost(bx, by, bz);

  return 0;
}

void Particle::print() const {
  std::cout << "Indice del tipo di particella: " << typeIndex_ << '\n'
            << "Nome della particella: " << particleTypes_[typeIndex_].getName() << '\n'
            << "Componenti dell'impulso:\n"
            << "\tx: " << px_ << '\n'
            << "\ty: " << py_ << '\n'
            << "\tz: " << pz_ << '\n';
}

void Particle::printParticles() {
  for (int i = 0; i < nParticleTypes_; i++) {
    particleTypes_[i].print();
    std::cout << '\n';
  }
}


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
    ParticleType(),// non inizializzati
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

int Particle::nParticleTypes_ = 0;

void Particle::boost(const SimpleVector<double> &other) {

  double energy = getEnergy();

  //boost this Lorentz vector
  double b2 = other * other;
  double gamma = 1.0 / sqrt(1.0 - b2);
  double bp = p_ * other;
  double gamma2 = b2 > 0 ? (gamma - 1.0) / b2 : 0.0;

  p_ = p_ + other * (gamma2 * bp + gamma * energy);

  //  px_ += gamma2 * bp * bx + gamma * bx * energy;
  //  py_ += gamma2 * bp * by + gamma * by * energy;
  //  pz_ += gamma2 * bp * bz + gamma * bz * energy;
}

Particle::Particle(Type type, double px, double py, double pz) : typeIndex_{type}, p_{SimpleVector<double>::createCartesian(px, py, pz)} {
  if (typeIndex_ >= nParticleTypes_ || type < 0) {
    throw std::runtime_error("Il tipo di particella richiesta non esiste");
  }
}

Particle::Particle() : typeIndex_{undefined}, p_{SimpleVector<double>::createEmpty()} {}

Type Particle::getType() const { return typeIndex_; }

bool Particle::setType(Type type) {
  if (type >= nParticleTypes_) {
    return false;
  }
  typeIndex_ = type;
  return true;
}

const SimpleVector<double> &Particle::getP() const { return p_; }

void Particle::setP(const SimpleVector<double> &p) {
  p_ = p;
}

double Particle::getMass() const { return particleTypes_[typeIndex_].getMass(); }

double Particle::getEnergy() const {
  return std::sqrt(getMass() * getMass() + p_ * p_);
}

double Particle::InvMass(const Particle &other) const {
  SimpleVector totalImpulse = p_ + other.p_;
  double totalEnergy = getEnergy() + other.getEnergy();

  return std::sqrt(totalEnergy * totalEnergy + totalImpulse * totalImpulse);
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

  if (particleTypes_[typeIndex_].isResonance()) {// add width effect

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
  double theta = rand() * norm * 0.5 - M_PI / 2.;                                             // [- pi/2, pi/2]
//  dau1.setP(pout * sin(theta) * cos(phi), pout * sin(theta) * sin(phi), pout * cos(theta));   // z > 0
//  dau2.setP(-pout * sin(theta) * cos(phi), -pout * sin(theta) * sin(phi), -pout * cos(theta));// z < 0

  dau1.setP(SimpleVector<double>::createPolar(phi, theta, pout));
  dau2.setP(SimpleVector<double>::createPolar(phi, theta, pout) * -1);

  double energy = sqrt(p_ * p_ + massMot * massMot);

//  double bx = px_ / energy;
//  double by = py_ / energy;
//  double bz = pz_ / energy;

  SimpleVector<double> b = p_ / energy;

  dau1.boost(b);
  dau2.boost(b);

  return 0;
}

void Particle::print() const {
  std::cout << "Indice del tipo di particella: " << typeIndex_ << '\n'
            << "Nome della particella: " << particleTypes_[typeIndex_].getName() << '\n'
            << "Componenti dell'impulso:\n"
            << "\tx: " << p_.x() << '\n'
            << "\ty: " << p_.y() << '\n'
            << "\tz: " << p_.z() << '\n';
}

void Particle::printParticles() {
  for (int i = 0; i < nParticleTypes_; i++) {
    particleTypes_[i].print();
    std::cout << '\n';
  }
}

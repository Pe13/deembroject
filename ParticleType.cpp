//
// Created by paolo on 02/11/2023.
//

#include "ParticleType.h"

#include <iostream>
#include <stdexcept>

ParticleType::ParticleType()
    : name_{"undefined"}, isResonance_{false}, mass_{-1}, charge_{1},
      width_{0} {}

ParticleType::ParticleType(const std::string &name, double mass, int charge)
    : name_{name}, isResonance_{false}, mass_{mass}, charge_{charge} {
  if (mass <= 0) {
    throw std::runtime_error("La massa non può essere negativa!");
  }
}

ParticleType::ParticleType(const std::string &name, double mass, int charge,
                           double width)
    : name_{name}, isResonance_{true}, mass_{mass}, charge_{charge},
      width_{width} {
  if (mass <= 0) {
    throw std::runtime_error("La massa non può essere negativa!");
  }
}

ParticleType &ParticleType::operator=(ParticleType const &particle_type) {
  auto name = const_cast<std::string *>(&name_);
  *name = particle_type.name_;

  auto mass = const_cast<double *>(&mass_);
  *mass = particle_type.mass_;

  auto isResonance = const_cast<bool *>(&isResonance_);
  *isResonance = particle_type.isResonance_;

  auto charge = const_cast<int *>(&charge_);
  *charge = particle_type.charge_;

  auto width = const_cast<double *>(&width_);
  *width = particle_type.width_;

  return *this;
}

bool ParticleType::isResonance() const { return isResonance_; }

const std::string &ParticleType::getName() const { return name_; }
// la & serve per passare solo il punto e non la stringa intera, che di per sè
// è molto pesante

double ParticleType::getMass() const { return mass_; }

int ParticleType::getCharge() const { return charge_; }

double ParticleType::getWidth() const { return width_; }

void ParticleType::print() const {
  std::cout << "Tipo: " << name_ << '\n'
            << "Massa: " << mass_ << '\n'
            << "Carica: " << charge_ << '\n';

  if (isResonance_) {
    std::cout << "Larghezza" << width_ << '\n';
  }
}

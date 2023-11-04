//
// Created by paolo on 02/11/2023.
//

#include "Particle.h"

#include <array>
#include <algorithm>
#include <execution>
#include <TRandom.h>

int main() {
  Particle::addParticleType("\\pi +", 0.139657, 1);  // pione positivo
  Particle::addParticleType("\\pi -", 0.139657, 1);  // pione negativo
  Particle::addParticleType("K+", 0.49367, 1);       // kaone positivo
  Particle::addParticleType("K-", 0.49367, -1);      // kaone negativo
  Particle::addParticleType("P+", 0.93827, 1);       // protone
  Particle::addParticleType("P-", 0.93827, -1);      // anti protone
  Particle::addParticleType("K*", 0.89166, 0, 0.050);// kaone *

  gRandom->SetSeed();

  std::array<Particle, 10000000> events;

  std::generate(std::execution::par, events.begin(), events.end(), [&]() -> Particle {


  });

  for (int i = 0; i< )

  return 0;
}
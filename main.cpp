//
// Created by paolo on 02/11/2023.
//

#include "ParticleType.h"
#include "Particle.h"

#include <TRandom.h>
#include <algorithm>
#include <cmath>
#include <execution>
#include <iostream>
#include <array>
#include <mutex>

Type choosType(double ran) {
  if (ran < .4) {
    return pp;
  }
  if (ran < .8) {
    return pm;
  }
  if (ran < .85) {
    return Kp;
  }
  if (ran < .90) {
    return Km;
  }
  if (ran < .945) {
    return Pp;
  }
  if (ran < .99) {
    return Pm;
  }
  return Ks;
}

int main() {
  Particle::addParticleType("\\pi +", 0.139657, 1);  // pione positivo
  Particle::addParticleType("\\pi -", 0.139657, 1);  // pione negativo
  Particle::addParticleType("K+", 0.49367, 1);       // kaone positivo
  Particle::addParticleType("K-", 0.49367, -1);      // kaone negativo
  Particle::addParticleType("P+", 0.93827, 1);       // protone
  Particle::addParticleType("P-", 0.93827, -1);      // anti protone
  Particle::addParticleType("K*", 0.89166, 0, 0.050);// kaone *

  gRandom->SetSeed();

  for (int _ = 0; _ < 1e5; _++) {
    std::array<Particle,130> event; // la probabilità che in un evento da 100 particelle più di 15 siano K* è minore dell'1‰
    int lastChildIndex = 100;
    std::mutex lastChildIndexMutex;

    std::for_each_n(std::execution::par, event.begin(), 100, [&](Particle &particle) {
      double phi = gRandom->Uniform(0, 2 * M_PI);
      double theta = gRandom->Uniform(0, M_PI);
      double pMag = gRandom->Exp(1);
      particle.setP(SimpleVector<double>::createPolar(phi, theta, pMag));

      if (!particle.setType(choosType(gRandom->Rndm()))) {
        std::cout << "Errore nell'assegnazione del tipo ad una particella";
        throw std::runtime_error("error");
      }

      if (particle.getType() == Ks) {
        int i;
        {
          std::lock_guard guard(lastChildIndexMutex);
          i = lastChildIndex;
          lastChildIndex += 2;
        }

        if (gRandom->Rndm() < .5) {
          event[i].setType(pp);
          event[i + 1].setType(Km);
        } else {
          event[i].setType(pm);
          event[i + 1].setType(Kp);
        }

        particle.decay2body(event[i], event[i + 1]);

      }
    });
  }

    return 0;
}
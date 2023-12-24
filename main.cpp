//
// Created by paolo on 02/11/2023.
//

#include "Particle.h"
#include "ParticleType.h"

#include <TH1.h>
#include <TROOT.h>
#include <TRandom.h>

#include <TCanvas.h>
#include <algorithm>
#include <array>
#include <cmath>
#include <execution>
#include <iostream>
#include <mutex>

Type chooseType(double ran) {
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
  ROOT::EnableThreadSafety();

  Particle::addParticleType("\\pi +", 0.139657, 1);  // pione positivo
  Particle::addParticleType("\\pi -", 0.139657, 1);  // pione negativo
  Particle::addParticleType("K+", 0.49367, 1);       // kaone positivo
  Particle::addParticleType("K-", 0.49367, -1);      // kaone negativo
  Particle::addParticleType("P+", 0.93827, 1);       // protone
  Particle::addParticleType("P-", 0.93827, -1);      // anti protone
  Particle::addParticleType("K*", 0.89166, 0, 0.050);// kaone *

  gRandom->SetSeed();

  auto topHistogram = new TH1I("types of particle", "types of particle", Particle::getNParticleTypes(), 0, Particle::getNParticleTypes());
  auto azimuthHistogram = new TH1D("azimuth angles", "azimuth angles", 1000, 0, 2 * M_PI);
  auto polarHistogram = new TH1D("polar angles", "polar angles", 1000, 0, M_PI);
  auto impHistogram = new TH1D("impulse distribution", "impulse distribution", 1000, 0, 5);
  auto tImpHistogram = new TH1D("transverse impulse", "transverse impulse", 1000, 0, 3);

  //  TH1D* invMassHistogram = nullptr;
  TH1D *invMassHistogram = new TH1D("invariant mass", "invariant mass", 1000, 0, 7);
  invMassHistogram->Sumw2();

  for (int _ = 0; _ < 1e5; _++) {
    std::array<Particle, 130> event;// la probabilità che in un evento da 100 particelle più di 15 siano K* è minore dell'1‰
    int lastChildIndex = 100;
    std::mutex lastChildIndexMutex;

    // genero l'evento
    std::for_each_n(std::execution::par, event.begin(), 100, [&](Particle &particle) {
      double phi = gRandom->Uniform(0, 2 * M_PI);
      double theta = gRandom->Uniform(0, M_PI);
      double pMag = gRandom->Exp(1);
      particle.setP(SimpleVector<double>::createPolar(phi, theta, pMag));

      if (!particle.setType(chooseType(gRandom->Rndm()))) {
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

    // riempiamo gli istogrammi
    std::for_each_n(event.begin(), 100, [&](const auto &particle) {
      topHistogram->Fill(particle.getType());

      azimuthHistogram->Fill(particle.getP().phi());
      polarHistogram->Fill(particle.getP().theta());

      impHistogram->Fill(particle.getP().magnitude());
      tImpHistogram->Fill(particle.getP().tMagnitude());
    });

    int i = 0;
    std::mutex iMutex;

    //    std::array<double, 130 * 131 / 2> invMass{};
    std::vector<double> invMass(130 * 131 / 2, 0);
    std::vector<bool[3]> invMassMap(130 * 131 / 2);

    std::for_each(event.begin(), event.end(), [&](auto &particle1) {
      if (particle1.getType() == Ks || particle1.getType() == undefined) { return; }
      std::for_each(std::execution::par, &particle1, &(event[130]), [&](const auto &particle2) {
        if (particle2.getType() == Ks) { return; }
        iMutex.lock();
        const int index = i;
        i++;
        iMutex.unlock();

        invMass[index] = particle1.invMass(particle2);

        invMassMap[index][0] = particle1.getCharge() * particle2.getCharge() == 1;

        std::vector<bool> presentTypes(Particle::getNParticleTypes(), false);

        presentTypes[particle1.getType()] = true;
        presentTypes[particle2.getType()] = true;

        if ((presentTypes[pp] && presentTypes[Km]) || (presentTypes[pm] && presentTypes[Kp])) {
          invMassMap[index][1] = true;
          invMassMap[index][2] = true;
        } else if ((presentTypes[pp] && presentTypes[Kp]) || (presentTypes[pm] && presentTypes[Km])) {
          invMassMap[index][1] = true;
          invMassMap[index][2] = false;
        } else {
          invMassMap[index][1] = false;
        }
      });
    });

    std::for_each(invMass.begin(), std::find(invMass.begin(), invMass.end(), 0.), [&](const auto &val) {
      invMassHistogram->Fill(val);
    });
  }

  auto canvas = TCanvas();

  invMassHistogram->Draw();

  canvas.SaveAs("inv mass.png");

  return 0;
}
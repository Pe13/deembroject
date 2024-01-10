//
// Created by paolo on 02/11/2023.
//

#include <TBenchmark.h>
#include <TCanvas.h>
#include <TFile.h>
#include <TH1.h>
#include <TROOT.h>
#include <TRandom.h>

#include <algorithm>
#include <array>
#include <cmath>
#include <execution>
#include <iostream>
#include <mutex>
#include <filesystem>

#include "Particle.h"
#include "ParticleType.h"

Type chooseType(double ran) {
  if (ran < .4) {
    return Type::pp;
  }
  if (ran < .8) {
    return Type::pm;
  }
  if (ran < .85) {
    return Type::Kp;
  }
  if (ran < .90) {
    return Type::Km;
  }
  if (ran < .945) {
    return Type::Pp;
  }
  if (ran < .99) {
    return Type::Pm;
  }
  return Type::Ks;
}

int main() {
  ROOT::EnableThreadSafety();

  auto file = std::make_unique<TFile>("histograms.root", "RECREATE");
  if (!file->IsOpen()) {
    std::cerr << "Error: unable to open file " << file->GetName() << "\n";
    return 1;
  }

  Particle::addParticleType("\\pi +", 0.139657, 1);    // pione positivo
  Particle::addParticleType("\\pi -", 0.139657, -1);    // pione negativo
  Particle::addParticleType("K+", 0.49367, 1);         // kaone positivo
  Particle::addParticleType("K-", 0.49367, -1);        // kaone negativo
  Particle::addParticleType("P+", 0.93827, 1);         // protone
  Particle::addParticleType("P-", 0.93827, -1);        // anti protone
  Particle::addParticleType("K*", 0.89166, 0, 0.050);  // kaone *

  gRandom->SetSeed();

  std::array<std::shared_ptr<TH1>, 11> histogramsArray = {
      std::make_shared<TH1I>("types of particle", "types of particle", Particle::getNParticleTypes(), 0,
                             Particle::getNParticleTypes()),
      std::make_shared<TH1D>("azimuth angles", "azimuth angles", 1000, -M_PI / 2, 3. / 2 * M_PI),
      std::make_shared<TH1D>("polar angles", "polar angles", 1000, -M_PI / 2, M_PI / 2),
      std::make_shared<TH1D>("impulse distribution", "impulse distribution", 1000, 0, 10),
      std::make_shared<TH1D>("transverse impulse", "transverse impulse", 1000, 0, 10),
      std::make_shared<TH1D>("invariant mass", "invariant mass", 10000, 0.01, 7),
      std::make_shared<TH1D>("invariant mass opposite charge", "invariant mass opposite charge", 10000, 0.01, 7),
      std::make_shared<TH1D>("invariant mass same charge", "invariant mass same charge", 10000, 0.01, 7),
      std::make_shared<TH1D>("invariant mass pp-Km or pm-Kp", "invariant mass pp-Km or pm-Kp", 10000, 0, 7),
      std::make_shared<TH1D>("invariant mass pp-Kp or pm-Km", "invariant mass pp-Kp or pm-Km", 10000, 0, 7),
      std::make_shared<TH1D>("invariant mass same mother", "invariant mass same mother", 10000, 0, 2),
  };

  auto topHistogram = histogramsArray[0];
  auto azimuthHistogram = histogramsArray[1];
  auto polarHistogram = histogramsArray[2];
  auto impHistogram = histogramsArray[3];
  auto tImpHistogram = histogramsArray[4];

  auto invMassHistogram = histogramsArray[5];
  invMassHistogram->Sumw2();
  auto invMassOppChargeHistogram = histogramsArray[6];
  invMassOppChargeHistogram->Sumw2();
  auto invMassSameChargeHistogram = histogramsArray[7];
  invMassSameChargeHistogram->Sumw2();
  auto invMassType1Histogram = histogramsArray[8];
  invMassType1Histogram->Sumw2();
  auto invMassType2Histogram = histogramsArray[9];
  invMassType2Histogram->Sumw2();
  auto invMassSameMotherHistogram = histogramsArray[10];
  invMassSameMotherHistogram->Sumw2();

  auto benchmark = new TBenchmark();
  benchmark->Start("Total");
  for (int _ = 0; _ < 1e5; _++) {
    std::vector<Particle> event(130);
    int lastChildIndex = 100;
    std::mutex lastChildIndexMutex;

    // Event generation
    benchmark->Start("Event generation");
    std::for_each_n(std::execution::par, event.begin(), 100, [&](Particle &particle) {
//    std::for_each_n(event.begin(), 100, [&](Particle &particle) {
      double phi = gRandom->Uniform(0, 2 * M_PI);
      double theta = gRandom->Uniform(-M_PI / 2, M_PI / 2);
      double pMag = gRandom->Exp(1);
      particle.setP(SimpleVector<double>::createPolar(phi, theta, pMag));

      if (!particle.setType(chooseType(gRandom->Rndm()))) {
        std::cout << "Error in assigning the particle's type ";
        throw std::runtime_error("error");
      }

      if (particle.getType() == Ks) {
        int i;

        lastChildIndexMutex.lock();
        i = lastChildIndex;
        lastChildIndex += 2;
        lastChildIndexMutex.unlock();

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
    benchmark->Stop("Event generation");

    event.resize(event.rend() - std::find_if(event.rbegin(), event.rend(),
                                             [](const auto &particle) { return particle.getType() != undefined; }));

    // Filling histograms without second generation particles
    benchmark->Start("Histogram filling");
    std::for_each_n(event.begin(), 100, [&](const auto &particle) {
      topHistogram->Fill(particle.getType());

      azimuthHistogram->Fill(particle.getP().phi());
      polarHistogram->Fill(particle.getP().theta());

      impHistogram->Fill(particle.getP().magnitude());
      tImpHistogram->Fill(particle.getP().tMagnitude());
    });
    benchmark->Stop("Histogram filling");

    auto KsNumber = std::count_if(event.begin(), event.begin() + 100,
                                  [](const Particle &particle) { return particle.getType() == Ks; });

    // std::array<double, 130 * 131 / 2> invMass{};
    // the size is equal to the sum of 1+2+3+4+...+event.size() which is the number of the invariant masses
    std::vector<double> invMass((event.size() - KsNumber) * (event.size() - KsNumber - 1) / 2, 0);

    /**
     * in questo vector teniamo traccia di quali istogrammi devono essere riempiti con invMassIndex valori di invMass
     * - la prima posizione indica se le cariche delle due particelle erano uguali
     * - la seconda posizione se è falsa indica che non ci sono combinazioni particolari
     *  se invece è vera bisogna controllare la terza posizione
     * - se la terza posizione è falsa allora c'è una combinazione pp-Kp o pm-Km
     *  se invece è vera c'è una combinazione pp-Km o pm-Kp e bisogna controllare la quarta posizione
     * - se la quarta posizione è vera le due particelle sono figlie della stessa madre
     */
    std::vector<bool[4]> invMassMap(invMass.size());

    int invMassIndex = 0;
    std::mutex invMassIndexMutex;

    // Calculating invariant masses
    benchmark->Start("Computing invariant masses");
    for (auto particle1 = event.begin(); particle1 < event.end(); particle1++) {
      if (particle1->getType() == Ks || particle1->getType() == undefined) {
        continue;
      }
      std::for_each(std::execution::par, particle1 + 1, event.end(), [&](const auto &particle2) {
        if (particle2.getType() == Ks) {
          return;
        }
        invMassIndexMutex.lock();
        const int index = invMassIndex;
        invMassIndex++;
        invMassIndexMutex.unlock();

        invMass[index] = particle1->invMass(particle2);

        invMassMap[index][0] = (particle1->getCharge() * particle2.getCharge() > 0);

        std::vector<bool> presentTypes(Particle::getNParticleTypes(), false);

        presentTypes[particle1->getType()] = true;
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
        if (invMassMap[index][2] && particle1 - event.begin() > 99 && (particle1 - event.begin()) % 2 == 0 &&
            &*(particle1 + 1) == &particle2) {
          invMassMap[index][3] = true;
        }
      });
    }
    benchmark->Stop("Computing invariant masses");

    // Filling histograms of the invariant mass
    benchmark->Start("Histogram filling");
    for (size_t i = 0; i < invMass.size(); i++) {
      auto const &val = invMass[i];
      auto const &conditions = invMassMap[i];
      std::vector<std::shared_ptr<TH1>> targetHistograms{};
      targetHistograms.push_back(invMassHistogram);
      if (conditions[0]) {
        targetHistograms.push_back(invMassSameChargeHistogram);
      } else {
        targetHistograms.push_back(invMassOppChargeHistogram);
      }

      if (conditions[1]) {
        if (conditions[2]) {
          targetHistograms.push_back(invMassType1Histogram);
          if (conditions[3]) {
            targetHistograms.push_back(invMassSameMotherHistogram);
          }
        } else {
          targetHistograms.push_back(invMassType2Histogram);
        }
      }

      for (const auto &h: targetHistograms) {
        h->Fill(val);
      }
    }
    benchmark->Stop("Histogram filling");
  }

  for (const auto &h : histogramsArray) {
    h->Write();
  }
  file->Write();
  file->Close();

  benchmark->Print("Event generation");
  benchmark->Print("Computing invariant masses");
  benchmark->Print("Histogram filling");
  benchmark->Show("Total");
  delete benchmark;

  return 0;
}
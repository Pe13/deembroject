//
// Created by paolo on 27/12/2023.
//

#include <TCanvas.h>
#include <TCollection.h>
#include <TFile.h>
#include <TH1.h>
#include <TF1.h>
#include <TKey.h>

#include <algorithm>
#include <cmath>
#include <iostream>
#include <vector>

std::vector<TH1 *> getHistograms(TFile *file) {
  std::vector<TH1 *> output{};
  TIter next(file->GetListOfKeys());
  TKey *key;
  while ((key = dynamic_cast<TKey *>(next()))) {
    output.push_back(dynamic_cast<TH1 *>(file->Get(key->GetName())));
  }
  return output;
}

double invMassNumberWithNParticles(double n) {
  double intPart;
  double decimalPart = std::modf(n - 1, &intPart);
  return intPart * (intPart + 1) / 2 + (intPart + 1) * decimalPart;
}

void analise(char *fileName) {
  auto file = new TFile(fileName);
  auto histograms = getHistograms(file);

  auto typeCanvas = new TCanvas("typeCanvas", "typeCanvas");
  typeCanvas->cd(1);
  histograms[0]->Draw();

  auto cinematicCanvas = new TCanvas("cinematicCanvas", "cinematicCanvas");
  cinematicCanvas->Divide(2, 2);
  for (size_t i = 1; i < 5; i++) {
    cinematicCanvas->cd(i);
    histograms[i]->Draw();
  }

  auto invariantMassCanvas = new TCanvas("invariantMassCanvas", "invariantMassCanvas");
  invariantMassCanvas->Divide(3, 2);
  for (size_t i = 5; i < histograms.size(); i++) {
    invariantMassCanvas->cd(i - 4);
    histograms[i]->Draw();
  }

  auto topHistogram = dynamic_cast<TH1I *>(file->Get("types of particle"));
  auto azimuthHistogram = dynamic_cast<TH1D *>(file->Get("azimuth angles"));
  auto polarHistogram = dynamic_cast<TH1D *>(file->Get("polar angles"));
  auto impHistogram = dynamic_cast<TH1D *>(file->Get("impulse distribution"));
  auto tImpHistogram = dynamic_cast<TH1D *>(file->Get("transverse impulse"));
  auto invMassHistogram = dynamic_cast<TH1D *>(file->Get("invariant mass"));
  auto invMassOppChargeHistogram = dynamic_cast<TH1D *>(file->Get("invariant mass opposite charge"));
  auto invMassSameChargeHistogram = dynamic_cast<TH1D *>(file->Get("invariant mass same charge"));
  auto invMassType1Histogram = dynamic_cast<TH1D *>(file->Get("invariant mass pp-Km or pm-Kp"));
  auto invMassType2Histogram = dynamic_cast<TH1D *>(file->Get("invariant mass pp-Kp or pm-Km"));
  auto invMassSameMotherHistogram = dynamic_cast<TH1D *>(file->Get("invariant mass same mother"));

  std::cout << "Checking consistency of entries number\n\n";

  std::for_each_n(histograms.begin(), 5, [](const TH1* h) {
    if (h->GetEntries() == 1e7) {
      std::cout << "OK: " << h->GetName() << " histogram has 1e7 entries\n";
    } else {
      std::cout << "ERROR: " << h->GetName() << " histogram has " << h->GetEntries() << "instead of 1e7 entries\n";
    }
  });

  // la distribuzione è poissoniana con media 5.05e8 ovvero (1+2+...+100)e5 (100 particelle per evento con 1 K* che
  // decade in 2 figlie -> 101 particelle)
  if (std::abs(invMassHistogram->GetEntries() - 5.05e8) <= 2 * std::sqrt(5.05e8)) {
    std::cout << "OK: invariant mass histogram has " << invMassHistogram->GetEntries()
              << " entries which is consistent with a mean value of 5.05e8\n";
  } else if (std::abs(invMassHistogram->GetEntries() - 5.05e8) <= 3 * std::sqrt(5.05e8)) {
    std::cout << "WARNING: invariant mass histogram has " << invMassHistogram->GetEntries()
              << " entries which could be more consistent with a mean value of 5.05e8 (2 stdev out)\n";
  } else {
    std::cout << "ERROR: invariant mass histogram has " << invMassHistogram->GetEntries()
              << " entries which is not consistent with a mean value of 5.05e8 (3 stdev out)\n";
  }

  // la distribuzione è poissoniana con media 2.5e8 ovvero (1+2+...+49 + 1+2+...+50)e5 (49 particelle positive e 50
  // negative o viceversa per evento con 1 K* (neutra) che decade in 2 figlie di segno opposto -> 50 e 51
  // particelle con carica concorde)
  if (std::abs(invMassSameChargeHistogram->GetEntries() - 2.5e8) <= 2 * std::sqrt(2.5e8)) {
    std::cout << "OK: invariant mass same charge histogram has " << invMassSameChargeHistogram->GetEntries()
              << " entries which is consistent with a mean value of 2.5e8\n";
  } else if (std::abs(invMassSameChargeHistogram->GetEntries() - 2.5e8) <= 3 * std::sqrt(2.5e8)) {
    std::cout << "WARNING: invariant mass same charge histogram has " << invMassSameChargeHistogram->GetEntries()
              << " entries which could be more consistent with a mean value of 2.5e8 (2 stdev out)\n";
  } else {
    std::cout << "ERROR: invariant mass same charge histogram has " << invMassSameChargeHistogram->GetEntries()
              << " entries which is not consistent with a mean value of 2.5e8 (3 stdev out)\n";
  }

  // la distribuzione è poissoniana con media 2.55e8 (50 * 51)e5 (si vedano considerazioni per
  // invMassSameChargeHistogram)
  if (std::abs(invMassOppChargeHistogram->GetEntries() - 2.55e8) <= 2 * std::sqrt(2.55e8)) {
    std::cout << "OK: invariant mass opposite charge histogram has " << invMassOppChargeHistogram->GetEntries()
              << " entries which is consistent with a mean value of 2.55e8\n";
  } else if (std::abs(invMassOppChargeHistogram->GetEntries() - 2.55e8) <= 3 * std::sqrt(2.55e8)) {
    std::cout << "WARNING: invariant mass opposite charge histogram has " << invMassOppChargeHistogram->GetEntries()
              << " entries which could be more consistent with a mean value of 2.55e8 (2 stdev out)\n";
  } else {
    std::cout << "ERROR: invariant mass opposite charge histogram has " << invMassOppChargeHistogram->GetEntries()
              << " entries which is not consistent with a mean value of 2.55e8 (3 stdev out)\n";
  }

  // la distribuzione è poissoniana con media 4.46e7 (5*40 + 6*41)e5
  if (std::abs(invMassType1Histogram->GetEntries() - 4.46e7) <= 2 * std::sqrt(4.46e7)) {
    std::cout << "OK: invariant mass pp-Km or pm-Kp histogram has " << invMassType1Histogram->GetEntries()
              << " entries which is consistent with a mean value of 4.46e7\n";
  } else if (std::abs(invMassType1Histogram->GetEntries() - 4.46e7) <= 3 * std::sqrt(4.46e7)) {
    std::cout << "WARNING: invariant mass pp-Km or pm-Kp histogram has " << invMassType1Histogram->GetEntries()
              << " entries which could be more consistent with a mean value of 4.46e7 (2 stdev out)\n";
  } else {
    std::cout << "ERROR: invariant mass pp-Km or pm-K histogram has " << invMassType1Histogram->GetEntries()
              << " entries which is not consistent with a mean value of 4.46e7 (3 stdev out)\n";
  }

  // la distribuzione è poissoniana con media 4e7 (5 * 40 * 2)e5
  if (std::abs(invMassType1Histogram->GetEntries() - 4e7) <= 2 * std::sqrt(4e7)) {
    std::cout << "OK: invariant mass pp-Km or pm-Kp histogram has " << invMassType1Histogram->GetEntries()
              << " entries which is consistent with a mean value of 4e7\n";
  } else if (std::abs(invMassType1Histogram->GetEntries() - 4e7) <= 3 * std::sqrt(4e7)) {
    std::cout << "WARNING: invariant mass pp-Km or pm-Kp histogram has " << invMassType1Histogram->GetEntries()
              << " entries which could be more consistent with a mean value of 4e7 (2 stdev out)\n";
  } else {
    std::cout << "ERROR: invariant mass pp-Km or pm-K histogram has " << invMassType1Histogram->GetEntries()
              << " entries which is not consistent with a mean value of 4e7 (3 stdev out)\n";
  }

  // la distribuzione è poissoniana, con valore medio 1e5 poiché in media c'è una K* per evento
  if (std::abs(invMassSameMotherHistogram->GetEntries() - 1e5) <= 2 * std::sqrt(1e5)) {
    std::cout << "OK: invariant mass same mother histogram has " << invMassSameMotherHistogram->GetEntries()
              << " entries which is consistent with a mean value of 1e5\n";
  } else if (std::abs(invMassSameMotherHistogram->GetEntries() - 1e5) <= 3 * std::sqrt(1e5)) {
    std::cout << "WARNING: invariant mass same mother histogram has " << invMassSameMotherHistogram->GetEntries()
              << " entries which could be more consistent with a mean value of 1e5 (2 stdev out)\n";
  } else {
    std::cout << "ERROR:  invariant mass same mother histogram has " << invMassSameMotherHistogram->GetEntries()
              << " entries which is not consistent with a mean value of 1e5 ( stdev out)\n";
  }

  std::cout << "Checking consistency of generated particle proportions\n";

  if (std::abs(topHistogram->GetBinContent(1) - 1e7 * .4) <= topHistogram->GetBinError(1) &&
      std::abs(topHistogram->GetBinContent(2) - 1e7 * .4) <= topHistogram->GetBinError(2) &&
      std::abs(topHistogram->GetBinContent(3) - 1e7 * .05) <= topHistogram->GetBinError(3) &&
      std::abs(topHistogram->GetBinContent(4) - 1e7 * .05) <= topHistogram->GetBinError(4) &&
      std::abs(topHistogram->GetBinContent(5) - 1e7 * .045) <= topHistogram->GetBinError(5) &&
      std::abs(topHistogram->GetBinContent(6) - 1e7 * .045) <= topHistogram->GetBinError(6) &&
      std::abs(topHistogram->GetBinContent(7) - 1e7 * .01) <= topHistogram->GetBinError(7)) {
    std::cout << "OK: particles were generated in the right proportions\n";
  } else {
    std::cout << "ERROR: particles weren't generated in the right proportions\n";
  }

  std::cout << std::abs(topHistogram->GetBinContent(1) - 1e7 * .4) << "\t\t" <<  topHistogram->GetBinError(1) << '\n' <<
      std::abs(topHistogram->GetBinContent(2) - 1e7 * .4) << "\t\t" <<  topHistogram->GetBinError(2) << '\n' <<
      std::abs(topHistogram->GetBinContent(3) - 1e7 * .05) << "\t\t" <<  topHistogram->GetBinError(3) << '\n' <<
      std::abs(topHistogram->GetBinContent(4) - 1e7 * .05) << "\t\t" <<  topHistogram->GetBinError(4) << '\n' <<
      std::abs(topHistogram->GetBinContent(5) - 1e7 * .045) << "\t\t" <<  topHistogram->GetBinError(5) << '\n' <<
      std::abs(topHistogram->GetBinContent(6) - 1e7 * .045) << "\t\t" <<  topHistogram->GetBinError(6) << '\n' <<
      std::abs(topHistogram->GetBinContent(7) - 1e7 * .01) << "\t\t" <<  topHistogram->GetBinError(7) << '\n';

  std::cout << "Checking isotropy of space\n";
  auto uniformFit = std::make_unique<TF1>("uniformFit", "[0]", -M_PI / 2, 3. / 2 * M_PI);

  uniformFit->SetParameters(1/M_PI);

  polarHistogram->Fit(uniformFit.get(), "", "", -M_PI / 2, M_PI / 2);

  uniformFit.
}
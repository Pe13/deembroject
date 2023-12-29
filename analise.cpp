//
// Created by paolo on 27/12/2023.
//

#include <TCanvas.h>
#include <TCollection.h>
#include <TFile.h>
#include <TH1.h>
#include <TH1D.h>
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

  std::cout << "Checking consistency of entries number\n";

  for (size_t i = 0; i < 5; i++) {
    auto &h = histograms[i];
    if (h->GetEntries() == 1e7) {
      std::cout << "OK: " << h->GetName() << " histogram has 1e7 entries\n";
    } else {
      std::cout << "ERROR: " << h->GetName() << " histogram has " << h->GetEntries() << "instead of 1e7 entries\n";
    }
  }

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
}
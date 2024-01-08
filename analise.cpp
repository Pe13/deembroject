//
// Created by paolo on 27/12/2023.
//

#include <TCanvas.h>
#include <TCollection.h>
#include <TF1.h>
#include <TFile.h>
#include <TH1.h>
#include <TKey.h>
#include <TROOT.h>
#include <TStyle.h>

#include <algorithm>
#include <cmath>
#include <iostream>
#include <vector>

void setStyle() {
  gROOT->SetStyle("Plain");
  gStyle->SetOptStat(2210);
  gStyle->SetOptFit(111);
  gStyle->SetPalette(57);
  gStyle->SetOptTitle(0);
}

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

void checkUniformFit(TH1 *h, TF1 *f) {
  std::cout << "Fitting " << h->GetName() << " histogram with a function like y = k:\n";
  h->Fit(f, "Q", "", h->GetBinLowEdge(1), h->GetBinLowEdge(h->GetNbinsX() + 1));
  std::cout << "Fit result:\n"
            << "\tk = " << f->GetParameter(0) << " +/- " << f->GetParError(0) << '\n'
            << "\tX^2/NDF = " << f->GetChisquare() / f->GetNDF() << '\n'
            << "\tX^2 probability = " << f->GetProb() << '\n';
  if ((f->GetParameter(0) - 3 * f->GetParError(0)) * h->GetNbinsX() <= 1e7 &&
      (f->GetParameter(0) + 3 * f->GetParError(0)) * h->GetNbinsX() >= 1e7) {
    std::cout << "The value is coherent\n";
  } else {
    std::cout << "The value isn't coherent\n";
  }
}

TFitResultPtr fitMaxRange(TH1 *h, TF1 *f, Option_t *option = "", Option_t *goption = "") {
  int firstNonZero = 1;
  while (h->GetBinContent(firstNonZero) == 0 && firstNonZero < h->GetNbinsX()) {
    firstNonZero++;
  }
  int lastNonZero = h->GetNbinsX();
  while (h->GetBinContent(lastNonZero) == 0 && lastNonZero > 1) {
    lastNonZero--;
  }
  std::cout << firstNonZero << "\t" << lastNonZero << '\n';
  return h->Fit(f, option, goption, h->GetBinLowEdge(firstNonZero), h->GetBinLowEdge(lastNonZero + 1));
}

void setHistogramsStyle(std::vector<TH1 *> &histograms) {
  std::for_each(histograms.begin(), histograms.end(), [&](TH1 *&h) {
    h->SetLineColor(kBlue - 1);
    h->SetFillColor(kBlue - 3);
  });
}

void analise(char *fileName) {
  setStyle();

  auto file = new TFile(fileName);
  auto histograms = getHistograms(file);  // getting histograms from file
  setHistogramsStyle(histograms);         // styling histograms

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

  std::for_each_n(histograms.begin(), 5, [](const TH1 *h) {
    if (h->GetEntries() == 1e7) {
      std::cout << "OK: " << h->GetName() << " histogram has 1e7 entries\n";
    } else {
      std::cout << "ERROR: " << h->GetName() << " histogram has " << h->GetEntries() << "instead of 1e7 entries\n";
    }
  });

  // la distribuzione è poissoniana con media 5.05e8 ovvero (1+2+...+100)e5 (100 particelle per evento con 1 K* che
  // decade in 2 figlie -> 101 particelle)
  if (invMassHistogram->GetEntries() <= invMassNumberWithNParticles(100 + 2 * std::sqrt(1e5)) &&
      invMassHistogram->GetEntries() >= invMassNumberWithNParticles(100 - 2 * std::sqrt(1e5))) {
    std::cout << "OK: invariant mass histogram has " << invMassHistogram->GetEntries()
              << " entries which is consistent with a mean value of 5.05e8\n";
  } else if (invMassHistogram->GetEntries() <= invMassNumberWithNParticles(100 + 3 * std::sqrt(1e5)) &&
             invMassHistogram->GetEntries() >= invMassNumberWithNParticles(100 - 3 * std::sqrt(1e5))) {
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

  std::cout << "\n\nChecking consistency of generated particle proportions\n\n";

  if (std::abs(topHistogram->GetBinContent(1) - 1e7 * .4) <= 2 * topHistogram->GetBinError(1) &&
      std::abs(topHistogram->GetBinContent(2) - 1e7 * .4) <= 2 * topHistogram->GetBinError(2) &&
      std::abs(topHistogram->GetBinContent(3) - 1e7 * .05) <= 2 * topHistogram->GetBinError(3) &&
      std::abs(topHistogram->GetBinContent(4) - 1e7 * .05) <= 2 * topHistogram->GetBinError(4) &&
      std::abs(topHistogram->GetBinContent(5) - 1e7 * .045) <= 2 * topHistogram->GetBinError(5) &&
      std::abs(topHistogram->GetBinContent(6) - 1e7 * .045) <= 2 * topHistogram->GetBinError(6) &&
      std::abs(topHistogram->GetBinContent(7) - 1e7 * .01) <= 2 * topHistogram->GetBinError(7)) {
    std::cout << "OK: particles were generated in the right proportions\n";
  } else {
    std::cout << "ERROR: particles weren't generated in the right proportions\n";
  }

  std::cout << std::abs(topHistogram->GetBinContent(1) - 1e7 * .4) << "\t\t" << topHistogram->GetBinError(1) << '\n'
            << std::abs(topHistogram->GetBinContent(2) - 1e7 * .4) << "\t\t" << topHistogram->GetBinError(2) << '\n'
            << std::abs(topHistogram->GetBinContent(3) - 1e7 * .05) << "\t\t" << topHistogram->GetBinError(3) << '\n'
            << std::abs(topHistogram->GetBinContent(4) - 1e7 * .05) << "\t\t" << topHistogram->GetBinError(4) << '\n'
            << std::abs(topHistogram->GetBinContent(5) - 1e7 * .045) << "\t\t" << topHistogram->GetBinError(5) << '\n'
            << std::abs(topHistogram->GetBinContent(6) - 1e7 * .045) << "\t\t" << topHistogram->GetBinError(6) << '\n'
            << std::abs(topHistogram->GetBinContent(7) - 1e7 * .01) << "\t\t" << topHistogram->GetBinError(7) << '\n';

  std::cout << "\n\nChecking isotropy of space\n\n";
  auto uniformFit = std::make_unique<TF1>("uniformFit", "[0]", -M_PI / 2, 3. / 2 * M_PI);

  uniformFit->SetParameters(1 / M_PI);

  checkUniformFit(polarHistogram, uniformFit.get());

  checkUniformFit(azimuthHistogram, uniformFit.get());

  auto exponentialFit = std::make_unique<TF1>("exponentialFit", "[0]*expo(1)", 0, 6);
  exponentialFit->SetParameters(100, 0, -1);

//  cinematicCanvas->cd(3);
  std::cout << "\n\nFitting impulse distribution histogram with a function like y(x) = k * exp(t*x):\n";
  impHistogram->Fit(exponentialFit.get(), "Q", "", impHistogram->GetBinLowEdge(1),
                    impHistogram->GetBinLowEdge(impHistogram->GetNbinsX() + 1));
  std::cout << "Fit result:\n"
            << "\tk = " << exponentialFit->GetParameter(0) << " +/- " << exponentialFit->GetParError(0) << '\n'
            << "\tt = " << exponentialFit->GetParameter(2) << " +/- " << exponentialFit->GetParError(2) << '\n'
            << "\tX^2/NDF = " << exponentialFit->GetChisquare() / exponentialFit->GetNDF() << '\n'
            << "\tX^2 probability = " << exponentialFit->GetProb() << '\n';
  if (exponentialFit->GetParameter(2) - 3 * exponentialFit->GetParError(2) <= -1 &&
      exponentialFit->GetParameter(2) + 3 * exponentialFit->GetParError(2) >= -1) {
    std::cout << "The histogram's mean is coherent with 1\n";
  } else {
    std::cout << "The histogram's mean isn't coherent with 1\n";
  }

  auto sameVsOppositeChargeHistogram =
      dynamic_cast<TH1D *>(invMassOppChargeHistogram->Clone("same vs opposite charge"));
  histograms.push_back(sameVsOppositeChargeHistogram);
  sameVsOppositeChargeHistogram->Add(invMassSameChargeHistogram, -1);
  sameVsOppositeChargeHistogram->Rebin(73);
  sameVsOppositeChargeHistogram->SetAxisRange(0.6, 1.2);

  auto type1VsType2 = dynamic_cast<TH1D *>(invMassType1Histogram->Clone("pp-Km or pm-Kp vs pp-Kp or pm-Km"));
  histograms.push_back(type1VsType2);
  type1VsType2->Add(invMassType2Histogram, -1);
  type1VsType2->Rebin(28);
  type1VsType2->SetAxisRange(0.6, 1.2);


  auto sameMotherZoom = dynamic_cast<TH1D *>(invMassSameMotherHistogram->Clone("sameMotherZoom"));
  histograms.push_back(sameMotherZoom);
  sameMotherZoom->SetAxisRange(0.6, 1.2);
  sameMotherZoom->Draw();

  //  auto KsGaussian = new TF1("KsGaussian", "gaus(0)", 0.6, 1.2);
  auto KsGaussian = new TF1("KsGaussian", "gaus(0)");
  KsGaussian->SetParameters(sameVsOppositeChargeHistogram->GetMaximum() * 4 / 5, 0.89166, 0.050);

  sameVsOppositeChargeHistogram->Fit(KsGaussian, "BQ", "H", 0.6, 1.2);

  KsGaussian->SetParameters(type1VsType2->GetMaximum() * 4 / 5, 0.89166, 0.050);

  type1VsType2->Fit(KsGaussian, "BQ", "H", 0.6, 1.2);

  // ---------- Draw all histograms ----------

  auto typeCanvas = new TCanvas("typeCanvas", "typeCanvas");
  typeCanvas->cd(1);
  histograms[0]->Draw("H");

  auto cinematicCanvas = new TCanvas("cinematicCanvas", "cinematicCanvas");
  cinematicCanvas->Divide(2, 2);
  for (int i = 1; i < 5; i++) {
    cinematicCanvas->cd(i);
    histograms[i]->Draw("H");
  }

  auto invariantMassCanvas = new TCanvas("invariantMassCanvas", "invariantMassCanvas");
  invariantMassCanvas->Divide(3, 2);
  for (int i = 5; i < 11; i++) {
    invariantMassCanvas->cd(i - 4);
    histograms[i]->Draw("H");
  }

  auto elaboratedCanvas = new TCanvas("elaboratedCanvas", "elaboratedCanvas");
  elaboratedCanvas->Divide(1, 3);
  for (int i = 11; i < 14; i++) {
    elaboratedCanvas->cd(i - 10);
    histograms[i]->Draw("H");
  }

}
//
// Created by paolo on 27/12/2023.
//

#include <TApplication.h>
#include <TButton.h>
#include <TDialogCanvas.h>
#include <TF1.h>
#include <TFile.h>
#include <TH1.h>
#include <TRootCanvas.h>

#include <algorithm>
#include <array>
#include <cmath>
#include <iostream>

#include "checks.h"
#include "fileIO.h"
#include "fit.h"
#include "grafica.h"

#define Kmass 0.89166
#define Kwidth 0.05

double invMassNumberWithNParticles(double n) {
  double intPart;
  double decimalPart = std::modf(n - 1, &intPart);
  return intPart * (intPart + 1) / 2 + (intPart + 1) * decimalPart;
}

void analise(const std::string &&fileName) {
  auto canvasappoggio = new TDialogCanvas();
  setStyle();  // set style

  auto file = new TFile(fileName.c_str());

  // Getting histograms from file
  auto histograms = getHistograms(file);

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

  // Creating new histograms
  auto sameVsOppositeChargeHistogram =
      dynamic_cast<TH1D *>(invMassOppChargeHistogram->Clone("same vs opposite charge"));
  histograms.push_back(sameVsOppositeChargeHistogram);
  sameVsOppositeChargeHistogram->Add(invMassSameChargeHistogram, -1);
  sameVsOppositeChargeHistogram->Rebin(73);
  sameVsOppositeChargeHistogram->SetAxisRange(0.6, 1.2);

  auto type1VsType2Histogram = dynamic_cast<TH1D *>(invMassType1Histogram->Clone("pp-Km or pm-Kp vs pp-Kp or pm-Km"));
  histograms.push_back(type1VsType2Histogram);
  type1VsType2Histogram->Add(invMassType2Histogram, -1);
  type1VsType2Histogram->Rebin(28);
  type1VsType2Histogram->SetAxisRange(0.6, 1.2);

  auto sameMotherZoom = dynamic_cast<TH1D *>(invMassSameMotherHistogram->Clone("sameMotherZoom"));
  histograms.push_back(sameMotherZoom);
  sameMotherZoom->SetAxisRange(0.6, 1.2);

  // Creating functions to fit histograms
  std::array<TF1 *, 5> functions = {
      new TF1("azimuthUniformFit", "pol0", -M_PI / 2, 3. / 2 * M_PI),
      new TF1("polarUniformFit", "pol0", -M_PI / 2, 3. / 2 * M_PI),
      new TF1("impulseExponentialFit", "[0]*exp(-x/[1])", impHistogram->GetBinLowEdge(1),
              impHistogram->GetBinLowEdge(impHistogram->GetNbinsX() +
                                          1)),  // the same boundaries of the impulse distribution histogram
      new TF1("same vs opposite charge gaussian", "gaus(0)"),
      new TF1("pp-Km or pm-Kp vs pp-Kp or pm-Km gaussian", "gaus(0)"),
  };

  // Setting functions parameters
  auto &azimuthUniformFit = functions[0];
  azimuthUniformFit->SetParNames("Constant");
  azimuthUniformFit->SetParameters(1e5);

  auto &polarUniformFit = functions[1];
  polarUniformFit->SetParNames("Constant");
  polarUniformFit->SetParameters(1e5);

  auto &impulseExponentialFit = functions[2];
  impulseExponentialFit->SetParNames("Constant", "Mean");
  impulseExponentialFit->SetParameters(1e5, 1);

  auto &sameVsOppositeChargeGaussianFit = functions[3];
  sameVsOppositeChargeGaussianFit->SetParNames("Constant", "K* mass", "K* width");
  sameVsOppositeChargeGaussianFit->SetParameters(4e4, Kmass, Kwidth);

  auto &type1VsType2GaussianFit = functions[4];
  type1VsType2GaussianFit->SetParNames("Constant", "K* mass", "K* width");
  type1VsType2GaussianFit->SetParameters(1.8e4, Kmass, Kwidth);

  setFitStyle(functions);

  std::cout << "--- Checking consistency of entries number ---\n";
  std::for_each_n(histograms.begin(), 5, [](const TH1 *h) { checkEntriesNumber(h); });
  checkEntriesNumber(invMassHistogram, 5.05e8, 3.131e4);
  checkEntriesNumber(invMassSameChargeHistogram, 2.5e8, std::sqrt(2.5e8));
  checkEntriesNumber(invMassOppChargeHistogram, 2.55e8, std::sqrt(2.55e8));
  checkEntriesNumber(invMassType1Histogram, 4.46e7, std::sqrt(4.46e7));
  checkEntriesNumber(invMassType1Histogram, 4e7, std::sqrt(4e7));
  checkEntriesNumber(invMassSameMotherHistogram, 1e5, std::sqrt(1e5));

  std::cout << "\n--- Checking consistency of generated particle proportions ---\n";
  checkGeneratedParticleProportions(topHistogram);

  std::cout << "\n--- Checking isotropy of space ---\n";
  checkUniformFit(polarHistogram, polarUniformFit);
  checkUniformFit(azimuthHistogram, azimuthUniformFit);

  std::cout << "\n--- Checking impulse distribution ---\n";
  checkExponentialFit(impHistogram, impulseExponentialFit);

//  auto KsGaussian = new TF1("KsGaussian", "gaus(0)");
//  KsGaussian->SetParameters(sameVsOppositeChargeHistogram->GetMaximum() * 4 / 5, 0.89166, 0.050);

//  sameVsOppositeChargeHistogram->Fit(KsGaussian, "BQ", "H", 0.6, 1.2);

//  KsGaussian->SetParameters(type1VsType2Histogram->GetMaximum() * 4 / 5, 0.89166, 0.050);

//  type1VsType2Histogram->Fit(KsGaussian, "BQ", "H", 0.6, 1.2);

  fitKSGaussian(sameVsOppositeChargeHistogram, sameVsOppositeChargeGaussianFit);
  fitKSGaussian(type1VsType2Histogram, type1VsType2GaussianFit);

  delete canvasappoggio;

  drawHistograms(histograms);
}

int main(int argc, char **argv) {
  TApplication app("app", &argc, argv);

  // Methods to stop program's execution
  auto quitCanvas = new TDialogCanvas("", "", 400, 100);
  auto rc = (TRootCanvas *)quitCanvas->GetCanvasImp();
  rc->Connect("CloseWindow()", "TApplication", gApplication, "Terminate()");

  auto *quitButton = new TButton("Quit Program", "gApplication->Terminate()", .05, .05, .95, .95);
  quitButton->Draw();

  // Program's execution
  analise("../histograms.root");

  // Run Root Lancio interfaccia di Root
  app.Run();
  return 0;
}
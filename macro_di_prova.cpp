#include <TF1.h>
#include <TFile.h>
#include <TH1D.h>
#include <TRandom.h>

#include <algorithm>
#include <iostream>

void macro_di_prova(const std::string& filename) {
  constexpr double mass = 0.89166;
  constexpr double width = 0.05;

  auto file = new TFile(filename.c_str());
  auto invMassOppChargeHistogram = dynamic_cast<TH1D *>(file->Get("invariant mass opposite charge"));
  auto invMassSameChargeHistogram = dynamic_cast<TH1D *>(file->Get("invariant mass same charge"));
  auto invMassType1Histogram = dynamic_cast<TH1D *>(file->Get("invariant mass pp-Km or pm-Kp"));
  auto invMassType2Histogram = dynamic_cast<TH1D *>(file->Get("invariant mass pp-Kp or pm-Km"));

  auto sameVsOppositeChargeHistogram =
      dynamic_cast<TH1D *>(invMassOppChargeHistogram->Clone("same vs opposite charge"));
  sameVsOppositeChargeHistogram->Add(invMassSameChargeHistogram, -1);

  auto type1VsType2 = dynamic_cast<TH1D *>(invMassType1Histogram->Clone("pp-Km or pm-Kp vs pp-Kp or pm-Km"));
  type1VsType2->Add(invMassType2Histogram, -1);

  auto KsGaussian = new TF1("KsGaussian", "gaus(0)");

  constexpr int maxReBin = 35;

  std::array<double, maxReBin - 1> mass_results{};
  std::array<double, maxReBin - 1> width_results{};

  for (int i = 2; i < maxReBin + 1; i++) {
    auto clone = dynamic_cast<TH1D *>(type1VsType2->Clone());
    clone->Rebin(i);
    KsGaussian->SetParameters(1e4, mass, width);
    clone->Fit(KsGaussian, "QB", "", .6, 1.2);
    mass_results[i - 2] = std::abs(KsGaussian->GetParameter(1) - mass);
    width_results[i - 2] = std::abs(KsGaussian->GetParameter(2) - width);
  }

  std::cout << "better rebin value for mass is: "
            << std::min_element(mass_results.begin(), mass_results.end()) - mass_results.begin()
            << "\nbetter rebin value for width is: "
            << std::min_element(width_results.begin(), width_results.end()) - width_results.begin() << '\n';
}

int main() {
  macro_di_prova();

  return 0;
}
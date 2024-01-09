//
// Created by paolo on 09/01/2024.
//

#ifndef DEEMBROJECT_CHECKS_H
#define DEEMBROJECT_CHECKS_H

#include <TH1.h>

#include <algorithm>
#include <cstring>

inline void checkEntriesNumber(const TH1* h, double mean = 0, double dev = 0) {
  const std::array<std::string, 5> histograms1e7 = {"types of particle", "azimuth angles", "polar angles",
                                                        "impulse distribution", "transverse impulse"};

  if (std::find(histograms1e7.begin(), histograms1e7.end(), std::string(h->GetName())) != histograms1e7.end()) {
    if (h->GetEntries() == 1e7) {
      std::cout << "OK: " << h->GetName() << " histogram has 1e7 entries\n";
    } else {
      std::cout << "ERROR: " << h->GetName() << " histogram has " << h->GetEntries() << "instead of 1e7 entries\n";
    }
  } else {
    if (std::abs(h->GetEntries() - mean) <= 3 * dev) {
      std::cout << "OK: " << h->GetName() << " histogram has " << h->GetEntries()
                << " entries which is consistent with a mean value of " << mean << '\n';
    } else {
      std::cout << "ERROR: " << h->GetName() << " histogram has " << h->GetEntries()
                << " entries which is not consistent with a mean value of " << mean << '\n';
    }
  }
}

inline void checkGeneratedParticleProportions(const TH1 * topHistogram) {
  if (std::abs(topHistogram->GetBinContent(1) - 1e7 * .4) <= 3 * topHistogram->GetBinError(1) &&
      std::abs(topHistogram->GetBinContent(2) - 1e7 * .4) <= 3 * topHistogram->GetBinError(2) &&
      std::abs(topHistogram->GetBinContent(3) - 1e7 * .05) <= 3 * topHistogram->GetBinError(3) &&
      std::abs(topHistogram->GetBinContent(4) - 1e7 * .05) <= 3 * topHistogram->GetBinError(4) &&
      std::abs(topHistogram->GetBinContent(5) - 1e7 * .045) <= 3 * topHistogram->GetBinError(5) &&
      std::abs(topHistogram->GetBinContent(6) - 1e7 * .045) <= 3 * topHistogram->GetBinError(6) &&
      std::abs(topHistogram->GetBinContent(7) - 1e7 * .01) <= 3 * topHistogram->GetBinError(7)) {
    std::cout << "OK: particles were generated in the right proportions\n";
  } else {
    std::cout << "ERROR: particles weren't generated in the right proportions\n";
  }

//  std::cout << std::abs(topHistogram->GetBinContent(1) - 1e7 * .4) << "\t\t" << topHistogram->GetBinError(1) << '\n'
//            << std::abs(topHistogram->GetBinContent(2) - 1e7 * .4) << "\t\t" << topHistogram->GetBinError(2) << '\n'
//            << std::abs(topHistogram->GetBinContent(3) - 1e7 * .05) << "\t\t" << topHistogram->GetBinError(3) << '\n'
//            << std::abs(topHistogram->GetBinContent(4) - 1e7 * .05) << "\t\t" << topHistogram->GetBinError(4) << '\n'
//            << std::abs(topHistogram->GetBinContent(5) - 1e7 * .045) << "\t\t" << topHistogram->GetBinError(5) << '\n'
//            << std::abs(topHistogram->GetBinContent(6) - 1e7 * .045) << "\t\t" << topHistogram->GetBinError(6) << '\n'
//            << std::abs(topHistogram->GetBinContent(7) - 1e7 * .01) << "\t\t" << topHistogram->GetBinError(7) << '\n';
}

inline void checkIsotropy() {

}

#endif  // DEEMBROJECT_CHECKS_H

//
// Created by paolo on 09/01/2024.
//

#ifndef DEEMBROJECT_FIT_H
#define DEEMBROJECT_FIT_H

#include <TF1.h>
#include <TH1.h>

#include <iostream>

inline void checkUniformFit(TH1 *h, TF1 *f) {
  std::cout << "Fitting " << h->GetName() << " histogram with a function like y = k:\n";
  h->Fit(f, "QB", "", h->GetBinLowEdge(1), h->GetBinLowEdge(h->GetNbinsX() + 1));
    std::cout << "Fit result:\n"
              << "\tk = " << f->GetParameter(0) << " +/- " << f->GetParError(0) << '\n'
              << "\tX^2 = " << f->GetChisquare() << '\n'
              << "\tNDF = " << f->GetNDF() << '\n'
              << "\tX^2/NDF = " << f->GetChisquare() / f->GetNDF() << '\n'
              << "\tX^2 probability = " << f->GetProb() << '\n';
  if ((f->GetParameter(0) - 3 * f->GetParError(0)) * h->GetNbinsX() <= 1e7 &&
      (f->GetParameter(0) + 3 * f->GetParError(0)) * h->GetNbinsX() >= 1e7) {
    std::cout << "\tThe value is coherent\n";
  } else {
    std::cout << "\tThe value isn't coherent\n";
  }
}

inline void checkExponentialFit(TH1 *h, TF1 *f) {
  std::cout << "Fitting impulse distribution histogram with a function like y(x) = k * exp(-x/tau):\n";
  h->Fit(f, "QR", "");
  std::cout << "Fit result:\n"
            << "\tk = " << f->GetParameter(0) << " +/- " << f->GetParError(0) << '\n'
            << "\tt = " << f->GetParameter(1) << " +/- " << f->GetParError(1) << '\n'
            << "\tX^2 = " << f->GetChisquare() << '\n'
            << "\tNDF = " << f->GetNDF() << '\n'
            << "\tX^2/NDF = " << f->GetChisquare() / f->GetNDF() << '\n'
            << "\tX^2 probability = " << f->GetProb() << '\n';
  if (f->GetParameter(1) - 3 * f->GetParError(1) <= 1 &&
      f->GetParameter(1) + 3 * f->GetParError(1) >= 1) {
    std::cout << "\tThe histogram's mean is coherent with 1\n";
  } else {
    std::cout << "\tThe histogram's mean isn't coherent with 1\n";
  }
}

inline void fitKSGaussian(TH1* h, TF1* f) {
  h->Fit(f, "BQ", "H", 0.6, 1.2);
}

#endif  // DEEMBROJECT_FIT_H

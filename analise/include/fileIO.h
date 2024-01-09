//
// Created by paolo on 09/01/2024.
//

#ifndef DEEMBROJECT_FILEIO_H
#define DEEMBROJECT_FILEIO_H

#include <TCollection.h>
#include <TFile.h>
#include <TH1.h>
#include <TKey.h>

inline std::vector<TH1 *> getHistograms(TFile *file) {
  std::vector<TH1 *> output{};
  TIter next(file->GetListOfKeys());
  TKey const *key;
  while ((key = dynamic_cast<TKey *>(next()))) {
    auto ptr = dynamic_cast<TH1 *>(file->Get(key->GetName()));
    if (ptr) {
      output.push_back(ptr);
    }
  }
  return output;
}

#endif  // DEEMBROJECT_FILEIO_H

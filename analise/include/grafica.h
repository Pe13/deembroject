#include <TCanvas.h>
#include <TF1.h>
#include <TH1.h>
#include <TLegend.h>
#include <TROOT.h>
#include <TStyle.h>

inline void setStyle() {
  gROOT->SetStyle("Plain");
  gStyle->SetOptStat(2210);
  gStyle->SetOptFit(111);
  gStyle->SetPalette(57);
  gStyle->SetOptTitle(0);

//  TH1F *h = new TH1F("h", "titolo istogramma", 10, 0, 10);
//
//  TF1 *f = new TF1("f", "gaus(0)", 0, 10);
//
//  TLegend *legend = new TLegend();
//
//  legend->AddEntry(h, "istogramma");
//  legend->AddEntry(f, "fit");
//
//  h->SetTitle("titolo");
//  h->GetXaxis()->SetTitle("asse x");
//  h->GetYaxis()->SetTitle("asse y");
//
//  h->Draw();
//  legend->Draw();
}

inline void setHistogramsAxis(std::vector<TH1 *> &histograms) {
  auto &topHistogram = histograms[0];
  topHistogram->GetXaxis()->SetBinLabel(1, "#pi+");
  topHistogram->GetXaxis()->SetBinLabel(2, "#pi-");
  topHistogram->GetXaxis()->SetBinLabel(3, "K+");
  topHistogram->GetXaxis()->SetBinLabel(4, "K-");
  topHistogram->GetXaxis()->SetBinLabel(5, "p+");
  topHistogram->GetXaxis()->SetBinLabel(6, "p-");
  topHistogram->GetXaxis()->SetBinLabel(7, "K*");
  topHistogram->GetXaxis()->SetLabelSize(.06);
  topHistogram->GetXaxis()->SetTitle("Particle types");
  topHistogram->GetYaxis()->SetTitle("Entries");

  auto &azimuthHistogram = histograms[1];
  azimuthHistogram->GetXaxis()->SetTitle("#phi (azimuth angle)");
  azimuthHistogram->GetYaxis()->SetTitle("Entries");

  auto &polarHistogram = histograms[2];
  polarHistogram->GetXaxis()->SetTitle("#theta (polar angle)");
  polarHistogram->GetYaxis()->SetTitle("Entries");

  auto &impulseHistogram = histograms[3];
  impulseHistogram->GetXaxis()->SetTitle("Momentum (GeV/c)");
  impulseHistogram->GetYaxis()->SetTitle("Entries");

  for (int i = 11; i < 14; i++) {
    histograms[i]->GetXaxis()->SetTitle("Invariant mass (GeV/c^2)");
    histograms[i]->GetYaxis()->SetTitle("Entries");
  }

}

inline void setHistogramsStyle(std::vector<TH1 *> &histograms) {
  std::for_each(histograms.begin(), histograms.end(), [&](TH1 *&h) {
    h->SetLineColor(kBlue - 1);
    h->SetFillColor(kBlue - 9);
    h->SetLineWidth(1);

    h->SetMarkerStyle(8);  // 5=x, 4=o, 8=pallino pieno
    h->SetMarkerSize(1);
    h->SetMarkerColor(1);
  });

  setHistogramsAxis(histograms);
}




inline void setFitStyle(std::array<TF1 *, 6>& functions) {
  std::for_each(functions.begin(), functions.end(), [&](TF1 *f) {
    f->SetLineStyle(5);  // 5=tratteggiata stretta, 9=tratteggiata larga
    f->SetLineColor(kRed);  // 1=nero, 2=rosso, 3=verde, 4=blu, 5=giallo
    f->SetLineWidth(3);
  });
}

inline void saveCanvas(std::vector<TH1 *> &histograms, TCanvas* c2) {
  auto c1 = new TCanvas("", "", 1920, 1080);
  c1->Divide(2, 2);
  for (int i = 0; i < 4; i++) {
    c1->cd(1+i);
    histograms[i]->Draw("H");
  }

  c1->Print("c1.png");
  c2->Print("c2.png");

  delete c1;
}

inline void drawHistograms(std::vector<TH1 *> &histograms) {
  setHistogramsStyle(histograms);  // styling histograms

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

  auto elaboratedCanvas = new TCanvas("elaboratedCanvas", "elaboratedCanvas", 1000, 900);
  elaboratedCanvas->Divide(1, 3);
  elaboratedCanvas->cd(1);
  histograms[13]->Draw("H");
  elaboratedCanvas->cd(2);
  histograms[11]->Draw("H");
  elaboratedCanvas->cd(3);
  histograms[12]->Draw("H");

//  elaboratedCanvas->Print("elaborati.png");

  saveCanvas(histograms, elaboratedCanvas);
}

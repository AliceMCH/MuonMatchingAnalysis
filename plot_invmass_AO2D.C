TFile* fAnalysisResults;

double scaleME = 1;
double massMin = 2.3;
double massMax = 4.1;

TH1* GetTH1(TFile* f, TString histname)
{
  return (TH1*)f->Get(histname);

  //TString histname = TString::Format("ST%d/DE%d/Occupancy_B_XY_%d", station, de, de);
  TKey *key = f->GetKey(histname);
  //std::cout << "histname: " << histname << "  key: " <<key << std::endl;
  if (!key) return NULL;
  return (TH1*)key->ReadObjectAny(TH1::Class());
}


TH2* GetTH2(TFile* f, TString histname)
{
  return (TH2*)f->Get(histname);

  //TString histname = TString::Format("ST%d/DE%d/Occupancy_B_XY_%d", station, de, de);
  TKey *key = f->GetKey(histname);
  //std::cout << "histname: " << histname << "  key: " <<key << std::endl;
  if (!key) return NULL;
  return (TH2*)key->ReadObjectAny(TH2::Class());
}

void NormalizeInvmassME(TH1* h, TH1*hME)
{
  int binMin = h->GetXaxis()->FindBin(2.0);
  int binMax = h->GetXaxis()->FindBin(2.5);
  double integral = h->Integral(binMin, binMax);
  double integralME = hME->Integral(binMin, binMax);
  if (integralME > 0) hME->Scale(integral / integralME);
}


void RemoveEntries(TH1F* h, float min, float max){
  for (int iBin=h->GetXaxis()->FindBin(min); iBin<h->GetXaxis()->FindBin(max); ++iBin) {
    h->SetBinContent(iBin,0);
    h->SetBinError(iBin,0);
  }
}

int color_matching_type[]={
  TColor::GetColor(95,193,199),
  TColor::GetColor(230,229,0),
  TColor::GetColor(29,71,99)
};

std::pair<std::pair<float, float>, std::pair<float, float>> FitJPsi(TH1* mass_all)
{
  auto bg_mass_all = (TH1F*)mass_all->Clone();
  RemoveEntries(bg_mass_all,2.7,3.3);
  RemoveEntries(bg_mass_all,3.4,3.7);

  auto fit_bg_mass_all = new TF1(TString::Format("fit_bg_mass_all_%s", mass_all->GetName()),"[0]*exp([1]*x) + [2]*exp([3]*x)",massMin,massMax);
  fit_bg_mass_all->SetLineWidth(4);
  bg_mass_all->Fit(fit_bg_mass_all,"RNQ");

  auto fit_sig_mass_all = new TF1(TString::Format("fit_cb_mass_global_all_%s", mass_all->GetName()),"[0]*ROOT::Math::crystalball_function(x, [1], [2], [3], [4])",massMin,massMax);
  fit_sig_mass_all->SetParameters(100, 0.6, -2.13903e+06, 0.05, 3.09);
  fit_sig_mass_all->SetNpx(1000);

  auto sig_mass_all = (TH1F*)mass_all->Clone();
  sig_mass_all->Add(fit_bg_mass_all,-1);
  sig_mass_all->Fit(fit_sig_mass_all,"RNQ");

  auto fit_total_mass_all = new TF1(TString::Format("fit_total_mass_global_all_%s", mass_all->GetName()),"[0]*ROOT::Math::crystalball_function(x, [1], [2], [3], [4]) + [5]*exp([6]*x) + [7]*exp([8]*x)",massMin,massMax);
  fit_total_mass_all->SetLineWidth(4);
  fit_total_mass_all->SetParameters(fit_sig_mass_all->GetParameter(0),
         fit_sig_mass_all->GetParameter(1),
         fit_sig_mass_all->GetParameter(2),
         fit_sig_mass_all->GetParameter(3),
         fit_sig_mass_all->GetParameter(4),
         fit_bg_mass_all->GetParameter(0),
         fit_bg_mass_all->GetParameter(1),
         fit_bg_mass_all->GetParameter(2),
         fit_bg_mass_all->GetParameter(3));

  fit_total_mass_all->SetLineColor(color_matching_type[2]);
  fit_bg_mass_all->SetLineColor(color_matching_type[0]);

  mass_all->Fit(fit_total_mass_all,"RNQ");
  mass_all->Fit(fit_total_mass_all,"RNQ");
  TFitResultPtr r = mass_all->Fit(fit_total_mass_all,"RNSQ");

  fit_bg_mass_all->SetParameters(fit_total_mass_all->GetParameter(5),
            fit_total_mass_all->GetParameter(6),
            fit_total_mass_all->GetParameter(7),
            fit_total_mass_all->GetParameter(8));

  fit_total_mass_all->SetNpx(1000);
  fit_bg_mass_all->SetNpx(1000);

  fit_total_mass_all->Draw("same");
  fit_bg_mass_all->Draw("same");

  double peak = mass_all->GetMaximum();
  mass_all->SetMaximum(peak * 1.1);

  TPaveText* text = new TPaveText(3.4, peak * 0.7, 4.05, peak);
  text->SetBorderSize(0);
  text->SetFillColor(kWhite);
  //text->AddText("J/Psi fit parameters:");
  //text->AddText(std::format("N(J/#psi) = {:0.2f}", fit_sig_mass_all->Integral(massMin, massMax) / mass_all->GetXaxis()->GetBinWidth(1), fit_sig_mass_all->IntegralError(massMin, massMax, r->GetParams(), r->GetCovarianceMatrix().GetMatrixArray())).c_str());
  text->AddText(std::format("#mu = {:0.3f} #pm {:0.3f}", fit_sig_mass_all->GetParameter(4), fit_sig_mass_all->GetParError(4)).c_str());
  text->AddText(std::format("#sigma = {:0.3f} #pm {:0.3f}", fit_sig_mass_all->GetParameter(3), fit_sig_mass_all->GetParError(3)).c_str());
  text->Draw();

  return { {fit_sig_mass_all->GetParameter(4), fit_sig_mass_all->GetParError(4)}, {fit_sig_mass_all->GetParameter(3), fit_sig_mass_all->GetParError(3)} };
}


void plot_invmass_AO2D()
{
  //fAnalysisResults = new TFile("AnalysisResults.root");
  //fAnalysisResults = new TFile("AnalysisResults/AnalysisResults-2.root");
  fAnalysisResults = new TFile("AnalysisResults/AnalysisResultsFull.root");
  //fAnalysisResults = new TFile("AnalysisResults-LHC25c3b-qa-MFT-misaligned/AnalysisResultsFull.root");
  //fAnalysisResults = new TFile("AnalysisResults-LHC24am-qa-with-MFT-realignment/AnalysisResultsFull.root");
  //fAnalysisResults = new TFile("AnalysisResults-LHC24am-qa-no-MFT-realignment/AnalysisResultsFull.root");
  //fAnalysisResults = new TFile("AnalysisResults-LHC24am-qa-5/AnalysisResultsFull.root");
  //fAnalysisResults = new TFile("AnalysisResults-LHC24l7/AnalysisResultsFull.root");
  //fAnalysisResults = new TFile("AnalysisResults-LHC22p/AnalysisResultsFull.root");
  //fAnalysisResults = new TFile("AnalysisResults-LHC23zk/AnalysisResultsFull.root");
  //fAnalysisResults = new TFile("AnalysisResults-LHC24am-4/AnalysisResultsFull.root");
  //fAnalysisResults = new TFile("AnalysisResults-LHC24am-7/AnalysisResultsFull.root");
  //fAnalysisResults = new TFile("AnalysisResults-LHC24aq-apass1_muon_matching-qa/AnalysisResultsFull.root");
  //fAnalysisResults = new TFile("AnalysisResults-LHC23h-apass4_skimmed-qa/AnalysisResultsFull.root");

  int rebin = 1;

  gStyle->SetOptStat(0);
  //gStyle->SetOptStat(1111111);
  gStyle->SetOptFit(1111);
  
  TCanvas c("c", "c", 1200, 800);

  TH1* invmassMuonTracks1 = GetTH1(fAnalysisResults, "muon-qa/dimuon/same-event/invariantMass_MuonKine_MuonCuts");
  TH1* invmassMuonTracks1ME = GetTH1(fAnalysisResults, "muon-qa/dimuon/mixed-event/invariantMass_MuonKine_MuonCuts");
  std::cout << std::format("invariantMass_MuonKine_MuonCuts: {} {}", (void*)invmassMuonTracks1, (void*)invmassMuonTracks1ME) << std::endl;
  NormalizeInvmassME(invmassMuonTracks1, invmassMuonTracks1ME);
  invmassMuonTracks1->SetLineColor(kRed);
  invmassMuonTracks1->Rebin(rebin);
  invmassMuonTracks1->GetXaxis()->SetRangeUser(massMin, massMax);
  invmassMuonTracks1->SetTitle(TString::Format("%s (muon kine, muon cuts)", invmassMuonTracks1->GetTitle()));
  invmassMuonTracks1->Draw("E");
  //invmassMuonTracks1ME->Draw("same");
  FitJPsi(invmassMuonTracks1);
  c.SaveAs("invmass_AO2D.pdf(");

  TH1* invmassMuonTracks2 = GetTH1(fAnalysisResults, "muon-qa/dimuon/same-event/invariantMass_MuonKine_GlobalMuonCuts");
  TH1* invmassMuonTracks2ME = GetTH1(fAnalysisResults, "muon-qa/dimuon/mixed-event/invariantMass_MuonKine_GlobalMuonCuts");
  NormalizeInvmassME(invmassMuonTracks2, invmassMuonTracks2ME);
  invmassMuonTracks2->SetLineColor(kBlue - 2);
  invmassMuonTracks2->Rebin(rebin);
  invmassMuonTracks2->GetXaxis()->SetRangeUser(massMin, massMax);
  invmassMuonTracks2->SetTitle(TString::Format("%s (muon kine, global muon cuts)", invmassMuonTracks2->GetTitle()));
  invmassMuonTracks2->Draw("E");
  //invmassMuonTracks2ME->Draw("same");
  FitJPsi(invmassMuonTracks2);
  c.SaveAs("invmass_AO2D.pdf");

  TH1* invmassMuonTracks3 = GetTH1(fAnalysisResults, "muon-qa/dimuon/same-event/invariantMass_MuonKine_GlobalMatchesCuts");
  TH1* invmassMuonTracks3ME = GetTH1(fAnalysisResults, "muon-qa/dimuon/mixed-event/invariantMass_MuonKine_GlobalMatchesCuts");
  NormalizeInvmassME(invmassMuonTracks3, invmassMuonTracks3ME);
  invmassMuonTracks3->SetLineColor(kGreen - 2);
  invmassMuonTracks3->Rebin(rebin);
  invmassMuonTracks3->GetXaxis()->SetRangeUser(massMin, massMax);
  invmassMuonTracks3->SetTitle(TString::Format("%s (muon kine, global matches cuts)", invmassMuonTracks3->GetTitle()));
  invmassMuonTracks3->Draw("E");
  //invmassMuonTracks3ME->Draw("same");
  FitJPsi(invmassMuonTracks3);
  c.SaveAs("invmass_AO2D.pdf");

  TH1* invmassMuonTracks4 = GetTH1(fAnalysisResults, "muon-qa/dimuon/same-event/invariantMass_GlobalMuonKine_GlobalMatchesCuts");
  TH1* invmassMuonTracks4ME = GetTH1(fAnalysisResults, "muon-qa/dimuon/mixed-event/invariantMass_GlobalMuonKine_GlobalMatchesCuts");
  NormalizeInvmassME(invmassMuonTracks4, invmassMuonTracks4ME);
  invmassMuonTracks4->SetLineColor(kMagenta);
  invmassMuonTracks4->Rebin(rebin);
  invmassMuonTracks4->GetXaxis()->SetRangeUser(massMin, massMax);
  invmassMuonTracks4->SetTitle(TString::Format("%s (global muon kine, global matches cuts)", invmassMuonTracks4->GetTitle()));
  invmassMuonTracks4->Draw("E");
  //invmassMuonTracks4ME->Draw("same");
  FitJPsi(invmassMuonTracks4);
  //c.SaveAs("invmass_AO2D.pdf");

  TH1* invmassMuonTracks5 = GetTH1(fAnalysisResults, "muon-qa/dimuon/same-event/invariantMass_ScaledMftKine_GlobalMatchesCuts");
  TH1* invmassMuonTracks5ME = GetTH1(fAnalysisResults, "muon-qa/dimuon/mixed-event/invariantMass_ScaledMftKine_GlobalMatchesCuts");
  NormalizeInvmassME(invmassMuonTracks5, invmassMuonTracks5ME);
  invmassMuonTracks5->SetLineColor(kBlack);
  invmassMuonTracks5->Rebin(rebin);
  invmassMuonTracks5->GetXaxis()->SetRangeUser(massMin, massMax);
  invmassMuonTracks5->SetTitle(TString::Format("%s (scaled MFT kine, global matches cuts)", invmassMuonTracks5->GetTitle()));
  invmassMuonTracks5->Draw("E");
  //invmassMuonTracks4ME->Draw("same");
  FitJPsi(invmassMuonTracks5);
  c.SaveAs("invmass_AO2D.pdf");


  TGraphErrors *gr1 = new TGraphErrors();
  TGraphErrors *gr2 = new TGraphErrors();
  gr2->SetLineStyle(kDashed);
  TMultiGraph* mgr = new TMultiGraph;
  mgr->Add(gr1,"p");
  mgr->Add(gr2,"p");

  TGraphErrors* grv[2]{ gr1, gr2 };

  std::vector<std::string> trackCombinations[2]{ {"TT", "BB"}, {"TB", "BT"} };
  TH2* invmassVsPt[2][2][2];
  for (int i = 0; i < 2; i++) {
    for (int j = 0; j < 2; j++) {
      std::cout << std::format("Loading \"muon-qa/dimuon/same-event/invariantMass_pT_ScaledMftKine_GlobalMatchesCuts_{}\"", trackCombinations[i][j]) << std::endl;
      invmassVsPt[i][j][0] = GetTH2(fAnalysisResults, TString::Format("muon-qa/dimuon/same-event/invariantMass_pT_ScaledMftKine_GlobalMatchesCuts_%s", trackCombinations[i][j].c_str()));
      invmassVsPt[i][j][1] = GetTH2(fAnalysisResults, TString::Format("muon-qa/dimuon/mixed-event/invariantMass_pT_ScaledMftKine_GlobalMatchesCuts_%s", trackCombinations[i][j].c_str()));
    }
  }
  for (int i = 0; i < 2; i++) {
    invmassVsPt[i][0][0]->Add(invmassVsPt[i][1][0]);
    invmassVsPt[i][0][1]->Add(invmassVsPt[i][1][1]);

    invmassVsPt[i][0][0]->Draw("col");
    c.SaveAs("invmass_AO2D.pdf");

    float dpT = 2;
    for (float pt = 0; pt < 20; pt += dpT) {
      int bin1 = invmassVsPt[i][0][0]->GetYaxis()->FindBin(pt + 0.1f);
      int bin2 = invmassVsPt[i][0][0]->GetYaxis()->FindBin(pt + dpT - 0.1f);

      //std::cout << std::format("Bin range: {} -> {}", bin1, bin2) << std::endl;

      TH1* h = invmassVsPt[i][0][0]->ProjectionX(TString::Format("%s_bin_%d_%d", invmassVsPt[i][0][0]->GetName(), i, bin1), bin1, bin2);
      TH1* hME = invmassVsPt[i][0][1]->ProjectionX(TString::Format("%s_bin_%d_%d_ME", invmassVsPt[i][0][1]->GetName(), i, bin1), bin1, bin2);
      NormalizeInvmassME(h, hME);
      h->Rebin(rebin);
      h->GetXaxis()->SetRangeUser(massMin, massMax);
      h->Draw("E");
      auto par = FitJPsi(h);
      std::cout << std::format("J/psi  mass = {} +/- {}  width = {} +/- {}", par.first.first, par.first.second, par.second.first, par.second.second) << std::endl;
      //c.SaveAs("invmass_AO2D.pdf");

      grv[i]->AddPoint(pt + dpT / 2, par.first.first);
      grv[i]->SetPointError(grv[i]->GetN() - 1, dpT / 2, par.first.second);
    }
  }
  c.Clear();
  grv[0]->SetMinimum(2.9);
  grv[1]->SetMinimum(2.9);
  grv[0]->SetMaximum(3.2);
  grv[1]->SetMaximum(3.2);
  mgr->SetMinimum(2.9);
  mgr->SetMaximum(3.2);
  mgr->Draw("A");
  c.SaveAs("invmass_AO2D.pdf");

  invmassMuonTracks5->Draw("hist");
  invmassMuonTracks3->Draw("same");
  c.SaveAs("invmass_AO2D.pdf");

  c.Clear();
  c.Divide(2,2);
  // top-top
  c.cd(1);
  TH1* invmass = GetTH1(fAnalysisResults, "muon-qa/dimuon/same-event/invariantMass_ScaledMftKine_GlobalMatchesCuts_TT");
  TH1* invmassME = GetTH1(fAnalysisResults, "muon-qa/dimuon/mixed-event/invariantMass_ScaledMftKine_GlobalMatchesCuts_TT");
  if (invmass) {
  NormalizeInvmassME(invmass, invmassME);
  invmass->SetLineColor(kBlack);
  invmass->Rebin(rebin);
  invmass->GetXaxis()->SetRangeUser(massMin, massMax);
  invmass->Draw("E");
  //invmassMuonTracks4ME->Draw("same");
  FitJPsi(invmass);
  }
  // top-bottom
  c.cd(2);
  invmass = GetTH1(fAnalysisResults, "muon-qa/dimuon/same-event/invariantMass_ScaledMftKine_GlobalMatchesCuts_TB");
  invmassME = GetTH1(fAnalysisResults, "muon-qa/dimuon/mixed-event/invariantMass_ScaledMftKine_GlobalMatchesCuts_TB");
  if (invmass) {
  NormalizeInvmassME(invmass, invmassME);
  invmass->SetLineColor(kBlack);
  invmass->Rebin(rebin);
  invmass->GetXaxis()->SetRangeUser(massMin, massMax);
  invmass->Draw("E");
  //invmassMuonTracks4ME->Draw("same");
  FitJPsi(invmass);
  }
  // bottom-top
  c.cd(3);
  invmass = GetTH1(fAnalysisResults, "muon-qa/dimuon/same-event/invariantMass_ScaledMftKine_GlobalMatchesCuts_BT");
  invmassME = GetTH1(fAnalysisResults, "muon-qa/dimuon/mixed-event/invariantMass_ScaledMftKine_GlobalMatchesCuts_BT");
  if (invmass) {
  NormalizeInvmassME(invmass, invmassME);
  invmass->SetLineColor(kBlack);
  invmass->Rebin(rebin);
  invmass->GetXaxis()->SetRangeUser(massMin, massMax);
  invmass->Draw("E");
  //invmassMuonTracks4ME->Draw("same");
  FitJPsi(invmass);
  }
  // bottom-bottom
  c.cd(4);
  invmass = GetTH1(fAnalysisResults, "muon-qa/dimuon/same-event/invariantMass_ScaledMftKine_GlobalMatchesCuts_BB");
  invmassME = GetTH1(fAnalysisResults, "muon-qa/dimuon/mixed-event/invariantMass_ScaledMftKine_GlobalMatchesCuts_BB");
  if (invmass) {
  NormalizeInvmassME(invmass, invmassME);
  invmass->SetLineColor(kBlack);
  invmass->Rebin(rebin);
  invmass->GetXaxis()->SetRangeUser(massMin, massMax);
  invmass->Draw("E");
  //invmassMuonTracks4ME->Draw("same");
  FitJPsi(invmass);
  }

  c.SaveAs("invmass_AO2D.pdf");

  c.Clear();

  TH1* histogram;
  TH2* histogram2;

  histogram2 = GetTH2(fAnalysisResults, "muon-qa/dimuon/same-event/invariantMass_MuonKine_vs_GlobalMuonKine");
  histogram2->GetXaxis()->SetRangeUser(massMin, massMax);
  histogram2->GetYaxis()->SetRangeUser(massMin, massMax);
  histogram2->Draw("col");

  c.SaveAs("invmass_AO2D.pdf");

  histogram2 = GetTH2(fAnalysisResults, "muon-qa/dimuon/same-event/invariantMass_ScaledMftKine_vs_GlobalMuonKine");
  histogram2->GetXaxis()->SetRangeUser(massMin, massMax);
  histogram2->GetYaxis()->SetRangeUser(massMin, massMax);
  histogram2->Draw("col");

  c.SaveAs("invmass_AO2D.pdf");

  c.SetLogy(kTRUE);
  histogram = GetTH1(fAnalysisResults, "muon-qa/dimuon/same-event/invariantMassFull_MuonKine_MuonCuts");
  histogram->Draw();
  c.SaveAs("invmass_AO2D.pdf");
  histogram->GetXaxis()->SetRangeUser(1.0, 15.0);
  c.SaveAs("invmass_AO2D.pdf");

  c.SetLogy(kTRUE);
  histogram = GetTH1(fAnalysisResults, "muon-qa/dimuon/same-event/invariantMassFull_ScaledMftKine_GlobalMatchesCuts");
  histogram->Draw();
  c.SaveAs("invmass_AO2D.pdf");
  histogram->GetXaxis()->SetRangeUser(1.0, 15.0);
  c.SaveAs("invmass_AO2D.pdf");

  TH1* histogramLSL = GetTH1(fAnalysisResults, "muon-qa/dimuon/same-event/invariantMass_GlobalMuonKine_GlobalMatchesCuts_leading_subleading");
  histogramLSL->Rebin(4);
  histogramLSL->Draw();

  c.SaveAs("invmass_AO2D.pdf");

  TH1* histogramSLL = GetTH1(fAnalysisResults, "muon-qa/dimuon/same-event/invariantMass_GlobalMuonKine_GlobalMatchesCuts_subleading_leading");
  histogramSLL->Rebin(4);
  histogramSLL->Draw();

  c.SaveAs("invmass_AO2D.pdf");

  histogramLSL->Add(histogramSLL);
  histogramLSL->SetTitle("M_{#mu^{+}#mu^{-}} - leading + sub-leading matches");
  histogramLSL->Draw();

  c.SaveAs("invmass_AO2D.pdf");

  histogram = GetTH1(fAnalysisResults, "muon-qa/dimuon/same-event/invariantMass_GlobalMuonKine_GlobalMatchesCuts_subleading_subleading");
  histogram->SetTitle("M_{#mu^{+}#mu^{-}} - sub-leading matches");
  histogram->Rebin(4);
  histogram->Draw();

  c.SaveAs("invmass_AO2D.pdf");

  histogram2 = GetTH2(fAnalysisResults, "muon-qa/dimuon/same-event/invariantMass_GlobalMuonKine_subleading_vs_leading");
  histogram2->GetXaxis()->SetRangeUser(1, 5);
  histogram2->GetYaxis()->SetRangeUser(1, 5);
  histogram2->Draw("col");

  c.SaveAs("invmass_AO2D.pdf)");
}

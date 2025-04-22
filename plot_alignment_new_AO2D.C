#include <MFTTracking/Constants.h>

TFile* fAnalysisResults;

constexpr int nPoints = 5;

TH1* GetTH1(TFile* f, TString histname)
{
  return (TH1*)f->Get(histname);

  //TString histname = TString::Format("ST%d/DE%d/Occupancy_B_XY_%d", station, de, de);
  TKey *key = f->GetKey(histname);
  std::cout << "histname: " << histname << "  key: " <<key << std::endl;
  if (!key) return NULL;
  return (TH1*)key->ReadObjectAny(TH1::Class());
}


TH2* GetTH2(TFile* f, TString histname)
{
  return (TH2*)f->Get(histname);

  //TString histname = TString::Format("ST%d/DE%d/Occupancy_B_XY_%d", station, de, de);
  TKey *key = f->GetKey(histname);
  std::cout << "histname: " << histname << "  key: " <<key << std::endl;
  if (!key) return NULL;
  return (TH2*)key->ReadObjectAny(TH2::Class());
}


void PlotDXYProjection(const char* fullHistName, const char* fullHistNameME, TH2* histogram2, TH2* histogram2ME, double scaleME, float yMin, float yMax, int projRebin, TCanvas& c, bool printFits = false)
{
  c.Clear();

  TH1* histogramMean = histogram2->ProjectionX();
  histogramMean->SetName(TString::Format("%s-mean", fullHistName));
  histogramMean->SetTitle(TString::Format("%s (mean)", histogram2->GetTitle()));
  histogramMean->SetMinimum(yMin);
  histogramMean->SetMaximum(yMax);
  TH1* histogramSigma = histogram2->ProjectionX();
  histogramSigma->SetName(TString::Format("%s-sigma", fullHistName));
  histogramSigma->SetTitle(TString::Format("%s (sigma)", histogram2->GetTitle()));
  histogramSigma->SetMinimum(0);
  histogramSigma->SetMaximum(yMax);

  int entriesMin = histogram2->GetEntries() / 20;
  for (int bin = 1; bin <= histogram2->GetXaxis()->GetNbins(); bin++) {
    TH1* proj = histogram2->ProjectionY(TString::Format("%s-%d", fullHistName, bin), bin, bin);
    TH1* projME = histogram2ME->ProjectionY(TString::Format("%s-%d", fullHistNameME, bin), bin, bin);

    proj->Rebin(projRebin);
    projME->Rebin(projRebin);

    projME->Scale(scaleME);
    proj->SetTitle(TString::Format("%s - bin %d", histogram2->GetTitle(), bin));
    double mean = -10000;
    double meanErr = 0;
    double sigma = -10000;
    double sigmaErr = 0;
    if (proj->GetEntries() > entriesMin) {

      if (printFits) {
        proj->SetLineColor(kRed);
        proj->Draw("E");
        projME->Draw("L same");
        c.SaveAs("alignment_AO2D.pdf");
      }

      proj->Add(projME, -1);

      //proj->Rebin(projRebin);
      int valuePeak = proj->GetMaximum();
      int binPeak = proj->GetMaximumBin();
      double xPeak = proj->GetXaxis()->GetBinCenter(binPeak);
      if (printFits) {
        std::cout << std::format("[\"{}\"] bin {} max={}  binPeak={}  xPeak={}\n", fullHistName, bin, valuePeak, binPeak, xPeak);
        //std::cout << std::format("[\"{}\"] bin {}  max={:0.2f}  binPeak={}  xPeak={:0.2f}",
        //    fullHistName, bin, valuePeak, binPeak, xPeak) << std::endl;
      }
      //std::cout << std::format("Bin {}  max={}  binPeak={}  xPeak={}\n", bin, valuePeak, binPeak, xPeak);
      TF1 fgaus("fgaus", "gaus", xPeak - 1, xPeak + 1);
      //fgaus.SetParLimits(0, 0, 100000000000.0);
      fgaus.FixParameter(1, xPeak);
      fgaus.SetParameter(2, 1);
      fgaus.SetParLimits(2, 0, 100);
      //proj->Fit("fgaus", "RQBN");
      fgaus.ReleaseParameter(1);
      //proj->Fit("fgaus", "RQBN");

      TF1 fgaus2("fgaus2", "gaus(0)+pol0(3)");
      fgaus2.SetNpx(1000);
      fgaus2.SetLineColor(kBlack);
      fgaus2.SetParameter(0, valuePeak / 10);
      fgaus2.SetParLimits(0, 0, valuePeak * 10);
      fgaus2.SetParameter(1, xPeak);
      fgaus2.FixParameter(1, xPeak);
      fgaus2.SetParameter(2, 1);
      fgaus2.SetParLimits(2, 0, 10);
      fgaus2.SetParameter(3, 0);
      fgaus2.SetParameter(4, 0);
      proj->Fit("fgaus2", "BQN");
      fgaus2.ReleaseParameter(1);
      fgaus2.SetParLimits(1, -10, 10);
      proj->Fit("fgaus2", "BQN");

      double m = fgaus2.GetParameter(1);
      double s2 = fgaus2.GetParameter(2) * 2;
      TF1 fgaus3("fgaus3", "gaus(0)", m - s2, m + s2);
      fgaus3.SetNpx(1000);
      fgaus3.SetLineColor(kBlack);
      fgaus3.SetParameters(
          fgaus2.GetParameter(0),
          fgaus2.GetParameter(1),
          fgaus2.GetParameter(2)
      );
      proj->Fit("fgaus3", "RBQ");

      if (printFits) {
        proj->Draw("E");
        c.SaveAs("alignment_AO2D.pdf");
      }

      mean = fgaus3.GetParameter(1);
      meanErr = fgaus3.GetParError(1);
      sigma = fgaus2.GetParameter(2);
      sigmaErr = fgaus2.GetParError(2);
    }

    histogramMean->SetBinContent(bin, mean);
    histogramMean->SetBinError(bin, meanErr);

    histogramSigma->SetBinContent(bin, sigma);
    histogramSigma->SetBinError(bin, sigmaErr);
  }
  histogramMean->Draw("E");
  c.SaveAs("alignment_AO2D.pdf");
  histogramSigma->Draw("E");
  c.SaveAs("alignment_AO2D.pdf");
}

std::pair<double, double> PlotDXY(std::string histName, TCanvas& c, bool printFits = false)
{
  c.Clear();
  c.Divide(2, 2);

  std::string fullHistName = std::string("qa-muon/alignment/") + histName;
  TH2* histogram2 = GetTH2(fAnalysisResults, fullHistName);
  std::cout << fullHistName << " -> " << histogram2 << std::endl;

  c.cd(1);
  histogram2->Draw("col");

  std::string fullHistNameME = std::string("qa-muon/alignment/mixed-events/") + histName;
  TH2* histogram2ME = GetTH2(fAnalysisResults, fullHistNameME);

  c.cd(2);
  histogram2ME->Draw("col");

  TH1* proj = histogram2->ProjectionY((fullHistName + "_py").c_str());
  proj->Rebin(2);
  TH1* projME = histogram2ME->ProjectionY((fullHistNameME + "_py").c_str());
  projME->Rebin(2);

  double integral1 = proj->Integral(1, proj->GetXaxis()->FindBin(-5));
  double integral2 = proj->Integral(proj->GetXaxis()->FindBin(5), proj->GetXaxis()->GetNbins());
  double integralME1 = projME->Integral(1, projME->GetXaxis()->FindBin(-5));
  double integralME2 = projME->Integral(projME->GetXaxis()->FindBin(5), projME->GetXaxis()->GetNbins());
  double scaleME = (integral1 + integral2) / (integralME1 + integralME2);

  c.cd(3);
  projME->Scale(scaleME);
  proj->SetLineColor(kRed);
  proj->Draw("E");
  projME->Draw("E same");

  c.cd(4);
  TH1* proj2 = (TH1*)proj->Clone();
  proj2->Add(projME, -1);
  proj2->SetLineColor(kRed);
  int valuePeak = proj2->GetMaximum();
  int binPeak = proj2->GetMaximumBin();
  double xPeak = proj2->GetXaxis()->GetBinCenter(binPeak);
  TF1 fgaus("fgaus", "gausn(0)+pol0(3)");
  fgaus.SetNpx(1000);
  fgaus.SetLineColor(kBlack);
  fgaus.SetParameter(0, valuePeak / 10);
  fgaus.SetParLimits(0, 0, valuePeak * 10);
  fgaus.SetParameter(1, xPeak);
  fgaus.SetParameter(2, 1);
  fgaus.SetParLimits(2, 0, 10);
  fgaus.SetParameter(3, 0);
  fgaus.SetParameter(4, 0);
  proj2->Fit("fgaus", "BQN");
  proj2->Draw("E");

  double m = fgaus.GetParameter(1);
  double s2 = fgaus.GetParameter(2) * 2;
  TF1 fgaus2("fgaus2", "gaus(0)", m - s2, m + s2);
  fgaus2.SetNpx(1000);
  fgaus2.SetLineColor(kBlack);
  fgaus2.SetParameters(
      fgaus.GetParameter(0),
      fgaus.GetParameter(1),
      fgaus.GetParameter(2)
  );
  proj2->Fit("fgaus2", "BQR");
  proj2->GetXaxis()->SetRangeUser(-10.0, 10.0);
  proj2->Draw("E");

  c.SaveAs("alignment_AO2D.pdf");

  PlotDXYProjection(fullHistName.c_str(), fullHistNameME.c_str(), histogram2, histogram2ME, scaleME, -5.0, 5.0, 4, c, printFits);

  return {fgaus2.GetParameter(1), fgaus2.GetParError(1)};
}

std::pair<double, double> PlotDTheta(std::string histName, TCanvas& c, bool printFits = false)
{
  c.Clear();
  c.Divide(2, 2);

  std::string fullHistName = std::string("qa-muon/alignment/") + histName;
  TH2* histogram2 = GetTH2(fAnalysisResults, fullHistName);
  //std::cout << fullHistName << " -> " << histogram2 << std::endl;

  c.cd(1);
  histogram2->Draw("col");

  std::string fullHistNameME = std::string("qa-muon/alignment/mixed-events/") + histName;
  TH2* histogram2ME = GetTH2(fAnalysisResults, fullHistNameME);

  c.cd(2);
  histogram2ME->Draw("col");

  TH1* proj = histogram2->ProjectionY((fullHistName + "_py").c_str());
  TH1* projME = histogram2ME->ProjectionY((fullHistNameME + "_py").c_str());

  double integral1 = proj->Integral(1, proj->GetXaxis()->FindBin(-2));
  double integral2 = proj->Integral(proj->GetXaxis()->FindBin(2), proj->GetXaxis()->GetNbins());
  double integralME1 = projME->Integral(1, projME->GetXaxis()->FindBin(-2));
  double integralME2 = projME->Integral(projME->GetXaxis()->FindBin(2), projME->GetXaxis()->GetNbins());
  double scaleME = (integral1 + integral2) / (integralME1 + integralME2);

  c.cd(3);
  projME->Scale(scaleME);
  proj->SetLineColor(kRed);
  proj->Draw("E");
  projME->Draw("E same");

  c.cd(4);
  TH1* proj2 = (TH1*)proj->Clone();
  proj2->Add(projME, -1);
  proj2->SetLineColor(kRed);
  int valuePeak = proj2->GetMaximum();
  int binPeak = proj2->GetMaximumBin();
  double xPeak = proj2->GetXaxis()->GetBinCenter(binPeak);
  TF1 fgaus("fgaus", "gausn(0)+pol0(3)");
  fgaus.SetNpx(1000);
  fgaus.SetLineColor(kBlack);
  fgaus.SetParameter(1, xPeak);
  fgaus.SetParameter(2, 0.5);
  fgaus.SetParLimits(2, 0, 100);
  proj2->Fit("fgaus", "BQ");
  proj2->Draw("E");

  c.SaveAs("alignment_AO2D.pdf");

  PlotDXYProjection(fullHistName.c_str(), fullHistNameME.c_str(), histogram2, histogram2ME, scaleME, -0.5, 0.5, 2, c, printFits);

  return {fgaus.GetParameter(1), fgaus.GetParError(1)};
}

Double_t DoubleSidedCB2(double x, double mu, double width, double a1, double p1, double a2, double p2)
{
  double u   = (x-mu)/width;
  double A1  = TMath::Power(p1/TMath::Abs(a1),p1)*TMath::Exp(-a1*a1/2);
  double A2  = TMath::Power(p2/TMath::Abs(a2),p2)*TMath::Exp(-a2*a2/2);
  double B1  = p1/TMath::Abs(a1) - TMath::Abs(a1);
  double B2  = p2/TMath::Abs(a2) - TMath::Abs(a2);

  double result(1);
  if      (u<-a1) result *= A1*TMath::Power(B1-u,-p1);
  else if (u<a2)  result *= TMath::Exp(-u*u/2);
  else            result *= A2*TMath::Power(B2+u,-p2);
  return result;
}


double DoubleSidedCB(double* x, double *par)
{
  return(par[0] * DoubleSidedCB2(x[0], par[1],par[2],par[3],par[4],par[5],par[6]));
}

std::pair<double, double> PlotDCAMFT(std::string histName)
{
  std::string fullHistName = std::string("qa-muon/alignment/") + histName;
  TH1* histogram = GetTH1(fAnalysisResults, fullHistName);
  //std::cout << fullHistName << " -> " << histogram << std::endl;

  std::string fullHistNameME = std::string("qa-muon/alignment/mixed-events/") + histName;
  TH1* histogramME = GetTH1(fAnalysisResults, fullHistNameME);

  double integral1 = histogram->Integral(1, histogram->GetXaxis()->FindBin(-0.2));
  double integral2 = histogram->Integral(histogram->GetXaxis()->FindBin(0.2), histogram->GetXaxis()->GetNbins());
  double integralME1 = histogramME->Integral(1, histogramME->GetXaxis()->FindBin(-0.2));
  double integralME2 = histogramME->Integral(histogramME->GetXaxis()->FindBin(0.2), histogramME->GetXaxis()->GetNbins());
  double scaleME = (integral1 + integral2) / (integralME1 + integralME2);

  histogramME->Scale(scaleME);

  histogram->SetLineColor(kRed);
  histogram->GetXaxis()->SetRangeUser(-0.2, 0.2);

  int valuePeak = histogram->GetMaximum();
  int binPeak = histogram->GetMaximumBin();
  double xPeak = histogram->GetXaxis()->GetBinCenter(binPeak);
  TF1 fgaus("fgaus", "gaus(0)", -0.05, 0.05);
  fgaus.SetNpx(1000);
  fgaus.SetLineColor(kBlack);
  fgaus.SetParameter(0, valuePeak / 10);
  fgaus.SetParLimits(0, 0, valuePeak * 10);
  fgaus.SetParameter(1, xPeak);
  fgaus.SetParameter(2, 1);
  fgaus.SetParLimits(2, 0, 10);
  //fgaus.SetParameter(4, xPeak);
  //fgaus.SetParameter(5, 10);
  //histogram->Fit("fgaus", "BQR");
  //histogram->Draw("E");

  TF1 fcb("fcb", DoubleSidedCB, -0.2, 0.2, 7);
  fcb.SetNpx(1000);
  fcb.SetLineColor(kBlack);
  //fcb.SetParameter(0, valuePeak);
  double par[7];
  par[0]=35000;
  par[1]=xPeak;
  par[2]=0.01;
  par[3]=1;
  par[4]=1;
  par[5]=1;
  par[6]=1;
  fcb.SetParameters(&par[0]);

  fcb.FixParameter(1, xPeak);
  fcb.SetParLimits(2, 0.01, 0.1);
  histogram->Fit("fcb", "BRN");
  fcb.ReleaseParameter(1);
  histogram->Fit("fcb", "BR");

  histogram->Draw("E");
  //histogramME->Draw("E same");
  //fpow.Draw("same");

  return {fgaus.GetParameter(1), fgaus.GetParError(1)};
  //return {histogram->GetMean(), fgaus.GetParError(1)};
}

std::pair<double, double> PlotDCAMCH(std::string histName)
{
  std::string fullHistName = std::string("qa-muon/alignment/") + histName;
  TH1* histogram = GetTH1(fAnalysisResults, fullHistName);
  //std::cout << fullHistName << " -> " << histogram << std::endl;
  if (!histogram)
    return {};
  histogram->Rebin(4);

  histogram->Draw("E");

  int valuePeak = histogram->GetMaximum();
  int binPeak = histogram->GetMaximumBin();
  double xPeak = histogram->GetXaxis()->GetBinCenter(binPeak);
  TF1 fgaus("fgaus", "gausn(0)", -4, 4);
  fgaus.SetNpx(1000);
  fgaus.SetLineColor(kBlack);
  fgaus.SetParameter(0, valuePeak / 10);
  fgaus.SetParLimits(0, 0, valuePeak * 10);
  fgaus.SetParameter(1, xPeak);
  fgaus.SetParameter(2, 1);
  fgaus.SetParLimits(2, 0, 10);
  //histogram->Fit("fgaus", "RBQ");

  TF1 fcb("fcb", DoubleSidedCB, -6, 6, 7);
  fcb.SetNpx(1000);
  fcb.SetLineColor(kBlack);
  //fcb.SetParameter(0, valuePeak);
  double par[7];
  par[0]=35000;
  par[1]=xPeak;
  par[2]=0.5;
  par[3]=1;
  par[4]=1;
  par[5]=1;
  par[6]=1;
  fcb.SetParameters(&par[0]);

  fcb.FixParameter(1, xPeak);
  fcb.SetParLimits(2, 0.5, 10);
  histogram->Fit("fcb", "BRN");
  fcb.ReleaseParameter(1);
  histogram->Fit("fcb", "BRN");
  histogram->Fit("fcb", "BR");

  histogram->Draw("E");

  //return {histogram->GetMean(), fgaus.GetParError(1)};
  //return {fgaus.GetParameter(1), fgaus.GetParError(1)};
  return {fcb.GetParameter(1), fcb.GetParError(1)};
}

void PlotZTrend(int n, double* xv, std::array<std::array<std::pair<double, double>, 4>, nPoints + 1>& values, const char* title, double ymin, double ymax, TCanvas& c)
{
  double exv[nPoints + 1] = {0, 0, 0, 0, 0};
  double yv[nPoints + 1];
  double eyv[nPoints + 1];
  std::array<std::string, 4> quadrants = {"Q0", "Q1", "Q2", "Q3"};

  int colors[4] = {kBlue, kRed, kOrange, kCyan};
  int markers[4] = {kStar, kCircle, kMultiply, kFullDotLarge};

  c.Clear();

  TMultiGraph* mg = new TMultiGraph();
  mg->SetTitle(title);
  TLegend* legend = new TLegend(0.6, 0.8, 0.9, 0.9);
  legend->SetNColumns(4);
  for (int j = 0; j < quadrants.size(); j++) {
    for (int i = 0; i < n; i++) {
      yv[i] = values[i][j].first;
      eyv[i] = values[i][j].second;
    }
    TGraphErrors* gr = new TGraphErrors(n, xv, yv, exv, eyv);
    gr->SetLineColor(colors[j]);
    gr->SetMarkerColor(colors[j]);
    gr->SetMarkerStyle(markers[1]);
    gr->SetMarkerSize(2);
    mg->Add(gr,"lp");
    auto* entry = legend->AddEntry(gr, (quadrants[j] + " ").c_str(), "P");
    entry->SetTextColor(colors[j]);
  }
  mg->Draw("a");
  mg->SetMinimum(ymin);
  mg->SetMaximum(ymax);
  legend->Draw();
  c.SaveAs("alignment_AO2D.pdf");
}

void plot_alignment_new_AO2D()
{
  //fAnalysisResults = new TFile("AnalysisResults.root");
  fAnalysisResults = new TFile("AnalysisResults/AnalysisResultsFull.root");
  //fAnalysisResults = new TFile("AnalysisResults-LHC24l7/AnalysisResultsFull.root");
  //fAnalysisResults = new TFile("AnalysisResults-LHC22p/AnalysisResultsFull.root");
  //fAnalysisResults = new TFile("AnalysisResults-LHC23zk/AnalysisResultsFull.root");
  //fAnalysisResults = new TFile("AnalysisResults-LHC24am-4/AnalysisResultsFull.root");
  //fAnalysisResults = new TFile("AnalysisResults-LHC24am-7/AnalysisResultsFull.root");
  //fAnalysisResults = new TFile("AnalysisResults-LHC23h-apass4_skimmed-qa/AnalysisResultsFull.root");

  constexpr double firstMFTPlaneZ = o2::mft::constants::mft::LayerZCoordinate()[0];
  constexpr double lastMFTPlaneZ = o2::mft::constants::mft::LayerZCoordinate()[9];
  std::array<double, nPoints> zRefPlane{
      firstMFTPlaneZ,
      lastMFTPlaneZ,
      -90.0,
      -300.0,
      //-505.0,
      -520.0
  };

  std::vector<std::string> referencePlaneNames = {
      "MFT-begin",
      "MFT-end",
      "absorber-begin",
      "absorber-mid",
      //"absorber-end",
      "MCH-begin"
  };

  std::array<std::string, 4> quadrants = {"Q0", "Q1", "Q2", "Q3"};


  std::array<std::pair<double, double>, 4> DCAx;
  std::array<std::pair<double, double>, 4> DCAy;

  std::array<std::array<std::pair<double, double>, 4>, nPoints + 1> meanDx;
  std::array<std::array<std::pair<double, double>, 4>, nPoints + 1> meanDy;
  std::array<std::array<std::pair<double, double>, 4>, nPoints + 1> meanDThetax;
  std::array<std::array<std::pair<double, double>, 4>, nPoints + 1> meanDThetay;

  //gStyle->SetOptStat(0);
  //gStyle->SetOptStat(1111);
  gStyle->SetOptFit(1111);
  
  TCanvas c("c", "c", 1200, 800);
  c.SaveAs("alignment_AO2D.pdf(");

  c.Clear();
  c.Divide(2, 2);
  for (int j = 0; j < quadrants.size(); j++) {
    if (j == 0) c.cd(2);
    if (j == 1) c.cd(1);
    if (j == 2) c.cd(3);
    if (j == 3) c.cd(4);
    meanDx[0][j] = PlotDCAMFT(std::string("DCA/MFT/") + quadrants[j] + "/DCA_x");
  }
  c.SaveAs("alignment_AO2D.pdf");

  c.Clear();
  c.Divide(2, 2);
  for (int j = 0; j < quadrants.size(); j++) {
    if (j == 0) c.cd(2);
    if (j == 1) c.cd(1);
    if (j == 2) c.cd(3);
    if (j == 3) c.cd(4);
    meanDx[0][j] = PlotDCAMFT(std::string("DCA/MFT/") + quadrants[j] + "/DCA_y");
  }
  c.SaveAs("alignment_AO2D.pdf");

  c.Clear();
  c.Divide(2, 2);
  for (int j = 0; j < quadrants.size(); j++) {
    if (j == 0) c.cd(2);
    if (j == 1) c.cd(1);
    if (j == 2) c.cd(3);
    if (j == 3) c.cd(4);
    meanDx[0][j] = PlotDCAMCH(std::string("DCA/MCH/") + quadrants[j] + "/DCA_x");
  }
  c.SaveAs("alignment_AO2D.pdf");
  c.Clear();
  c.Divide(2, 2);
  for (int j = 0; j < quadrants.size(); j++) {
    if (j == 0) c.cd(2);
    if (j == 1) c.cd(1);
    if (j == 2) c.cd(3);
    if (j == 3) c.cd(4);
    PlotDCAMCH(std::string("DCA/MCH/") + quadrants[j] + "/DCA_x_pos");
  }
  c.SaveAs("alignment_AO2D.pdf");
  c.Clear();
  c.Divide(2, 2);
  for (int j = 0; j < quadrants.size(); j++) {
    if (j == 0) c.cd(2);
    if (j == 1) c.cd(1);
    if (j == 2) c.cd(3);
    if (j == 3) c.cd(4);
    PlotDCAMCH(std::string("DCA/MCH/") + quadrants[j] + "/DCA_x_neg");
  }
  c.SaveAs("alignment_AO2D.pdf");

  c.Clear();
  c.Divide(2, 2);
  for (int j = 0; j < quadrants.size(); j++) {
    if (j == 0) c.cd(2);
    if (j == 1) c.cd(1);
    if (j == 2) c.cd(3);
    if (j == 3) c.cd(4);
    meanDy[0][j] = PlotDCAMCH(std::string("DCA/MCH/") + quadrants[j] + "/DCA_y");
  }
  c.SaveAs("alignment_AO2D.pdf");
  c.Clear();
  c.Divide(2, 2);
  for (int j = 0; j < quadrants.size(); j++) {
    if (j == 0) c.cd(2);
    if (j == 1) c.cd(1);
    if (j == 2) c.cd(3);
    if (j == 3) c.cd(4);
    PlotDCAMCH(std::string("DCA/MCH/") + quadrants[j] + "/DCA_y_pos");
  }
  c.SaveAs("alignment_AO2D.pdf");
  c.Clear();
  c.Divide(2, 2);
  for (int j = 0; j < quadrants.size(); j++) {
    if (j == 0) c.cd(2);
    if (j == 1) c.cd(1);
    if (j == 2) c.cd(3);
    if (j == 3) c.cd(4);
    PlotDCAMCH(std::string("DCA/MCH/") + quadrants[j] + "/DCA_y_neg");
  }
  c.SaveAs("alignment_AO2D.pdf");

  for (int j = 0; j < quadrants.size(); j++) {
    for (int i = 0; i < referencePlaneNames.size(); i++) {
      meanDx[i + 1][j] = PlotDXY(referencePlaneNames[i] + "/" + quadrants[j] + "/dx_vs_y", c, false);
    }
  }

  for (int j = 0; j < quadrants.size(); j++) {
    for (int i = 0; i < referencePlaneNames.size(); i++) {
      PlotDXY(referencePlaneNames[i] + "/" + quadrants[j] + "/dx_vs_y", c, false);
    }
  }

  for (int j = 0; j < quadrants.size(); j++) {
    for (int i = 0; i < referencePlaneNames.size(); i++) {
      meanDy[i + 1][j] = PlotDXY(referencePlaneNames[i] + "/" + quadrants[j] + "/dy_vs_x", c, false);
    }
  }

  for (int j = 0; j < quadrants.size(); j++) {
    for (int i = 0; i < referencePlaneNames.size(); i++) {
      PlotDXY(referencePlaneNames[i] + "/" + quadrants[j] + "/dy_vs_x", c, false);
    }
  }

  for (int j = 0; j < quadrants.size(); j++) {
    for (int i = 0; i < referencePlaneNames.size(); i++) {
      meanDThetax[i][j] = PlotDTheta(referencePlaneNames[i] + "/" + quadrants[j] + "/dthetax_vs_x", c, false);
    }
  }

  for (int j = 0; j < quadrants.size(); j++) {
    for (int i = 0; i < referencePlaneNames.size(); i++) {
      PlotDTheta(referencePlaneNames[i] + "/" + quadrants[j] + "/dthetax_vs_y", c, false);
    }
  }

  for (int j = 0; j < quadrants.size(); j++) {
    for (int i = 0; i < referencePlaneNames.size(); i++) {
      PlotDTheta(referencePlaneNames[i] + "/" + quadrants[j] + "/dthetax_vs_thetax", c, false);
    }
  }

  for (int j = 0; j < quadrants.size(); j++) {
    for (int i = 0; i < referencePlaneNames.size(); i++) {
      PlotDTheta(referencePlaneNames[i] + "/" + quadrants[j] + "/dthetay_vs_x", c, false);
    }
  }

  for (int j = 0; j < quadrants.size(); j++) {
    for (int i = 0; i < referencePlaneNames.size(); i++) {
      meanDThetay[i][j] = PlotDTheta(referencePlaneNames[i] + "/" + quadrants[j] + "/dthetay_vs_y", c, false);
    }
  }

  for (int j = 0; j < quadrants.size(); j++) {
    for (int i = 0; i < referencePlaneNames.size(); i++) {
      PlotDTheta(referencePlaneNames[i] + "/" + quadrants[j] + "/dthetay_vs_thetay", c, false);
    }
  }


  double xv[nPoints + 1] = {0, -zRefPlane[0], -zRefPlane[1], -zRefPlane[2], -zRefPlane[3], -zRefPlane[4]};
  double xv2[nPoints] = {-zRefPlane[0], -zRefPlane[1], -zRefPlane[2], -zRefPlane[3], -zRefPlane[4]};

  PlotZTrend(nPoints + 1, xv, meanDx, "#Delta(x) vs. z;z (cm); #Delta(x) (cm)", -1.0, 1.0, c);
  PlotZTrend(nPoints, xv2, meanDThetax, "#Delta(#theta_{x}) vs. z;z (cm); #Delta(#theta_{x}) (cm)", -0.2, 0.2, c);

  PlotZTrend(nPoints + 1, xv, meanDy, "#Delta(y) vs. z;z (cm); #Delta(y) (cm)", -2.0, 2.0, c);
  PlotZTrend(nPoints, xv2, meanDThetay, "#Delta(#theta_{y}) vs. z;z (cm); #Delta(#theta_{y}) (cm)", -0.2, 0.2, c);

  c.Clear();
  c.SaveAs("alignment_AO2D.pdf)");
}

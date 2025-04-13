#include <MFTTracking/Constants.h>

TFile* fAnalysisResults;

TH1* histDxVsY;

double scaleME = 1;

constexpr int nPoints = 4;

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
  // skip first bin because it is biased?
  for (int bin = 2; bin <= histogram2->GetXaxis()->GetNbins(); bin++) {
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
        c.SaveAs("residuals_AO2D.pdf");
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

      TF1 fcb("fcb","[0]*ROOT::Math::crystalball_function(x, [1], [2], [3], [4]) + [5]");
      fcb.SetParameters(100, 0.6, -2.13903e+06, 1, xPeak);
      fcb.SetParLimits(3, 0, 100);
      fcb.SetNpx(1000);
      fcb.SetLineColor(kBlack);
      proj->Fit("fcb", "BNQ");
      proj->Fit("fcb", "BQ");

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
      //proj->Fit("fgaus2", "BQN");
      fgaus2.ReleaseParameter(1);
      fgaus2.SetParLimits(1, -10, 10);
      //proj->Fit("fgaus2", "BQ");

      if (printFits) {
        proj->Draw("E");
        c.SaveAs("residuals_AO2D.pdf");
      }

      mean = fcb.GetParameter(4);
      meanErr = fcb.GetParError(4);
      sigma = fcb.GetParameter(3);
      sigmaErr = fcb.GetParError(3);

      //mean = fgaus2.GetParameter(1);
      //meanErr = fgaus2.GetParError(1);
      //sigma = fgaus2.GetParameter(2);
      //sigmaErr = fgaus2.GetParError(2);
    }

    histogramMean->SetBinContent(bin, mean);
    histogramMean->SetBinError(bin, meanErr);

    histogramSigma->SetBinContent(bin, sigma);
    histogramSigma->SetBinError(bin, sigmaErr);
  }
  histogramMean->Draw("E");
  c.SaveAs("residuals_AO2D.pdf");
  histogramSigma->Draw("E");
  c.SaveAs("residuals_AO2D.pdf");

  histDxVsY = histogramMean;
}

std::pair<double, double> PlotDXY(std::string histName, TCanvas& c, bool printFits = false)
{
  c.Clear();
  c.Divide(2, 2);

  std::string fullHistName = std::string("qa-muon/alignment/residuals/") + histName;
  TH2* histogram2 = GetTH2(fAnalysisResults, fullHistName);
  std::cout << fullHistName << " -> " << histogram2 << std::endl;

  c.cd(1);
  histogram2->Draw("col");

  std::string fullHistNameME = std::string("qa-muon/alignment/mixed-events/residuals/") + histName;
  TH2* histogram2ME = GetTH2(fAnalysisResults, fullHistNameME);

  c.cd(2);
  histogram2ME->Draw("col");

  TH1* proj = histogram2->ProjectionY((fullHistName + "_py").c_str());
  TH1* projME = histogram2ME->ProjectionY((fullHistNameME + "_py").c_str());

  double integral1 = proj->Integral(1, proj->GetXaxis()->FindBin(-15));
  double integral2 = proj->Integral(proj->GetXaxis()->FindBin(15), proj->GetXaxis()->GetNbins());
  double integralME1 = projME->Integral(1, projME->GetXaxis()->FindBin(-15));
  double integralME2 = projME->Integral(projME->GetXaxis()->FindBin(15), projME->GetXaxis()->GetNbins());
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
  proj2->Fit("fgaus", "BQ");


  TF1 fcb("fcb","[0]*ROOT::Math::crystalball_function(x, [1], [2], [3], [4]) + [5]");
  fcb.SetParameters(100, 0.6, -2.13903e+06, 1, xPeak);
  fcb.SetParLimits(3, 0, 100);
  fcb.SetNpx(1000);
  fcb.SetLineColor(kBlack);
  //proj2->Fit("fcb", "BNQ");
  //proj2->Fit("fcb", "BQ");

  proj2->Draw("E");

  c.SaveAs("residuals_AO2D.pdf");

  PlotDXYProjection(fullHistName.c_str(), fullHistNameME.c_str(), histogram2, histogram2ME, scaleME, -5.0, 5.0, 4, c, printFits);

  return {fgaus.GetParameter(1), fgaus.GetParError(1)};
  //return {fcb.GetParameter(4), fgaus.GetParError(4)};
}

void PlotZTrend(int n, double* xv, std::array<std::array<std::pair<double, double>, 10>, 4>& values, const char* title, double ymin, double ymax, TCanvas& c)
{
  double exv[10] = {0, 0, 0, 0, 0, 0, 0, 0, 0, 0};
  double yv[10];
  double eyv[10];
  std::array<std::string, 4> quadrants = {"Q0", "Q1", "Q2", "Q3"};

  int colors[4] = {kBlue, kRed, kGreen - 2, kMagenta};
  int markers[4] = {kStar, kCircle, kMultiply, kFullDotLarge};

  c.Clear();

  TMultiGraph* mg = new TMultiGraph();
  mg->SetTitle(title);
  TLegend* legend = new TLegend(0.6, 0.8, 0.9, 0.9);
  legend->SetNColumns(4);
  for (int j = 0; j < quadrants.size(); j++) {
    for (int i = 0; i < 10; i++) {
      yv[i] = values[j][i].first;
      eyv[i] = values[j][i].second;
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
  c.SaveAs("residuals_AO2D.pdf");
}

void plot_residuals_AO2D()
{
  //fAnalysisResults = new TFile("AnalysisResults.root");
  fAnalysisResults = new TFile("AnalysisResults/AnalysisResultsFull.root");
  //fAnalysisResults = new TFile("AnalysisResults-LHC24l7/AnalysisResultsFull.root");
  //fAnalysisResults = new TFile("AnalysisResults-LHC22p/AnalysisResultsFull.root");
  //fAnalysisResults = new TFile("AnalysisResults-LHC23zk/AnalysisResultsFull.root");
  //fAnalysisResults = new TFile("AnalysisResults-LHC24am-4/AnalysisResultsFull.root");
  //fAnalysisResults = new TFile("AnalysisResults-LHC24am-7/AnalysisResultsFull.root");
  //fAnalysisResults = new TFile("AnalysisResults-LHC23h-apass4_skimmed-qa/AnalysisResultsFull.root");

  std::array<std::string, 4> quadrants = {"Q0", "Q1", "Q2", "Q3"};


  std::array<std::pair<double, double>, 4> DCAx;
  std::array<std::pair<double, double>, 4> DCAy;

  std::array<std::array<std::pair<double, double>, 10>, 4> meanDx;
  std::array<std::array<std::pair<double, double>, 10>, 4> meanDy;

  std::array<std::array<TH1*, 4>, 10> dxVsYhistograms;

  //gStyle->SetOptStat(0);
  //gStyle->SetOptStat(1111);
  gStyle->SetOptFit(1111);
  
  TCanvas c("c", "c", 1200, 800);
  c.SaveAs("residuals_AO2D.pdf(");

  c.Divide(2, 2);

  for (int i = 0; i < 10; i++) {
    for (int j = 0; j < quadrants.size(); j++) {
      bool print = false;
      if (i == 9 && j == 0) print = true;
      meanDx[j][i] = PlotDXY(quadrants[j] + "/CH" + std::to_string(i+1) + "/dx_vs_y", c, print);
      dxVsYhistograms[i][j] = histDxVsY;
      meanDy[j][i] = PlotDXY(quadrants[j] + "/CH" + std::to_string(i+1) + "/dy_vs_x", c, false);
    }
  }

  c.Clear();
  c.Divide(2, 2);

  for (int i = 0; i < 10; i++) {
    for (int j = 0; j < 4; j++) {
      if (j == 0) c.cd(2);
      if (j == 1) c.cd(1);
      if (j == 2) c.cd(3);
      if (j == 3) c.cd(4);

      dxVsYhistograms[i][j]->SetTitle(std::format("#Deltax vs. y - CH{} Q{}", i + 1, j).c_str());
      dxVsYhistograms[i][j]->Draw("E");
    }
    c.SaveAs("residuals_AO2D.pdf");
  }

  double xv[10] = {1, 2, 3, 4, 5, 6, 7, 8, 9, 10};

  PlotZTrend(10, xv, meanDx, "#Delta(x) vs. chamber;chamber; #Delta(x) (cm)", -5.0, 5.0, c);
  PlotZTrend(10, xv, meanDy, "#Delta(y) vs. chamber;chamber; #Delta(y) (cm)", -5.0, 5.0, c);

  c.Clear();
  c.SaveAs("residuals_AO2D.pdf)");
}

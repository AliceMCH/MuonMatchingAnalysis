#include <filesystem>
#include <fstream>
#include <algorithm>
#include <string>
#include <set>

#include "nlohmann/json.hpp"
using json = nlohmann::json;

void plot_DE_position()
{
  gStyle->SetOptStat(0);
  gStyle->SetOptFit(1111);
  gStyle->SetPalette(57, 0);
  gStyle->SetNumberContours(40);

  std::map<int, int> CH1L{{102, 0}, {101, 1}};
  std::map<int, int> CH1R{{103, 0}, {100, 1}};
  std::array<std::map<int, int>, 2> CH1{CH1L, CH1R};

  std::map<int, int> CH2L{{202, 0}, {201, 1}};
  std::map<int, int> CH2R{{203, 0}, {200, 1}};
  std::array<std::map<int, int>, 2> CH2{CH2L, CH2R};

  std::map<int, int> CH3L{{302, 0}, {301, 1}};
  std::map<int, int> CH3R{{303, 0}, {300, 1}};
  std::array<std::map<int, int>, 2> CH3{CH3L, CH3R};

  std::map<int, int> CH4L{{402, 0}, {401, 1}};
  std::map<int, int> CH4R{{403, 0}, {400, 1}};
  std::array<std::map<int, int>, 2> CH4{CH4L, CH4R};

  std::map<int, int> CH5L{{513, 0}, {512, 1}, {511, 2}, {510, 3}, {509, 4}, {508, 5}, {507, 6}, {506, 7}, {505, 8}};
  std::map<int, int> CH5R{{514, 0}, {515, 1}, {516, 2}, {517, 3}, {500, 4}, {501, 5}, {502, 6}, {503, 7}, {504, 8}};
  std::array<std::map<int, int>, 2> CH5{CH5L, CH5R};

  std::map<int, int> CH6L{{613, 0}, {612, 1}, {611, 2}, {610, 3}, {609, 4}, {608, 5}, {607, 6}, {606, 7}, {605, 8}};
  std::map<int, int> CH6R{{614, 0}, {615, 1}, {616, 2}, {617, 3}, {600, 4}, {601, 5}, {602, 6}, {603, 7}, {604, 8}};
  std::array<std::map<int, int>, 2> CH6{CH6L, CH6R};

  std::map<int, int> CH7L{{719, 0}, {718, 1}, {717, 2}, {716, 3}, {715, 4}, {714, 5}, {713, 6}, {712, 7}, {711, 8}, {710, 9}, {709, 10}, {708, 11}, {707, 12}};
  std::map<int, int> CH7R{{720, 0}, {721, 1}, {722, 2}, {723, 3}, {724, 4}, {725, 5}, {700, 6}, {701, 7}, {702, 8}, {703, 9}, {704, 10}, {705, 11}, {706, 12}};
  std::array<std::map<int, int>, 2> CH7{CH7L, CH7R};

  std::map<int, int> CH8L{{819, 0}, {818, 1}, {817, 2}, {816, 3}, {815, 4}, {814, 5}, {813, 6}, {812, 7}, {811, 8}, {810, 9}, {809, 10}, {808, 11}, {807, 12}};
  std::map<int, int> CH8R{{820, 0}, {821, 1}, {822, 2}, {823, 3}, {824, 4}, {825, 5}, {800, 6}, {801, 7}, {802, 8}, {803, 9}, {804, 10}, {805, 11}, {806, 12}};
  std::array<std::map<int, int>, 2> CH8{CH8L, CH8R};

  std::map<int, int> CH9L{{919, 0}, {918, 1}, {917, 2}, {916, 3}, {915, 4}, {914, 5}, {913, 6}, {912, 7}, {911, 8}, {910, 9}, {909, 10}, {908, 11}, {907, 12}};
  std::map<int, int> CH9R{{920, 0}, {921, 1}, {922, 2}, {923, 3}, {924, 4}, {925, 5}, {900, 6}, {901, 7}, {902, 8}, {903, 9}, {904, 10}, {905, 11}, {906, 12}};
  std::array<std::map<int, int>, 2> CH9{CH9L, CH9R};

  std::map<int, int> CH10L{{1019, 0}, {1018, 1}, {1017, 2}, {1016, 3}, {1015, 4}, {1014, 5}, {1013, 6}, {1012, 7}, {1011, 8}, {1010, 9}, {1009, 10}, {1008, 11}, {1007, 12}};
  std::map<int, int> CH10R{{1020, 0}, {1021, 1}, {1022, 2}, {1023, 3}, {1024, 4}, {1025, 5}, {1000, 6}, {1001, 7}, {1002, 8}, {1003, 9}, {1004, 10}, {1005, 11}, {1006, 12}};
  std::array<std::map<int, int>, 2> CH10{CH10L, CH10R};

  std::array<std::array<std::map<int, int>, 2>, 10> chambers{CH1, CH2, CH3, CH4, CH5, CH6, CH7, CH8, CH9, CH10};

  //------------------------------------------

  std::vector<std::array<double, 3>> CH1Lpos(CH1L.size());
  std::vector<std::array<double, 3>> CH1Rpos(CH1R.size());
  std::array<std::vector<std::array<double, 3>>, 2> CH1pos{CH1Lpos, CH1Rpos};

  std::vector<std::array<double, 3>> CH2Lpos(CH2L.size());
  std::vector<std::array<double, 3>> CH2Rpos(CH2R.size());
  std::array<std::vector<std::array<double, 3>>, 2> CH2pos{CH2Lpos, CH2Rpos};

  std::vector<std::array<double, 3>> CH3Lpos(CH3L.size());
  std::vector<std::array<double, 3>> CH3Rpos(CH3R.size());
  std::array<std::vector<std::array<double, 3>>, 2> CH3pos{CH3Lpos, CH3Rpos};

  std::vector<std::array<double, 3>> CH4Lpos(CH4L.size());
  std::vector<std::array<double, 3>> CH4Rpos(CH4R.size());
  std::array<std::vector<std::array<double, 3>>, 2> CH4pos{CH4Lpos, CH4Rpos};

  std::vector<std::array<double, 3>> CH5Lpos(CH5L.size());
  std::vector<std::array<double, 3>> CH5Rpos(CH5R.size());
  std::array<std::vector<std::array<double, 3>>, 2> CH5pos{CH5Lpos, CH5Rpos};

  std::vector<std::array<double, 3>> CH6Lpos(CH6L.size());
  std::vector<std::array<double, 3>> CH6Rpos(CH6R.size());
  std::array<std::vector<std::array<double, 3>>, 2> CH6pos{CH6Lpos, CH6Rpos};

  std::vector<std::array<double, 3>> CH7Lpos(CH7L.size());
  std::vector<std::array<double, 3>> CH7Rpos(CH7R.size());
  std::array<std::vector<std::array<double, 3>>, 2> CH7pos{CH7Lpos, CH7Rpos};

  std::vector<std::array<double, 3>> CH8Lpos(CH8L.size());
  std::vector<std::array<double, 3>> CH8Rpos(CH8R.size());
  std::array<std::vector<std::array<double, 3>>, 2> CH8pos{CH8Lpos, CH8Rpos};

  std::vector<std::array<double, 3>> CH9Lpos(CH9L.size());
  std::vector<std::array<double, 3>> CH9Rpos(CH9R.size());
  std::array<std::vector<std::array<double, 3>>, 2> CH9pos{CH9Lpos, CH9Rpos};

  std::vector<std::array<double, 3>> CH10Lpos(CH10L.size());
  std::vector<std::array<double, 3>> CH10Rpos(CH10R.size());
  std::array<std::vector<std::array<double, 3>>, 2> CH10pos{CH10Lpos, CH10Rpos};

  std::array<std::array<std::vector<std::array<double, 3>>, 2>, 10> chambersPos{CH1pos, CH2pos, CH3pos, CH4pos, CH5pos, CH6pos, CH7pos, CH8pos, CH9pos, CH10pos};
  std::array<std::array<std::vector<std::array<double, 3>>, 2>, 10> chambersPosIdeal{CH1pos, CH2pos, CH3pos, CH4pos, CH5pos, CH6pos, CH7pos, CH8pos, CH9pos, CH10pos};

  std::ifstream fGeom("geom.json");
  auto jGeom = json::parse(fGeom);

  // reference runs
  if (jGeom.count("alignables") > 0) {
    auto alignables = jGeom.at("alignables");
    for (const auto& alignable : alignables) {
      if (alignable.count("deid") < 1) continue;
      auto deId = alignable.at("deid").get<int>();
      auto& transform = alignable.at("transform");
      double tx = transform.at("tx").get<double>();
      std::cout << std::format("DE{} TX {}\n", deId, tx);
      double ty = transform.at("ty").get<double>();
      std::cout << std::format("DE{} TY {}\n", deId, ty);
      double tz = transform.at("tz").get<double>();
      std::cout << std::format("DE{} TZ {}\n", deId, tz);

      int chamberId = deId / 100 - 1;
      //if (chamberId < 0 || chamberId > 1) continue;
      if (chamberId < 0) continue;

      int LR = -1;
      if (chambers[chamberId][0].count(deId) > 0) LR = 0;
      if (chambers[chamberId][1].count(deId) > 0) LR = 1;

      if (LR < 0) continue;

      int deIndex = chambers[chamberId][LR][deId];

      chambersPos[chamberId][LR][deIndex][0] = tx;
      chambersPos[chamberId][LR][deIndex][1] = ty;
      chambersPos[chamberId][LR][deIndex][2] = tz;
    }
  } else {
    std::cout << "Key \"" << "alignables" << "\" not found in configuration" << std::endl;
  }

  std::ifstream fGeomIdeal("geom-ideal.json");
  jGeom = json::parse(fGeomIdeal);

  // reference runs
  if (jGeom.count("alignables") > 0) {
    auto alignables = jGeom.at("alignables");
    for (const auto& alignable : alignables) {
      if (alignable.count("deid") < 1) continue;
      auto deId = alignable.at("deid").get<int>();
      auto& transform = alignable.at("transform");
      double tx = transform.at("tx").get<double>();
      std::cout << std::format("DE{} TX {}\n", deId, tx);
      double ty = transform.at("ty").get<double>();
      std::cout << std::format("DE{} TY {}\n", deId, ty);
      double tz = transform.at("tz").get<double>();
      std::cout << std::format("DE{} TZ {}\n", deId, tz);

      int chamberId = deId / 100 - 1;
      //if (chamberId < 0 || chamberId > 1) continue;
      if (chamberId < 0) continue;

      int LR = -1;
      if (chambers[chamberId][0].count(deId) > 0) LR = 0;
      if (chambers[chamberId][1].count(deId) > 0) LR = 1;

      if (LR < 0) continue;

      int deIndex = chambers[chamberId][LR][deId];

      chambersPosIdeal[chamberId][LR][deIndex][0] = tx;
      chambersPosIdeal[chamberId][LR][deIndex][1] = ty;
      chambersPosIdeal[chamberId][LR][deIndex][2] = tz;
    }
  } else {
    std::cout << "Key \"" << "alignables" << "\" not found in configuration" << std::endl;
  }

  std::array<TGraph*, 3> grCH1L{new TGraph(), new TGraph(), new TGraph()};
  std::array<TGraph*, 3> grCH1Ld{new TGraph(), new TGraph(), new TGraph()};
  std::array<TGraph*, 3> grCH1R{new TGraph(), new TGraph(), new TGraph()};
  std::array<TGraph*, 3> grCH1Rd{new TGraph(), new TGraph(), new TGraph()};
  std::array<std::array<TGraph*, 3>, 4> grCH1{grCH1L, grCH1R, grCH1Ld, grCH1Rd};

  std::array<TGraph*, 3> grCH2L{new TGraph(), new TGraph(), new TGraph()};
  std::array<TGraph*, 3> grCH2Ld{new TGraph(), new TGraph(), new TGraph()};
  std::array<TGraph*, 3> grCH2R{new TGraph(), new TGraph(), new TGraph()};
  std::array<TGraph*, 3> grCH2Rd{new TGraph(), new TGraph(), new TGraph()};
  std::array<std::array<TGraph*, 3>, 4> grCH2{grCH2L, grCH2R, grCH2Ld, grCH2Rd};

  std::array<TGraph*, 3> grCH3L{new TGraph(), new TGraph(), new TGraph()};
  std::array<TGraph*, 3> grCH3Ld{new TGraph(), new TGraph(), new TGraph()};
  std::array<TGraph*, 3> grCH3R{new TGraph(), new TGraph(), new TGraph()};
  std::array<TGraph*, 3> grCH3Rd{new TGraph(), new TGraph(), new TGraph()};
  std::array<std::array<TGraph*, 3>, 4> grCH3{grCH3L, grCH3R, grCH3Ld, grCH3Rd};

  std::array<TGraph*, 3> grCH4L{new TGraph(), new TGraph(), new TGraph()};
  std::array<TGraph*, 3> grCH4Ld{new TGraph(), new TGraph(), new TGraph()};
  std::array<TGraph*, 3> grCH4R{new TGraph(), new TGraph(), new TGraph()};
  std::array<TGraph*, 3> grCH4Rd{new TGraph(), new TGraph(), new TGraph()};
  std::array<std::array<TGraph*, 3>, 4> grCH4{grCH4L, grCH4R, grCH4Ld, grCH4Rd};

  std::array<TGraph*, 3> grCH5L{new TGraph(), new TGraph(), new TGraph()};
  std::array<TGraph*, 3> grCH5Ld{new TGraph(), new TGraph(), new TGraph()};
  std::array<TGraph*, 3> grCH5R{new TGraph(), new TGraph(), new TGraph()};
  std::array<TGraph*, 3> grCH5Rd{new TGraph(), new TGraph(), new TGraph()};
  std::array<std::array<TGraph*, 3>, 4> grCH5{grCH5L, grCH5R, grCH5Ld, grCH5Rd};

  std::array<TGraph*, 3> grCH6L{new TGraph(), new TGraph(), new TGraph()};
  std::array<TGraph*, 3> grCH6Ld{new TGraph(), new TGraph(), new TGraph()};
  std::array<TGraph*, 3> grCH6R{new TGraph(), new TGraph(), new TGraph()};
  std::array<TGraph*, 3> grCH6Rd{new TGraph(), new TGraph(), new TGraph()};
  std::array<std::array<TGraph*, 3>, 4> grCH6{grCH6L, grCH6R, grCH6Ld, grCH6Rd};

  std::array<TGraph*, 3> grCH7L{new TGraph(), new TGraph(), new TGraph()};
  std::array<TGraph*, 3> grCH7Ld{new TGraph(), new TGraph(), new TGraph()};
  std::array<TGraph*, 3> grCH7R{new TGraph(), new TGraph(), new TGraph()};
  std::array<TGraph*, 3> grCH7Rd{new TGraph(), new TGraph(), new TGraph()};
  std::array<std::array<TGraph*, 3>, 4> grCH7{grCH7L, grCH7R, grCH7Ld, grCH7Rd};

  std::array<TGraph*, 3> grCH8L{new TGraph(), new TGraph(), new TGraph()};
  std::array<TGraph*, 3> grCH8Ld{new TGraph(), new TGraph(), new TGraph()};
  std::array<TGraph*, 3> grCH8R{new TGraph(), new TGraph(), new TGraph()};
  std::array<TGraph*, 3> grCH8Rd{new TGraph(), new TGraph(), new TGraph()};
  std::array<std::array<TGraph*, 3>, 4> grCH8{grCH8L, grCH8R, grCH8Ld, grCH8Rd};

  std::array<TGraph*, 3> grCH9L{new TGraph(), new TGraph(), new TGraph()};
  std::array<TGraph*, 3> grCH9Ld{new TGraph(), new TGraph(), new TGraph()};
  std::array<TGraph*, 3> grCH9R{new TGraph(), new TGraph(), new TGraph()};
  std::array<TGraph*, 3> grCH9Rd{new TGraph(), new TGraph(), new TGraph()};
  std::array<std::array<TGraph*, 3>, 4> grCH9{grCH9L, grCH9R, grCH9Ld, grCH9Rd};

  std::array<TGraph*, 3> grCH10L{new TGraph(), new TGraph(), new TGraph()};
  std::array<TGraph*, 3> grCH10Ld{new TGraph(), new TGraph(), new TGraph()};
  std::array<TGraph*, 3> grCH10R{new TGraph(), new TGraph(), new TGraph()};
  std::array<TGraph*, 3> grCH10Rd{new TGraph(), new TGraph(), new TGraph()};
  std::array<std::array<TGraph*, 3>, 4> grCH10{grCH10L, grCH10R, grCH10Ld, grCH10Rd};

  std::array<std::array<std::array<TGraph*, 3>, 4>, 10> graphs{grCH1, grCH2, grCH3, grCH4, grCH5, grCH6, grCH7, grCH8, grCH9, grCH10};

  std::array<std::string, 3> axes{"x", "y", "z"};
  for (int chamberId = 0; chamberId < chambersPos.size(); chamberId++) {
    for (int i = 0; i < 2; i++) {
      for (int k = 0; k < 3; k++) {
        TGraph* gr = graphs[chamberId][i][k];
        std::string title = std::format("CH{}{} DE {}", chamberId + 1, (i==0 ? "L" : "R"), axes[k]);
        gr->SetTitle(title.c_str());
        TGraph* grd = graphs[chamberId][i + 2][k];
        title = std::format("CH{}{} DE #Delta{}", chamberId + 1, (i==0 ? "L" : "R"), axes[k]);
        grd->SetTitle(title.c_str());
        for (int j = 0; j < chambersPos[chamberId][i].size(); j++) {
          double pos = chambersPos[chamberId][i][j][k];
          double posIdeal = chambersPosIdeal[chamberId][i][j][k];
          std::cout << std::format("CH{}-{}-{}  {}={}  {}Ideal={}\n", chamberId, i, j, axes[k], pos, axes[k], posIdeal);
          gr->AddPoint(j, pos);

          grd->AddPoint(j, pos - posIdeal);

          //if (j == 0) continue;
          //double tyPrev = chambersPos[chamberId][i][j - 1].second;
          //grd->AddPoint(j, ty - tyPrev);
        }
      }
    }
  }

  //return;

  TCanvas c("c", "c", 1200, 1200);
  c.Divide(2, 2);
  for (int chamberId = 0; chamberId < chambersPos.size(); chamberId++) {
    for (int k = 0; k < 3; k++) {
      for (int i = 0; i < 4; i++) {
        c.cd(i + 1);
        graphs[chamberId][i][k]->Draw("AL*");
        if (i > 1) {
          if (k < 2) {
            graphs[chamberId][i][k]->SetMinimum(-2);
            graphs[chamberId][i][k]->SetMaximum(2);
          //} else {
            //graphs[chamberId][i][k]->SetMinimum(-10);
            //graphs[chamberId][i][k]->SetMaximum(10);
          }
        }
      }
      if (chamberId == 0 && k == 0) c.SaveAs("c.pdf(");
      else c.SaveAs("c.pdf");
    }
  }

  c.Clear();
  c.SaveAs("c.pdf)");
}

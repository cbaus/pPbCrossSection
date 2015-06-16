#ifndef _include_modelInfo_h_
#define _include_modelInfo_h_

#include "TColor.h"
#include <vector>
#include <string>

class modelInfoClass{

 public:
  std::vector<std::string> models;
  std::vector<std::string> names;
  std::vector<int> colors;
  std::vector<int> styles;

  modelInfoClass()
    {
      models.push_back("PythiaMonash"); names.push_back("Pythia8 Monash Tune"); styles.push_back(10); colors.push_back(kGreen+2);
      models.push_back("PythiaZ2Star"); names.push_back("Pythia6 Z2*"); styles.push_back(9); colors.push_back(kBlue);
      models.push_back("PythiaMBR"); names.push_back("Pythia8 MBR Tune"); styles.push_back(8); colors.push_back(kRed);
      // models.push_back("Epos"); names.push_back("EPOS-LHC"); styles.push_back(7); colors.push_back(kRed);
      // models.push_back("QGSJetII"); names.push_back("QGSJETII-04"); styles.push_back(6); colors.push_back(kOrange);
    }
  virtual ~modelInfoClass() {}

  ClassDef(modelInfoClass, 1);
};
#endif

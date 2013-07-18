#include <iomanip>

void makeTable_noise()
{
  TFile f("plots/hf_cuts_noise.root");
  TVectorD* hf_m_cuts_light_noise = NULL;
  TVectorD* hf_p_cuts_light_noise = NULL;
  TVectorD* hf_m_cuts_medium_noise = NULL;
  TVectorD* hf_p_cuts_medium_noise = NULL;
  TVectorD* hf_m_cuts_heavy_noise = NULL;
  TVectorD* hf_p_cuts_heavy_noise = NULL;
  hf_m_cuts_light_noise = (TVectorD*)f.Get("hf_m_cuts_light_noise");
  hf_p_cuts_light_noise = (TVectorD*)f.Get("hf_p_cuts_light_noise");
  hf_m_cuts_medium_noise = (TVectorD*)f.Get("hf_m_cuts_medium_noise");
  hf_p_cuts_medium_noise = (TVectorD*)f.Get("hf_p_cuts_medium_noise");
  hf_m_cuts_heavy_noise = (TVectorD*)f.Get("hf_m_cuts_heavy_noise");
  hf_p_cuts_heavy_noise = (TVectorD*)f.Get("hf_p_cuts_heavy_noise");
  if(!hf_m_cuts_light_noise || !hf_p_cuts_light_noise) {cerr << "error1" << endl; return;}
  if(!hf_m_cuts_medium_noise || !hf_p_cuts_medium_noise) {cerr << "error2" << endl; return;}
  if(!hf_m_cuts_heavy_noise || !hf_p_cuts_heavy_noise) {cerr << "error3" << endl; return;}


  int ieta = -41;
  for (int i=0; i<12; i++)
    {
      cout << setprecision(3)
           << ieta++ << " & "
           << (*hf_m_cuts_light_noise)[i] << " & "
           << (*hf_m_cuts_medium_noise)[i] << " & "
           << (*hf_m_cuts_heavy_noise)[i] << " \\\\\\hline"
           << endl;
    }
  ieta = 30;
  for (int i=0; i<12; i++)
    {
      cout << setprecision(3)
           << ieta++ << " & "
           << (*hf_p_cuts_light_noise)[i] << " & "
           << (*hf_p_cuts_medium_noise)[i] << " & "
           << (*hf_p_cuts_heavy_noise)[i] << " \\\\\\hline"
           << endl;
    }

}

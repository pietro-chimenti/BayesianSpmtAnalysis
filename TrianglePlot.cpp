#include <iostream>

#include "TCanvas.h"
#include "TFile.h"
#include "TTree.h" 
#include "TStyle.h"

using namespace std;

#ifndef __CINT__
int main()
{
  cout << " Drawing triangle plot "<< endl;

  gStyle->SetOptStat(0);

  TCanvas c1("TP","Triangle Plot",700,700);

  TFile f("BayesianSpmtFit_mcmc.root");
  f.ls();

  TTree* tp;
  f.GetObject("BayesianSpmtFit_mcmc",tp);

  c1.Divide(3,3);
  c1.cd(1);

  tp->Print();
  tp->Draw("s2t12","phase==1");

  c1.cd(4);
  tp->Draw("DelM2_21:s2t12","phase==1");

  c1.cd(5);
  tp->Draw("DelM2_21","phase==1");

  c1.cd(7);
  tp->Draw("s2t13:s2t12","phase==1");

  c1.cd(8);
  tp->Draw("DelM2_21:s2t13","phase==1");

  c1.cd(9);
  tp->Draw("s2t13","phase==1");


  c1.Draw();
  c1.SaveAs("TrianglePlot.pdf");
  c1.SaveAs("TrianglePlot.jpg");

  return 0;
}
#endif

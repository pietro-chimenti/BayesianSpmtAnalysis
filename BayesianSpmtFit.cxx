#include "BayesianSpmtFit.h"

#include <iostream>

#include "TFile.h"
#include "TH1F.h"
#include "TCanvas.h"
#include "TMath.h"

#include <BAT/BCLog.h>

// ----------------------------------------------------------------------------
BayesianSpmtFit::BayesianSpmtFit(const std::string& name, BayesianSpmtConfig& config)
    : BCModel(name), book(), myConfig(config)
{
  // define parameters, observables, priors and set initial values

  if(myConfig.getInt(std::string("par_s2t12"))){
    AddParameter("s2t12", 0.29, 0.33, "#s2t12", "");
    GetParameter("s2t12").SetPriorConstant();
    book.push_back("s2t12");
  }

  if(myConfig.getInt(std::string("par_DelM2_21"))){
    AddParameter("DelM2_21",7.1e-5,7.5e-5,"#DelM2_21","");
    GetParameter("DelM2_21").SetPriorConstant();
    book.push_back("DelM2_21");
  }

  if(myConfig.getInt(std::string("par_s2t13"))){
    AddParameter("s2t13",0.01,0.04,"#s2t13","");
    GetParameter("s2t13").SetPriorConstant();
    book.push_back("s2t13");
  }

  AddObservable("s22t12", 0., 1., "#s22t12", "");


  // initialize oscilation parameters
  InitOscPars();
  LoadSimTree();
  LoadBins();
  LoadSpectrumExp();
//  PlotNuFit();

}

// ----------------------------------------------------------------------------
BayesianSpmtFit::~BayesianSpmtFit()
{
}

// ----------------------------------------------------------------------------
void BayesianSpmtFit::LoadParameters(const std::vector<double>& pars)
{

  for( unsigned int i=0; i < pars.size(); ++i ){ 
    if( book[i] == "s2t12" )    s2t12    = pars[i];
    if( book[i] == "DelM2_21" ) DelM2_21 = pars[i];
    if( book[i] == "s2t13" )    s2t13    = pars[i];
  }
}

// ----------------------------------------------------------------------------
double BayesianSpmtFit::LogLikelihood(const std::vector<double>& pars)
{
  double ll=0;

  LoadParameters(pars);
  LoadSpectrumTh();

  // calculate terms related to covariance matrix
  for(unsigned int i = 0; i<spectrum_exp.size(); ++i)
    for(unsigned int j = 0; j<spectrum_exp.size(); ++j)
      ll+=(spectrum_exp[i]-spectrum_th[i])*M_inv[i][j]*(spectrum_exp[j]-spectrum_th[j]);

  // now add pulls
  if(myConfig.getInt(std::string("par_s2t13"))){
   ll+=( (s2t13-pull_s2t13)*(s2t13-pull_s2t13) )/( pull_s2t13_err*pull_s2t13_err );
  }

  ll*=-0.5;

  return ll;
}

// ----------------------------------------------------------------------------
//double BayesianSpmtFit::LogAPrioriProbability(const std::vector<double>& pars)
//{
//  return 0;
//}

// ----------------------------------------------------------------------------
void BayesianSpmtFit::CalculateObservables(const std::vector<double>& pars)
{
  LoadParameters(pars);
  GetObservable("s22t12") = 4*s2t12*(1-s2t12);
}

// ----------------------------------------------------------------------------
void BayesianSpmtFit::LoadSimTree()
{
  std::string filename=myConfig.getString("simDataPath")+"/"+myConfig.getString("simDataFile");
  TFile f(filename.c_str());
  simTree = (TTree*) f.Get(myConfig.getString("simTreeName").c_str());
  //simTree->Print();

  Double_t NeutrinoEnergy_Th;
  Double_t NeutrinoDistance_Th;
  Double_t PromptEvisID;
  simTree->SetBranchAddress("myNeutrinoEnergy_Th",&NeutrinoEnergy_Th);
  simTree->SetBranchAddress("myNeutrinoDistance_Th",&NeutrinoDistance_Th);
  simTree->SetBranchAddress("myPromptEvisID",&PromptEvisID);
  for(int i=0; i<simTree->GetEntries(); i++){
    simTree->GetEntry(i);
    vectorELP.push_back(std::make_tuple(NeutrinoEnergy_Th,NeutrinoDistance_Th,PromptEvisID));
  }

}

// ----------------------------------------------------------------------------
void BayesianSpmtFit::LoadBins()
{
  int nBins = myConfig.getInt("nBins");  
  for(int i = 0; i<=nBins; i++){
    std::string binName = std::string("bin_")+std::to_string(i);
    double binLimit = myConfig.getDouble(binName);
    binLimits.push_back(binLimit); 
    if(i>0){
      spectrum_th.push_back(0);
      spectrum_exp.push_back(0);
    }
  }
}

// ----------------------------------------------------------------------------
double BayesianSpmtFit::Pee(double E, double L)
{
  double oscFact = myConfig.getDouble("oscFact");

  double prob = 1;

  double s22t12 = 4 * s2t12*(1-s2t12);
  double s22t13 = 4 * s2t13*(1-s2t13);
  double c4t13=(1-s2t13)*(1-s2t13);
  double c2t12=(1-s2t12);

  double Del21=oscFact*DelM2_21*L/E;
  double Del31=oscFact*DelM2_31*L/E;
  double Del32=oscFact*(DelM2_31-DelM2_21)*L/E;

  double fact1 = s22t12 * sin(Del21)*sin(Del21)*c4t13;
  double fact2 = s22t13 * (c2t12*sin(Del31)*sin(Del31)+s2t12*sin(Del32)*sin(Del32));
 
  prob = prob - fact1 - fact2; 

  return prob;
}

// ----------------------------------------------------------------------------
int BayesianSpmtFit::get_bin(double e)
{
  for(unsigned int i = 0; i <= binLimits.size(); ++i ) 
    if(e<binLimits[i]) return i;
  return binLimits.size()+1;
}

// ----------------------------------------------------------------------------
void BayesianSpmtFit::set_m_inv()
{
  set_m_total();

  M_inv.ResizeTo(spectrum_exp.size(),spectrum_exp.size()); 
  M_inv = M_total.Invert();
}

// ----------------------------------------------------------------------------
void BayesianSpmtFit::set_m_total()
{
  set_m_stat();
  set_m_norm();
  M_total.ResizeTo(spectrum_exp.size(),spectrum_exp.size());
  M_total = M_stat;
  if( myConfig.getInt("ll_norm") ) M_total += M_norm;
  if( myConfig.getInt("ll_b2b") ) M_total += M_b2b;
}

// ----------------------------------------------------------------------------
void BayesianSpmtFit::set_m_stat()
{
  M_stat.ResizeTo(spectrum_exp.size(),spectrum_exp.size());
  for(unsigned int i = 0; i<spectrum_exp.size(); ++i)
    for(unsigned int j = 0; j<spectrum_exp.size(); ++j)
      if(i==j)M_stat[i][j]=spectrum_exp[j];
      else M_stat[i][j]=0;
}

// ----------------------------------------------------------------------------
void BayesianSpmtFit::set_m_norm()
{
  double norm_error = myConfig.getDouble("norm_error");
  M_norm.ResizeTo(spectrum_exp.size(),spectrum_exp.size());
  for(unsigned int i = 0; i<spectrum_exp.size(); ++i)
    for(unsigned int j = 0; j<spectrum_exp.size(); ++j)
      M_norm[i][j]=norm_error*norm_error*spectrum_exp[i]*spectrum_exp[j];  
}

// ----------------------------------------------------------------------------
void BayesianSpmtFit::set_m_b2b()
{
  double b2b_error = myConfig.getDouble("b2b_error");
  M_b2b.ResizeTo(spectrum_exp.size(),spectrum_exp.size());
  for(unsigned int i = 0; i<spectrum_exp.size(); ++i)
    for(unsigned int j = 0; j<spectrum_exp.size(); ++j)
      if( i == j ) M_b2b[i][j] = b2b_error*b2b_error*spectrum_exp[i]*spectrum_exp[j];  
      else M_b2b[i][j] = 0;
}

// ----------------------------------------------------------------------------
void BayesianSpmtFit::LoadSpectrumTh()
{

  for(double &i: spectrum_th){ i=0;}
  double nuE;
  double nuL;
  double visE;
  int index=0;
  for(auto i : vectorELP){
    std::tie(nuE,nuL,visE) = i;
    index=get_bin(visE); 
    if(index==0||index==binLimits.size()+1) continue;
    spectrum_th[index]+=Pee(nuE,nuL);
  }

  // now normalize:
  for(unsigned int i =0; i<spectrum_th.size(); ++i) spectrum_th[i]*=normalization;
 

}

// ----------------------------------------------------------------------------
void BayesianSpmtFit::LoadSpectrumExp()
{

  for(double &i: spectrum_exp){ i=0;}
  double nuE;
  double nuL;
  double visE;
  int index=0;
  for(auto i : vectorELP){
    std::tie(nuE,nuL,visE) = i;
    index=get_bin(visE); 
    if(index==0||index==binLimits.size()+1) continue;
    spectrum_exp[index]+=Pee(nuE,nuL);
  }

  // now normalize to total number of events:
  tot_events = myConfig.getInt("tot_meas_events");
  double sum=0;
  for(unsigned int i =0; i<spectrum_exp.size(); ++i) sum+=spectrum_exp[i];
  normalization=tot_events/sum;
  for(unsigned int i =0; i<spectrum_exp.size(); ++i) spectrum_exp[i]*=normalization;
  for(unsigned int i =0; i<spectrum_exp.size(); ++i) if(spectrum_exp[i]==0)spectrum_exp[i]+=1e-3; // just to avoid a singular error matrix

  set_m_inv();
}

// ----------------------------------------------------------------------------
void BayesianSpmtFit::InitOscPars()
{
  s2t12=myConfig.getDouble("init_s2t12");
  s2t23=myConfig.getDouble("init_s2t23");
  s2t13=myConfig.getDouble("init_s2t13");
  DelM2_21=myConfig.getDouble("init_DelM2_21");
  DelM2_31=myConfig.getDouble("init_DelM2_31");

  pull_s2t13 = myConfig.getDouble("pull_s2t13");
  pull_s2t13_err = myConfig.getDouble("pull_s2t13_err");

}

// ----------------------------------------------------------------------------
void BayesianSpmtFit::PlotNuFit()
{
  // load file first
  TH1F NuFitOscPlot("NuFitOscPlot","NuFit 4.0 B.F. oscilation",spectrum_exp.size(),binLimits.data());
  TCanvas ca;

  for(unsigned int i=1; i<=spectrum_exp.size(); ++i){
    NuFitOscPlot.SetBinContent(i,spectrum_exp[i-1]);
  }  

  NuFitOscPlot.Draw();
  ca.SaveAs("NuFitOscPlot.pdf");
}

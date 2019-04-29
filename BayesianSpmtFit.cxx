#include "BayesianSpmtFit.h"

#include <iostream>

#include "TFile.h"
#include "TH1F.h"
#include "TCanvas.h"
#include "TMath.h"

#include <BAT/BCLog.h>

const std::string BayesianSpmtFit::configFile = "config.txt";

// ----------------------------------------------------------------------------
BayesianSpmtFit::BayesianSpmtFit(const std::string& name)
    : BCModel(name)
{

  AddParameter("s2t12", 0.29, 0.33, "#s2t12", "");
  GetParameters().Back().SetPriorConstant();
 
  AddParameter("DelM2_21",7.1e-5,7.5e-5,"#DelM2_21","");
  GetParameters().Back().SetPriorConstant();

  //AddParameter("s2t13",0.01,0.04,"#s2t13","");
  //GetParameters().Back().SetPriorConstant();

  AddObservable("s22t12", 0.86, 0.87, "#s22t12", "");


}

// ----------------------------------------------------------------------------
BayesianSpmtFit::~BayesianSpmtFit()
{
}

// ----------------------------------------------------------------------------
double BayesianSpmtFit::LogLikelihood(const std::vector<double>& pars)
{
  static int call_number = 0; // to monitor calls
  ++call_number;
  std::cout << " call LL " << call_number << std::endl; 
  std::cout << pars[0] << " " << pars[1] << " " << pars[2] << std::endl; 

  double ll=0;

  // set parameters and theoretical spectrum
  s2t12 = pars[0];
  DelM2_21 = pars[1];
  s2t13 = pars[2];
  LoadSpectrumTh();

  // calculate residuals
  std::vector<double> residuals(spectrum_exp.size());
  for(unsigned int i = 0; i<spectrum_exp.size(); ++i)
    residuals[i]=(spectrum_exp[i]-spectrum_th[i]);

  // calculate terms related to covariance matrix
  for(unsigned int i = 0; i<spectrum_exp.size(); ++i)
    for(unsigned int j = 0; j<spectrum_exp.size(); ++j)
      ll+=(residuals[i])*M_inv[i][j]*(residuals[j]);

  // now add pulls
  double s2t13_nf = getDouble("s2t13_nf");
  double s2t13_err_nf = getDouble("s2t13_err_nf");
  //ll+=( (pars[2]-s2t13_nf)*(pars[2]-s2t13_nf) )/( s2t13_err_nf*s2t13_err_nf );

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
  GetObservable(0) = 4*pars[0]*(1-pars[0]);
}

// ----------------------------------------------------------------------------
void BayesianSpmtFit::Setup()
{
  LoadSimTree();
  LoadBins();
}

// ----------------------------------------------------------------------------
void BayesianSpmtFit::LoadSimTree()
{
  std::string filename=getString("simDataPath")+"/"+getString("simDataFile");
  TFile f(filename.c_str());
  simTree = (TTree*) f.Get(getString("simTreeName").c_str());
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
  int nBins = getInt("nBins");  
  for(int i = 0; i<=nBins; i++){
    std::string binName = std::string("bin_")+std::to_string(i);
    double binLimit = getDouble(binName);
    binLimits.push_back(binLimit); 
    if(i>0){
      spectrum_th.push_back(0);
      spectrum_exp.push_back(0);
    }
  }
}

// ----------------------------------------------------------------------------
void BayesianSpmtFit::LoadConfig()
{
  BCLog::OutSummary("Loading Config");
  std::string line;
  std::ifstream inFile(configFile.c_str());
  if(inFile.is_open()){
    while(getline(inFile,line)){
      bool isBlanck = true;
      for( std::string::const_iterator it = line.begin(); it!= line.end(); ++it){
        if(!isspace(*it)){
          isBlanck = false;
          break;
        }
      }
      if(isBlanck) continue;
      //BCLOG_SUMMARY(line.c_str());
      BCLog::OutSummary(line.c_str());
      
      std::istringstream iss(line);
      std::string tag;
      std::string key;
      iss >> tag;
      if(tag=="#") continue;
      if(tag=="int"){
        int value;
        iss >> key;
        iss >> value;
        intParams.insert(std::pair<std::string,int>(key,value));
        BCLog::OutDebug("int param");
      }
      if(tag=="double"){
        double value;
        iss >> key;
        iss >> value;
        doubleParams.insert(std::pair<std::string,double>(key,value));
        BCLog::OutDebug("double param");
      }
      if(tag=="string"){
        std::string value;
        iss >> key;
        iss >> value;
        stringParams.insert(std::pair<std::string,std::string>(key,value));
        BCLog::OutDebug("int param");
      }
        
    }
  }
}

// ----------------------------------------------------------------------------
int BayesianSpmtFit::getInt(std::string name)
{
  if(intParams.find(name)==intParams.end()) throw std::out_of_range("Unable to fine"+name);
  return intParams.find(name)->second;
}

// ----------------------------------------------------------------------------
double BayesianSpmtFit::getDouble(std::string name)
{
  if(doubleParams.find(name)==doubleParams.end()) throw std::out_of_range("Unable to fine"+name);
  return doubleParams.find(name)->second;
}

// ----------------------------------------------------------------------------
std::string BayesianSpmtFit::getString(std::string name)
{
  if(stringParams.find(name)==stringParams.end()) throw std::out_of_range("Unable to fine"+name);
  return stringParams.find(name)->second;
}

// ----------------------------------------------------------------------------
double BayesianSpmtFit::Pee(double E, double L)
{
  double oscFact = getDouble("oscFact");

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
  for(unsigned int i=0; i< binLimits.size(); ++i ) 
    if(e<binLimits[i]) return i;
  return binLimits.size();
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
  //M_total = M_stat+M_norm;
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
  double norm_error = getDouble("norm_error");
  M_norm.ResizeTo(spectrum_exp.size(),spectrum_exp.size());
  for(unsigned int i = 0; i<spectrum_exp.size(); ++i)
    for(unsigned int j = 0; j<spectrum_exp.size(); ++j)
      M_norm[i][j]=norm_error*norm_error*spectrum_exp[i]*spectrum_exp[j];  
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
    spectrum_th[index]+=Pee(nuE,nuL);
  }

  // now normalize to total number of events:
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
    spectrum_exp[index]+=Pee(nuE,nuL);
  }

  // now normalize to total number of events:
  tot_events = getInt("tot_meas_events");
  double sum=0;
  for(unsigned int i =0; i<spectrum_exp.size(); ++i) sum+=spectrum_exp[i];
  normalization=tot_events/sum;
  for(unsigned int i =0; i<spectrum_exp.size(); ++i) spectrum_exp[i]*=normalization;
  for(unsigned int i =0; i<spectrum_exp.size(); ++i) if(spectrum_exp[i]==0)spectrum_exp[i]+=1e-3; // just to avoid a singular error matrix

  set_m_inv();
}

// ----------------------------------------------------------------------------
void BayesianSpmtFit::SetNuFit_NO()
{
  s2t12=getDouble("s2t12_nf");
  s2t23=getDouble("s2t23_nf");
  s2t13=getDouble("s2t13_nf");
  DelM2_21=getDouble("DelM2_21_nf");
  DelM2_31=getDouble("DelM2_3l_nf");
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

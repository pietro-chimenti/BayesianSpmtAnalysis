#include <BAT/BCLog.h>

#include "BayesianSpmtFit.h"
#include "BayesianSpmtConfig.h"

int main( int argc, char *argv[] )
{

  if( argc != 2 ){
    std::cout << " usage: " << argv[0] << " <filename>\n";
    return -1;
  }

  BayesianSpmtConfig myConfig;
  myConfig.LoadConfig(std::string(argv[1]));

  // open log file
  BCLog::OpenLog(myConfig.getString("outpath") + "log.txt", BCLog::detail, BCLog::detail);


  BayesianSpmtFit m("BayesianSpmtFit",myConfig);

  m.SetNChains(myConfig.getInt("NC"));
  m.SetNIterationsPreRunCheck(myConfig.getInt("NIPRC"));
  m.SetNIterationsPreRunMin(myConfig.getInt("NIPRMI"));
  m.SetNIterationsPreRunMax(myConfig.getInt("NIPRMA"));
  m.SetNIterationsRun(myConfig.getInt("NIR"));

  BCLog::OutSummary("Test model created");
  m.WriteMarkovChain(myConfig.getString("outpath") + m.GetSafeName() + "_mcmc.root", "RECREATE");
  m.MarginalizeAll(BCIntegrate::kMargMetropolis);
  BCLog::OutSummary("Marginalized done");
  m.FindMode(m.GetBestFitParameters());
  BCLog::OutSummary("FindMode done");

  m.PrintAllMarginalized(myConfig.getString("outpath") + m.GetSafeName() + "_plots.pdf");

  m.PrintParameterPlot(myConfig.getString("outpath") + m.GetSafeName() + "_parameters.pdf");
  m.PrintCorrelationPlot(myConfig.getString("outpath") + m.GetSafeName() + "_correlation.pdf");
  m.PrintCorrelationMatrix(myConfig.getString("outpath") + m.GetSafeName() + "_correlationMatrix.pdf");
  m.PrintKnowledgeUpdatePlots(myConfig.getString("outpath") + m.GetSafeName() + "_update.pdf");

  m.PrintSummary();

  BCLog::OutSummary("Exiting");
  BCLog::CloseLog();

  return 0;
}

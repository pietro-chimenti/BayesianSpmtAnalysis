#include <BAT/BCLog.h>

#include "BayesianSpmtFit.h"
#include "BayesianSpmtConfig.h"

int main()
{
    // open log file
    BCLog::OpenLog("log.txt", BCLog::detail, BCLog::detail);

    BayesianSpmtConfig myConfig;
    myConfig.LoadConfig(std::string("config.txt"));

    BayesianSpmtFit m("BayesianSpmtFit",myConfig);

    m.SetNChains(1);
    m.SetNIterationsPreRunCheck(30);
    m.SetNIterationsPreRunMin(90);
    m.SetNIterationsPreRunMax(300);
    m.SetNIterationsRun(300);

    BCLog::OutSummary("Test model created");
    m.WriteMarkovChain(m.GetSafeName() + "_mcmc.root", "RECREATE");
    m.MarginalizeAll(BCIntegrate::kMargMetropolis);
    BCLog::OutSummary("Marginalized done");
    m.FindMode(m.GetBestFitParameters());
    BCLog::OutSummary("FindMode done");

    m.PrintAllMarginalized(m.GetSafeName() + "_plots.pdf");

    m.PrintParameterPlot(m.GetSafeName() + "_parameters.pdf");
    m.PrintCorrelationPlot(m.GetSafeName() + "_correlation.pdf");
    m.PrintCorrelationMatrix(m.GetSafeName() + "_correlationMatrix.pdf");
    m.PrintKnowledgeUpdatePlots(m.GetSafeName() + "_update.pdf");

    m.PrintSummary();

    BCLog::OutSummary("Exiting");
    BCLog::CloseLog();

    return 0;
}

#include <BAT/BCLog.h>

#include "BayesianSpmtFit.h"

int main()
{
    // open log file
    BCLog::OpenLog("log.txt", BCLog::detail, BCLog::detail);

    BayesianSpmtFit m("BayesianSpmtFit");
    m.LoadConfig();
    m.Setup();

    // show oscillation plot
    m.SetNuFit_NO();

    m.LoadSpectrumExp();
    m.PlotNuFit();
    m.SetNChains(1);
    m.SetNIterationsPreRunCheck(300);
    m.SetNIterationsPreRunMin(900);
    m.SetNIterationsPreRunMax(3000);
    m.SetNIterationsRun(3000);

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

#ifndef __BAYESIANSPMTFIT__H
#define __BAYESIANSPMTFIT__H

/**
 * @class BayesianSpmtFit
 * @brief A class used to evaluate the sensitivity 
 * of the SPMT system to solar parameters 
 * @author Pietro Chimenti
 * @version 0.1
 * @date 04.2019
 * @details This class implement a BCModel of the BAT package and is used to estimate
 * the sensitivity of the JUNO SPMT system to solar oscilation parameter.
 */

/*
 * Copyright (C) 2019, Pietro Chimenti
 * All rights reserved.
 *
 * For the licensing terms see COPYING.
 */

// ----------------------------------------------------------------------------
#include <BAT/BCModel.h>

#include "TTree.h"
#include "TMatrixDSym.h"

#include <string>
#include <vector>
#include <map>
#include <tuple>

#include "BayesianSpmtConfig.h"

// ----------------------------------------------------------------------------
class BayesianSpmtFit : public BCModel
{

public:

    BayesianSpmtFit(const std::string& name, BayesianSpmtConfig& config);

    ~BayesianSpmtFit();

    double LogLikelihood(const std::vector<double>& pars);

    // double LogAPrioriProbability(const std::vector<double> & pars);

    void CalculateObservables(const std::vector<double> & pars);

    /**
     * finally save plots */ 
    void SavePlots();

  private:

    /**
     * calculate energy scale error matrix */
    void set_m_escale();  

    /**
     * calculate flux error */
    void set_m_flux();

    /**
     * calculate the shape error matrix: fully uncorrelated approx */ 
    void set_m_b2b();

    /**
     * calculate the normalization error matrix: fully correlated */ 
    void set_m_norm();

    /**
     * calculate the statistical error matrix: diag(N_i) */
    void set_m_stat();

    /**
     * calculate the total error matrix: stat + norm */
    void set_m_total();

    /**
     * calculate the inverse error matrix */ 
    void set_m_inv();

    /**
     * here we load oscillation parameters from fit parameters */ 
    void LoadParameters(const std::vector<double>& pars);

    /**
     * here we load the simulation tree */
    void LoadSimTree();  

    /**
     * here we load the infos necessary to perform the fit:
     * fractional matrices, pull term infos etc. */
     void LoadFitInfo(); 

    /**
     * given the visible energy e return the corresponding bin number */
    int get_bin(double e);

    /**
     * calculate the oscillation probability */
    double Pee(double E, double L);

    /**
     * here we load the theoretical spectrum 
     * to be compared with the experimental one */ 
    void LoadSpectrumTh();

    /**
     * here we load the Mock simulation to generate Exp data */
    void LoadMockTree();  

    /**
     * here we load the bin limits */ 
    void LoadBins();    

    /**
     * here we load the simulated experimental spectrum
     * for now with the Asimov aproximation (no statistical fluctuation) */
    void LoadSpectrumExp();

    /**
     * here we set the data used for fit: either simulated or measured */
     void SetData(); 

    /**
     *  set initial values for oscillation parameters */
    void InitPars();

    /**
     * here se set the fit parameters, priors and observables */
    void SetParameters();

// ----------------------------------------------------------------------------
// now fields
//

    /**
     * used for bookeeping parameters */
    std::vector<std::string> book; 

    /**
     *  ttree of simulated IBD interactions */ 
    TTree * simTree;

    /**
     * vector vector of neutrino energy, distance and p.e. number 
     * used to calculate both exp and th spectra */ 
    std::vector<std::tuple<double,double,double>> vectorELP;

    /**
     *  vector of bin limits */
    std::vector<double> binLimits;

    /**
     * theoretical oscillated spectrum */ 
    std::vector<double> spectrum_th;

    /**
     * simulated experimental spectrum */ 
    std::vector<double> spectrum_exp;

    /// now oscilation parameters
    double s2t12; /// Sin^2(t_12)
    double s2t23; /// Sin^2(t_23)
    double s2t13; /// Sin^2(t_13)
    double DelM2_21; /// Delta M^2_21
    double DelM2_31; /// Delta M^2_31 - positive or negative depending on M.H.
    double normalization; 

    double tot_events; 

    /// now error matrixes and pull terms
    TMatrixDSym M_inv;        /// inverse ot the total error matrix
    TMatrixDSym M_total;      /// total error matrix: stat + norm
    TMatrixDSym M_stat;       /// statistical error matrix: diagonal
    TMatrixDSym M_norm;       /// normalization error matrix: fully correlated
    TMatrixDSym M_b2b;        /// shape error in the b2b approx.
    TMatrixDSym M_flux;       /// flux error from Thiago (H.-M. approx)
    TMatrixDSym M_flux_frac;  /// fractional flux error matrix
    TMatrixDSym M_EScale;     /// energy scale matrix from Thiago
    TMatrixDSym M_EScale_frac;/// energy scale matrix from thiago, fractional

    double pull_s2t13;
    double pull_s2t13_err;

    // Configuration interface

    BayesianSpmtConfig& myConfig;


};
// ---------------------------------------------------------

#endif

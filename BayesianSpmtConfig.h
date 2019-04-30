#ifndef __BAYESIANSPMTCONFIG__H
#define __BAYESIANSPMTCONFIG__H

/**
 * @class BayesianSpmtConfig
 * @brief A class used to set the fit parameters 
 * of the SPMT system to solar parameters 
 * @author Pietro Chimenti
 * @version 0.1
 * @date 04.2019
 * @details This class implement an interface between a configuration file 
 * and the fit code.
 */

/*
 * Copyright (C) 2019, Pietro Chimenti
 * All rights reserved.
 *
 * For the licensing terms see COPYING.
 */

// ----------------------------------------------------------------------------

#include <string>
#include <map>

// ----------------------------------------------------------------------------
class BayesianSpmtConfig
{

public:

    BayesianSpmtConfig();

    ~BayesianSpmtConfig();


    /**
     *  read configuration file and fill configuration maps */
    void LoadConfig(std::string configFile);

    /**
     * get the int value of the name configuration variable */
    int getInt(std::string name);

    /**
     * get the double value of the name configuration variable */ 
    double getDouble(std::string name);

    /**
     * get the string value of the name configuration variable */ 
    std::string getString(std::string name);

  private:

    /**
     * map of integer parameters */ 
    std::map<std::string, int> intParams;

    /**
     * map of double parameters */ 
    std::map<std::string, double> doubleParams;

    /**
     * map of string parameters */ 
    std::map<std::string, std::string> stringParams;
};
// ---------------------------------------------------------

#endif

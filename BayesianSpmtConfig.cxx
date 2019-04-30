#include "BayesianSpmtConfig.h"

#include <iostream>
#include <sstream>

#include <BAT/BCLog.h>

// ----------------------------------------------------------------------------
BayesianSpmtConfig::BayesianSpmtConfig()
{
}

// ----------------------------------------------------------------------------
BayesianSpmtConfig::~BayesianSpmtConfig()
{
}

// ----------------------------------------------------------------------------
void BayesianSpmtConfig::LoadConfig(std::string configFile)
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
        BCLog::OutDebug("string param");
      }
        
    }
  }
}

// ----------------------------------------------------------------------------
int BayesianSpmtConfig::getInt(std::string name)
{
  if(intParams.find(name)==intParams.end()) throw std::out_of_range("Unable to fine"+name);
  return intParams.find(name)->second;
}

// ----------------------------------------------------------------------------
double BayesianSpmtConfig::getDouble(std::string name)
{
  if(doubleParams.find(name)==doubleParams.end()) throw std::out_of_range("Unable to fine"+name);
  return doubleParams.find(name)->second;
}

// ----------------------------------------------------------------------------
std::string BayesianSpmtConfig::getString(std::string name)
{
  if(stringParams.find(name)==stringParams.end()) throw std::out_of_range("Unable to fine"+name);
  return stringParams.find(name)->second;
}

//-------------------------------------------------------
// Description:
//    Class for set cuts & vetoes
// Authors:
//    Chiara Rovelli & Emanuele Di Marco
//    Universita' di Roma "La Sapienza" & INFN Roma
//-------------------------------------------------------

#ifndef SELECTION_H
#define SELECTION_H

#include <vector>
#include <map>
#include <string>

class Selection {
  
public:
  Selection(std::string fileCuts, std::string fileSwitches);
  virtual ~Selection();
  // add a cut to the selection
  void addCut(std::string name);
  // add a switch to enable something from file
  void addSwitch(std::string name);

  // retrieve selection configuration
  std::pair<float,float> getCut(std::string name);
  int getSwitch(std::string name);

  // apply cut on a single var
  bool passCut(std::string name, float var);
  // apply the same cut on n vars (ex. n electrons)
  bool passCut(std::string name, std::vector<float> vars);

  // get cut map
  std::map<std::string,std::pair<float,float> > getSelection() { return _cut; };
  void summary();

private:
  std::pair<float,float> readIntervalFromFile(std::string name);
  int readSwitchFromFile(std::string name);
  std::string _fileCut, _fileSwitch;
  std::map<std::string,std::pair<float,float> > _cut;
  std::map<std::string,int> _switch;

  
};
#endif // SELECTION_H

#ifndef NEUTRINOS_Tools_Transition_Maps_H
#define NEUTRINOS_Tools_Transition_Maps_H

#include "ATOOLS/Phys/Flavour.H"
#include "ATOOLS/Org/Exception.H"
#include <iostream>
#include <vector>
#include <algorithm>
#include <list>
#include <map>
#include <set>

using namespace ATOOLS;
using namespace std;

namespace NEUTRINOS {
  vector<int> getQuarkContent(int x);
  std::string compareQuarkContent(vector<int> x1, vector<int> x2);
}

#endif
#include "NEUTRINOS++/Tools/Transition_Maps.H"

using namespace NEUTRINOS;
using namespace ATOOLS;
using namespace std;

namespace NEUTRINOS {
  vector<int> getQuarkContent(int x)
  {
      //Pull out quark content from kf_code
      vector<int> digits = {};
      for(int di=0; di<7; di++) {
        digits.push_back(x % 10);
        x /= 10;
        if (x == 0) break;
      }
      int length = digits.size();
      return {digits[length-4],digits[length-3],digits[length-2]};
  }

  std::string compareQuarkContent(vector<int> x1, vector<int> x2)
  {
    sort(x1.begin(), x1.end());
    sort(x2.begin(), x2.end());

    //Get Intersection of quark content
    vector<int> inter;
    set_intersection(
      x1.begin(), x1.end(),
      x2.begin(), x2.end(),
      back_inserter(inter)
    );

    //Remove mutual quark content
    x1.erase(
      set_difference(
        x1.begin(), x1.end(),
        inter.begin(), inter.end(),
        x1.begin()
      ),
      x1.end()
    );

    x2.erase(
      set_difference(
        x2.begin(), x2.end(),
        inter.begin(), inter.end(),
        x2.begin()
      ),
      x2.end()
    );
    if ( x1.size() == 0 && x2.size() == 0 ) return "Vqq";
    if ( x1.size() != 1 | x2.size() != 1 ) THROW(fatal_error,"Flavour changing of more or less than one quarks not implemented");
    if ( x1[0] == 2 | x2[0] == 2 ){ // u
      if ( x1[0] == 1 | x2[0] == 1 ) return "Vud";
      if ( x1[0] == 3 | x2[0] == 3 ) return "Vus";
      if ( x1[0] == 5 | x2[0] == 5 ) return "Vub";
    } else if ( x1[0] == 4 | x2[0] == 4 ){ // c
      if ( x1[0] == 1 | x2[0] == 1 ) return "Vcd";
      if ( x1[0] == 3 | x2[0] == 3 ) return "Vcs";
      if ( x1[0] == 5 | x2[0] == 5 ) return "Vcb";
    } else if ( x1[0] == 6 | x2[0] == 6 ){ // c
      if ( x1[0] == 1 | x2[0] == 1 ) return "Vtd";
      if ( x1[0] == 3 | x2[0] == 3 ) return "Vts";
      if ( x1[0] == 5 | x2[0] == 5 ) return "Vtb";
    }
    return "Vuu_Vdd";
  }
}

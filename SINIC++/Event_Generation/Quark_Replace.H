#ifndef SINIC_Main_Quark_Replace_H
#define SINIC_Main_Quark_Replace_H

#include "SINIC++/Event_Generation/Ladder.H"
#include "SINIC++/Event_Generation/Simple_MEs.H"
#include "MODEL/Main/Strong_Coupling.H"
#include "ATOOLS/Math/Poincare.H"
#include <set>

#include "ATOOLS/Math/Histogram.H"
#include <map>

namespace SINIC {

  class Quark_Replace {
  private:
    Simple_MEs               m_me2s;
    MODEL::Strong_Coupling * p_as;

    double GetTripletWeight(LadderMap::iterator yiter1, const TPropList::iterator piter,
                  int & flav,bool & glgl);
    void   ConstructMandelstams(const ATOOLS::Vec4D & out1,const ATOOLS::Vec4D & out2,
		     const double & mu2,double & hats, double & hatt, double & hatu);
    bool   MassiveKinematics(const double m,
                  ATOOLS::Vec4D & outmom1, ATOOLS::Vec4D & outmom2, double & y1, double & y2,
                  const LadderMap::iterator yiter1, const LadderMap::iterator yfirst, 
                  const LadderMap::iterator ylast);

  public:
    Quark_Replace();
    ~Quark_Replace();

    void   ReplaceWithQuarks(Ladder * ladder);
  };
}
#endif

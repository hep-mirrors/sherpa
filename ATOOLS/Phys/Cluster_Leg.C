#include "ATOOLS/Phys/Cluster_Leg.H"

#include "ATOOLS/Org/STL_Tools.H"
#include "ATOOLS/Org/MyStrStream.H"
#include <iomanip>

namespace ATOOLS {

  std::ostream &operator<<(std::ostream &ostr,const ColorID &col)
  {
    return ostr<<'('<<col.m_i<<','<<col.m_j<<')';
  }

  std::ostream &operator<<(std::ostream &ostr,const Cluster_Leg &leg)
  {
    ostr<<std::right<<std::setw(12)<<ToString(ID(leg.Id()))
	<<std::setw(12)<<leg.Flav()
	<<" "<<std::left<<leg.Mom()<<" "<<leg.Col();
    ostr<<" ["<<leg.Stat()<<","<<leg.NMax()<<"]";
    if (leg.K()>0) ostr<<" <-> "<<ID(leg.K());
    if (leg.Q2Shower()>=0.0) ostr<<" "<<sqrt(dabs(leg.Q2Shower()));
    return ostr;
  }

}

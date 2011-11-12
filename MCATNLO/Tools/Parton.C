#include "MCATNLO/Tools/Parton.H"
#include "ATOOLS/Phys/Cluster_Leg.H"

using namespace MCATNLO;
using namespace std;

namespace MCATNLO {
  std::ostream& operator<<(std::ostream& str, const Parton &part) {
    str<<"  "<<(part.m_pst==pst::IS?"IS":"FS")<<" Parton "
       <<&part<<" ("<<part.m_kin<<")["<<ATOOLS::ID(part.m_id)
       <<"]: "<<part.m_flav<<" : "<<part.m_mom
       <<" ("<<part.GetFlow(1)<<","<<part.GetFlow(2)<<")"
       <<"["<<part.GetRFlow(1)<<","<<part.GetRFlow(2)<<"]"<<endl;
    str<<"  k_T start : "<<sqrt(part.m_kt_start);
    str<<"  k_T test : "<<sqrt(part.m_kt_test);
    str<<"  k_T veto : "<<sqrt(part.m_kt_veto)<<"("<<sqrt(part.m_kt_max)<<")";
    str<<"  x_B : "<<part.m_xBj<<std::endl;
    return str;
  }
}

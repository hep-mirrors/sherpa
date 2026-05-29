#include "METOOLS/Explicit/Form_Factor.H"
#include "METOOLS/Explicit/Current.H"
#include "METOOLS/Explicit/Vertex.H"

#include "MODEL/Main/Single_Vertex.H"

#include "ATOOLS/Org/Message.H"
#include "ATOOLS/Org/Exception.H"
#include <set>

namespace METOOLS {
  class FFVMD: public Form_Factor {
  private:
    int    m_pos;
  public:
    FFVMD(const Vertex_Key &key);
    Complex FF();
  };// end of class FFPoint
}// end of namespace METOOLS


using namespace METOOLS;
using namespace ATOOLS;
using namespace std;
/////////////////////////////////////////////////////////////////////
//
// VMD form factors, basically assuming some (interfering)
// vector-meson propagators.
//
/////////////////////////////////////////////////////////////////////
FFVMD::FFVMD(const Vertex_Key &key):
  Form_Factor("VMD",key),
  m_pos(-1) {
  // Find the position of the (elementary) particle:
  // Normally this would be the photon, but I'll assume it is anything
  // that is not a hadron.
  // Also: extract the flavours in the vertex to look up parameters
  // in the model look-up table.
  std::set<kf_code> kfs;
  for (size_t i(0);i<key.m_j.size();++i) {
    if (!key.m_j[i]->Flav().IsHadron()) m_pos=i;
    kfs.insert(key.m_j[i]->Flav().Kfcode());
  }
  msg_Out()<<METHOD<<":\n";
}

Complex FFVMD::FF() {
  return Complex(1.,0.);
}

DECLARE_GETTER(FFVMD,"VMD",Form_Factor,Vertex_Key);
Form_Factor *Getter<Form_Factor,Vertex_Key,FFVMD>::
operator()(const Vertex_Key &args) const {
  msg_Out()<<METHOD<<"\n";
  return new FFVMD(args);
}

void ATOOLS::Getter<Form_Factor,Vertex_Key,FFVMD>::
PrintInfo(std::ostream &str,const size_t width) const { str<<"VMD"; }


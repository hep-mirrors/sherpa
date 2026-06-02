#include "METOOLS/Explicit/Form_Factor.H"
#include "METOOLS/Explicit/Current.H"
#include "METOOLS/Explicit/Vertex.H"
#include "METOOLS/FormFactors/Propagator.H"
#include "METOOLS/FormFactors/Line_Shapes.H"

#include "MODEL/Main/Single_Vertex.H"

#include "ATOOLS/Org/Message.H"
#include "ATOOLS/Org/Exception.H"
#include <set>

namespace METOOLS {
  enum FF_0_PP_mode {
    pipi_plus       = 1,
    KK_plus         = 2,
    Kpi_plus        = 11,
    pipi_0          = 101,
    KK_0            = 102,
    unknown         = 999
  };

  class FFVMD: public Form_Factor {
  private:
    METOOLS::Summed_Propagator m_props;
    FF_0_PP_mode m_mode;
    int          m_pos;
    void FixMode(const Vertex_Key &key);
    void Construct();
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
  m_props(Summed_Propagator()), m_mode(FF_0_PP_mode::unknown),
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
  FixMode(key);
  Construct();
  msg_Out()<<METHOD<<":\n";
}

void FFVMD::FixMode(const Vertex_Key &key) {
  if (key.m_j[0]->Flav()==Flavour(kf_pi_plus) &&
      key.m_j[1]->Flav()==Flavour(kf_pi_plus).Bar() )     m_mode = FF_0_PP_mode::pipi_0;
  else if (key.m_j[0]->Flav()==Flavour(kf_K_plus) &&
	   key.m_j[1]->Flav()==Flavour(kf_K_plus).Bar() ) m_mode = FF_0_PP_mode::KK_0;
}

void FFVMD::Construct() {
  switch (int(m_mode)) {
  case int(FF_0_PP_mode::pipi_0):
    m_props.Add(new METOOLS::BreitWigner(METOOLS::LineShapes->Get(Flavour(kf_rho_770))),
		Complex(1.,0.));
    m_props.Add(new METOOLS::BreitWigner(METOOLS::LineShapes->Get(Flavour(kf_omega_782))),
		Complex(0.0,0.1));
    break;
  case int(FF_0_PP_mode::unknown):
  default:
    THROW(fatal_error,"Cannot construct form factor for unknown FF_0_PP mode.");
    break;
  }
}

Complex FFVMD::FF() {
  Current *j = m_pos<0?p_v->JC():p_v->J(m_pos);
  // there was a minus sign before in Q2.
  double Q2  = j->P().Abs2();
  /*
  for (map<Propagator_Base *,Complex>::iterator pit=m_props.GetAll().begin();
       pit!=m_props.GetAll().end();pit++) {
    msg_Out()<<METHOD<<"("<<Q2<<"): "<<pit->second<<" * "
	     <<(*pit->first)(Q2)<<".\n";
  }
  msg_Out()<<METHOD<<": F(Q^2 = "<<Q2<<") = "<<m_props(Q2)<<"\n";
  */
  return m_props(Q2);
}



DECLARE_GETTER(FFVMD,"VMD",Form_Factor,Vertex_Key);
Form_Factor *Getter<Form_Factor,Vertex_Key,FFVMD>::
operator()(const Vertex_Key &args) const {
  msg_Out()<<METHOD<<"\n";
  return new FFVMD(args);
}

void ATOOLS::Getter<Form_Factor,Vertex_Key,FFVMD>::
PrintInfo(std::ostream &str,const size_t width) const { str<<"VMD"; }


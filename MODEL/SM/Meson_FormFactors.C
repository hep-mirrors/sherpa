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
    void ConstructPionFormFactor();
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
  msg_Out()<<METHOD<<"("<<key.m_j.size()<<"): "<<key.m_j[0]->Flav()<<", "
	   <<key.m_j[1]->Flav()<<"\n";
  msg_Out()<<key.p_mv->in[0]<<", "<<key.p_mv->in[1]<<", "<<key.p_mv->in[2]<<"\n";
  for (size_t i(0);i<key.m_j.size();++i) {
    if (!key.m_j[i]->Flav().IsHadron()) m_pos=i;
    else kfs.insert(key.m_j[i]->Flav().Kfcode());
  }
  msg_Out()<<METHOD<<"("<<key.m_j.size()<<"): "<<key.m_j[0]->Flav()<<", "
	   <<key.m_j[1]->Flav()<<"\n"; //<<", "<<key.m_j[2]->Flav()
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
    ConstructPionFormFactor();
    break;
  case int(FF_0_PP_mode::unknown):
  default:
    msg_Out()<<METHOD<<" yields no form factor.\n";
    // THROW(fatal_error,"Cannot construct form factor for unknown FF_0_PP mode.");
    break;
  }
}

void FFVMD::ConstructPionFormFactor() {
  METOOLS::Summed_Propagator * rhofac = new METOOLS::Summed_Propagator();
  rhofac->Add(new METOOLS::Unity(),Complex(1.,0.));
  rhofac->Add(new METOOLS::WeightedBreitWigner(LineShapes->Get(Flavour(kf_omega_782))),
	      0.00158*Complex(cos(0.075),sin(0.075)) );
  // I will need to make a phi(1020) line shape or, better, produce a simplistic
  // lineshape with a fixed width (which boils down to make a fixed-width propagator).
  METOOLS::Multiplied_Propagator * rho = new METOOLS::Multiplied_Propagator();
  rho->Add(new METOOLS::GounarisSakurai(LineShapes->Get(Flavour(kf_rho_770))));
  rho->Add(rhofac);
  m_props = METOOLS::Summed_Propagator();
  m_props.Add(rho);
  m_props.Add(new METOOLS::GounarisSakurai(LineShapes->Get(Flavour(kf_rho_1450))),
	      0.14104*Complex(cos(3.7797),sin(3.7797)) );
  m_props.Add(new METOOLS::GounarisSakurai(LineShapes->Get(Flavour(kf_rho_1700))),
	      0.0614*Complex(cos(1.429),sin(1.429)) );
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
  msg_Out()<<METHOD<<"(size = "<<args.m_j.size()<<")\n";
  return new FFVMD(args);
}

void ATOOLS::Getter<Form_Factor,Vertex_Key,FFVMD>::
PrintInfo(std::ostream &str,const size_t width) const { str<<"VMD"; }


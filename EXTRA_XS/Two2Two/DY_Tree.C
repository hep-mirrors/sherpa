#include "EXTRA_XS/Main/ME_Base.H"
#include "ATOOLS/Org/Exception.H"
#include "METOOLS/Main/XYZFuncs.H"
#include "MODEL/Main/Model_Base.H"
#include "MODEL/Interaction_Models/Interaction_Model_Base.H"
#include "ATOOLS/Org/Run_Parameter.H"

using namespace EXTRAXS;
using namespace PHASIC;
using namespace ATOOLS;
using namespace METOOLS;

namespace EXTRAXS {
  class DY_Tree : public ME_Base {

    bool   m_Z_on, m_P_on;
    double m_mz, m_wz;
    double m_sin2tw, m_cos2tw, m_eq1, m_eq2, m_y3f1, m_y3f2, m_v1, m_a1, m_v2, m_a2;
    double m_aqed, m_pref_qed, m_pref_Z, m_colfac;
    const Flavour_Vector m_flavvec;

  public:
    DY_Tree(const Process_Info& pi, const Flavour_Vector& flavs, const Particle_Vector& parts) :
      ME_Base(pi, flavs, parts),
  m_Z_on(ATOOLS::Flavour(kf_Z).IsOn()), m_P_on(ATOOLS::Flavour(kf_photon).IsOn()),
  m_mz(ATOOLS::Flavour(kf_Z).Mass()),
  m_wz(ATOOLS::Flavour(kf_Z).Width()),
  m_sin2tw(MODEL::s_model->ScalarConstant(std::string("sin2_thetaW"))),m_cos2tw(1.-m_sin2tw),
  m_eq1(flavs[0].Charge()),
  m_eq2(flavs[2].Charge()),
  m_y3f1((2.*int(flavs[0].IsUptype())-1)/2.),
  m_y3f2((2.*int(flavs[2].IsUptype())-1)/2.),
  m_v1(m_y3f1-2.*m_eq1*m_sin2tw), m_a1(m_y3f1),
  m_v2(m_y3f2-2.*m_eq2*m_sin2tw), m_a2(m_y3f2),
  m_aqed(MODEL::s_model->GetInteractionModel()->ScalarFunction("alpha_QED",sqr(rpa.gen.Ecms()))),
  m_pref_qed(4.*M_PI*m_aqed),m_pref_Z((4.*M_PI*m_aqed)/(4.*m_sin2tw*m_cos2tw)), 
  m_flavvec(flavs)

 {
  m_colfac = 1.;

  if (flavs[2].IsQuark() && flavs[0].IsLepton()) {
        m_colfac = 3.;
  }
  if (flavs[0].IsQuark() && flavs[2].IsLepton())  {
      m_colfac  = 1./3.;
  }

    }

    ~DY_Tree() {

    }

    void Calc(const ATOOLS::Vec4D_Vector& momenta);
  };
}


void DY_Tree::Calc(const Vec4D_Vector& momenta) {

  const Vec4D *p_momenta = &momenta[0];
  const Flavour *p_flavs = &m_flavvec[0];

  XYZFunc F(4,p_momenta,p_flavs,1);

  Complex amp(0.0,0.0);
  Complex i(0.0,1.0);

  Complex cR(1.0,0.0);
  Complex cL(1.0,0.0);

  Complex prop=-i/(momenta[0]+momenta[1]).Abs2();
  Complex factor=-i*m_pref_qed*m_eq1*m_eq2*prop*sqrt(m_colfac)/2.;

  Complex cR1 = (m_v1-m_a1);
  Complex cL1 = (m_v1+m_a1);

  Complex cR2 = (m_v2-m_a2);
  Complex cL2 = (m_v2+m_a2);

  Complex propz=-i/((momenta[0]+momenta[1]).Abs2()-sqr(m_mz)+i*m_mz*m_wz);
  Complex factorz=i*m_pref_Z*propz*sqrt(m_colfac)/2.;

  vector<int> spins(4);
  for (spins[0]=0; spins[0]<2; ++spins[0]) {
    for (spins[1]=0; spins[1]<2; ++spins[1]) {
      for (spins[2]=0; spins[2]<2; ++spins[2]) {
        for (spins[3]=0; spins[3]<2; ++spins[3]) {
	  amp = Complex(0.0,0.0);
	  if (m_P_on){
	    amp+=factor*F.Z(1,spins[1], 0,spins[0], 2,spins[2], 3,spins[3], cR,cL,cR,cL);
	  }
	  if (m_Z_on){
	    amp+=factorz*F.Z(1,spins[1], 0,spins[0], 2,spins[2], 3,spins[3], cR1,cL1,cR2,cL2);
	  }
          m_res.Insert(amp, spins);
        }
      }
    }
  }

}



DECLARE_ME_GETTER(DY_Tree_Getter,"DY_Tree")
ME_Base *DY_Tree_Getter::operator()(const Process_Info &pi) const
{
  return NULL;
  if (pi.m_fi.m_nloqcdtype==nlo_type::lo && pi.m_fi.m_nloewtype==nlo_type::lo) {
    Flavour_Vector fl=pi.ExtractFlavours();
    if (fl.size()!=4) return NULL;
    if ((fl[2].IsFermion() && fl[3]==fl[2].Bar() &&
         fl[0].IsFermion()  && fl[1]==fl[0].Bar())  &&
        (fl[0]!=fl[2] && fl[0]!=fl[3])) {
      if ((pi.m_oqcd==0 || pi.m_oqcd==99) && (pi.m_oew==2 || pi.m_oew==99)) {

   Particle_Vector particles;
   Vec4D p_i(0.,0.,0.,0.);
   for(size_t i=0;i<4;i++){
     Particle part(i, fl[i], p_i,'a');
     Particle *p_part;
     p_part = &part;
     particles.push_back(p_part);
   }
std::cout<<"DY_Tree_Getter"<<endl;
       return new DY_Tree(pi, fl, particles);
      }
    }
  }
  return NULL;
}


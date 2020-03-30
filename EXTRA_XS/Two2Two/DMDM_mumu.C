#include "PHASIC++/Process/External_ME_Args.H"
#include "MODEL/Main/Running_AlphaS.H"
#include "MODEL/Main/Running_AlphaQED.H"
#include "ATOOLS/Phys/Flow.H"
#include "MODEL/Main/Model_Base.H"
#include "MODEL/UFO/UFO_Model.H"
#include "EXTRA_XS/Main/ME2_Base.H"

#define PropID(i,j) ((1<<i)|(1<<j))

using namespace EXTRAXS;
using namespace MODEL;
using namespace ATOOLS;
using namespace PHASIC;
using namespace std;

namespace EXTRAXS {
  class DMDM_mumu : public ME2_Base {
    /* Describing the annihilation of fermionic dark matter with its antiparticle 
       through a Z to produce a muon-antimuon pair. 
    */
  private:
  public:
    DMDM_mumu(const External_ME_Args& args);
    double operator()(const Vec4D_Vector& mom);
    bool SetColours(const Vec4D_Vector& mom);
  };
}

DMDM_mumu::DMDM_mumu(const External_ME_Args& args) :
  ME2_Base(args)
{
  m_sintt=1; //what should this be? initialised to 0, but most processes use 1
  for (short int i=0;i<4;i++) m_colours[i][0] = m_colours[i][1] = 0;
  m_oew = 2; m_oqcd = 0;
  m_cfls[3] = new Flavour_Vector;
  m_cfls[3]->push_back(kf_photon);
  m_cfls[3]->push_back(kf_Z);
  m_cfls[12] = new Flavour_Vector;
  m_cfls[12]->push_back(kf_photon);
  m_cfls[12]->push_back(kf_Z);
}

double DMDM_mumu::operator()(const Vec4D_Vector& mom)
{
  double s=(mom[0]+mom[1]).Abs2();

  // how to set these in the model from here/runcard?
  double V(-0.1), A(0.1);
  double Vtil = sqrt(4*M_PI/128.) * (-0.5 + 2*0.23)/ (2*sqrt(0.23)*0.88); //lepton vector coupling constant
  double Atil = sqrt(4*M_PI/128.) *0.5 / (2*0.23*0.88); //lepton axial c.c.

  double M = ATOOLS::Flavour(kf_Z).Mass();
  double gamma = ATOOLS::Flavour(kf_Z).Width();
  double m_DM = m_flavs[0].Mass();

  // // positive definite definition
  // double cos_theta = (s/4 - mom_dot)/sqrt(s*s/16 - s*sqr(m_DM)/4);
  double mom_dot = 0;
  for (short int i = 1; i < 4; i++) mom_dot += mom[0][i]*mom[2][i]; // just spatial parts

  double cos_theta = mom_dot/(Vec3<double>(mom[0]).Abs()*Vec3<double>(mom[2]).Abs()); // just spatial parts 
  // positive definite definition
  if (cos_theta < 0) cos_theta = -cos_theta;
  if (cos_theta > 1.) {
    msg_Error() << "Cosine greater than 1!" << endl;
  }

  double factor1 = 4/(sqr(s-M*M) + M*M*sqr(gamma));
  double part1 = (V*V+A*A)*(Vtil*Vtil+Atil*Atil)*(s*s/8 + s/2 * (s/4 - sqr(m_DM)) * sqr(cos_theta));
  double part2 = 2*V*A*Vtil*Atil * s* sqrt(s*s/4 - s*sqr(m_DM)) * cos_theta;
  double part3 = (V*V-A*A)*(Vtil*Vtil+Atil*Atil)*s*sqr(m_DM);

  // cout << 1/(16*M_PI * 2*s*sqrt(1-4*sqr(m_DM)/s)) * factor1*((V*V+A*A)*(Vtil*Vtil+Atil*Atil)*(sqr(s)/3-sqr(m_DM)*s/3) + 2*part3) << endl; //debugging

  //cout << factor1*(part1+part2+part3)/m_symfac << endl; //debugging
  return factor1*(part1+part2+part3)/m_symfac;
  // return 2*pow(4*M_PI,4)*factor1*(part1+part2+part3)/m_symfac; //test
}


bool DMDM_mumu::SetColours(const Vec4D_Vector& mom)
{
  return true;
}

DECLARE_TREEME2_GETTER(DMDM_mumu,"DMDM_mumu")
Tree_ME2_Base *ATOOLS::Getter<Tree_ME2_Base,External_ME_Args,DMDM_mumu>::
operator()(const External_ME_Args &args) const
{
  msg_Out()<<METHOD<<".\n";
  if (dynamic_cast<UFO::UFO_Model*>(MODEL::s_model)) return NULL;
  if (MODEL::s_model->Name()!="SMDM") return NULL;

  const Flavour_Vector fl=args.Flavours();
  if (fl.size()!= 4) return NULL;
  if (fl[0].Kfcode()==52 && fl[1].Kfcode()==52 &&
      fl[2].IsFermion() && fl[2].Charge() &&
      fl[3]==fl[2].Bar()) {
        if (args.m_orders[0]==0 && args.m_orders[1]==2) {
	  std::cout<<"   initialising ME.\n";
          return new DMDM_mumu(args);
        }
  }
  return NULL;
}

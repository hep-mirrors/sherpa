#include "EXTRA_XS/Main/ME_Base.H"
#include "Exception.H"
#include "HELICITIES/Main/XYZFuncs.H"
#include "MODEL/Main/Model_Base.H"
#include "ATOOLS/Org/Run_Parameter.H"
#include "ATOOLS/Phys/Flow.H"
#include "MODEL/Main/Running_AlphaQED.H"

using namespace EXTRAXS;
using namespace PHASIC;
using namespace ATOOLS;
using namespace HELICITIES;
using namespace std;


namespace EXTRAXS {
  class nunuee_EW_Real_ME : public ME_Base {

    bool   m_Z_on, m_P_on;
    double m_mz2, m_wz2;
    double m_sin2tw, m_cos2tw, m_eq1, m_eq2, m_y3f1, m_y3f2, m_v1, m_a1, m_v2, m_a2;
    double m_aqed, m_pref_qed, m_pref_Z, colfac;
    const Flavour_Vector m_flavvec;

  public:
    nunuee_EW_Real_ME(const Process_Info& pi, const Flavour_Vector& flavours, const Particle_Vector& parts) :
      ME_Base(pi, flavours, parts),
  m_Z_on(ATOOLS::Flavour(kf_Z).IsOn()), m_P_on(ATOOLS::Flavour(kf_photon).IsOn()),
  m_mz2(ATOOLS::sqr(ATOOLS::Flavour(kf_Z).Mass())),
  m_wz2(ATOOLS::sqr(ATOOLS::Flavour(kf_Z).Width())),
  m_sin2tw(MODEL::s_model->ScalarConstant(std::string("sin2_thetaW"))),m_cos2tw(1.-m_sin2tw),
  m_eq1(flavours[0].Charge()),
  m_eq2(flavours[3].Charge()),
  m_y3f1((2.*int(flavours[0].IsUptype())-1)/2.),
  m_y3f2((2.*int(flavours[3].IsUptype())-1)/2.),
  m_v1(m_y3f1-2.*m_eq1*m_sin2tw), m_a1(m_y3f1),
  m_v2(m_y3f2-2.*m_eq2*m_sin2tw), m_a2(m_y3f2),
//   m_aqed(MODEL::aqed->Aqed(m_mz2)),
  m_aqed(MODEL::aqed->Aqed((ATOOLS::sqr(ATOOLS::rpa.gen.Ecms())))),
  m_pref_qed(4.*M_PI*m_aqed),m_pref_Z((4.*M_PI*m_aqed)/(4.*m_sin2tw*m_cos2tw)),
  m_flavvec(pi.ExtractFlavours())

    {

if(flavours[3].IsQuark()) colfac = sqrt(3.);
else colfac = 1.;
//   ATOOLS::CMatrix* colour;
//   colour = new ATOOLS::CMatrix(1);
// 
//   (*colour)[0][0] = 4./3.;
// 
//  m_colres.SetColorMatrix(colour);


    }

    ~nunuee_EW_Real_ME() {
    }

    void Calc(const ATOOLS::Vec4D_Vector& momenta);

  };
}

void nunuee_EW_Real_ME::Calc(const Vec4D_Vector& momenta) {

  const Vec4D *p_momenta = &momenta[0];
  const Flavour *p_flavs = &m_flavvec[0];
  XYZFunc F(5,p_momenta,p_flavs,1);

  Complex amp(0.,0.);
  Complex i(0.0,1.0);

  Complex cR(1.0,0.0);
  Complex cL(1.0,0.0);

  Complex cL1 = (m_v1+m_a1);
  Complex cR1 = (m_v1-m_a1);
  Complex cL2 = (m_v2+m_a2);
  Complex cR2 = (m_v2-m_a2);

  double cpl = sqrt(m_pref_qed);

  Complex prop=-1./(momenta[0]+momenta[1]).Abs2();
  Complex factor1=m_pref_qed*m_eq1*m_eq2*prop*cpl*m_eq2*colfac/2./(momenta[3]+momenta[2]).Abs2();
  Complex factor2=m_pref_qed*m_eq1*m_eq2*prop*cpl*m_eq2*colfac/2./(momenta[4]+momenta[2]).Abs2();

  Complex propz=-1./((momenta[0]+momenta[1]).Abs2()-m_mz2+i*sqrt(m_mz2*m_wz2));
  Complex factorz1=m_pref_Z*propz*cpl*m_eq2*colfac/2./(momenta[3]+momenta[2]).Abs2();
  Complex factorz2=m_pref_Z*propz*cpl*m_eq2*colfac/2./(momenta[4]+momenta[2]).Abs2();

  //ppol: photon polarisation vector, TG diploma eq(3.29)
  ATOOLS:: Vec4C   ppol;
  vector<pair<int,int> > spins(5);

  for (size_t hf1=0; hf1<2; ++hf1) { // e-
    spins[0] = make_pair(0,hf1);
    for (size_t hf1b=0; hf1b<2; ++hf1b) { // e+
      spins[1] = make_pair(1,hf1b);
      for (size_t hg=0; hg<2; ++hg) { //photon
        spins[2] = make_pair(2,hg);
        ppol[1] = (momenta[2].CosTheta()*momenta[2].CosPhi()+(1.-2.*hg)*i*momenta[2].SinPhi())/sqrt(2.);
        ppol[2] = (momenta[2].CosTheta()*momenta[2].SinPhi()-(1.-2.*hg)*i*momenta[2].CosPhi())/sqrt(2.);
        ppol[3] = -(momenta[2].SinTheta())/sqrt(2.);
        for (size_t hf2=0; hf2<2; ++hf2) { // q
          spins[3] = make_pair(3,hf2);
            for (size_t hf2b=0; hf2b<2; ++hf2b) { // qb
              spins[4] = make_pair(4,hf2b);
              amp= Complex(0.0,0.0);
             if (hf1==hf1b && hf2==hf2b){
              for (size_t h=0; h<2; ++h) {
                if(m_P_on){
                amp+=factor2*(F.Z(1,hf1b,0,hf1, 3,hf2, 2,h, cR,cL,cR,cL)*F.X(2,h, ppol, 4,hf2b, cR,cL) + (momenta[4]*ppol)*F.Z(1,hf1b,0,hf1, 3,hf2, 4,hf2b, cR,cL,cR,cL));
                amp+=-factor1*((momenta[3]*ppol)*F.Z(1,hf1b,0,hf1, 3,hf2, 4,hf2b, cR,cL,cR,cL) + F.Z(1,hf1b,0,hf1, 2,h, 4,hf2b, cR,cL,cR,cL) * F.X(3, hf2, ppol, 2, h, cR,cL));
                }
                if(m_Z_on){
                amp+=factorz2*(F.Z(1,hf1b,0,hf1, 3,hf2, 2,h, cR1,cL1,cR2,cL2)*F.X(2,h, ppol, 4,hf2b, cR,cL) + (momenta[4]*ppol)*F.Z(1,hf1b,0,hf1, 3,hf2, 4,hf2b, cR1,cL1,cR2,cL2));
                amp+=-factorz1*((momenta[3]*ppol)*F.Z(1,hf1b,0,hf1, 3,hf2, 4,hf2b, cR1,cL1,cR2,cL2)+ F.Z(1,hf1b,0,hf1, 2,h, 4,hf2b, cR1,cL1,cR2,cL2)*F.X(3, hf2, ppol, 2, h, cR,cL));
                }
              }
             }
// std::cout<<" spins   =  "<<hf1<<hf1b<<hg<<hf2<<hf2b<<endl;
// std::cout<<"   norm amp    =  "<<norm(amp)<<endl;
          m_res.Insert(amp, spins);
          }
        }
      }
    }
  }
// exit(0);
  return; 
}


DECLARE_ME_GETTER(nunuee_EW_Real_ME_Getter,"nunuee_EW_Real_ME")
ME_Base *nunuee_EW_Real_ME_Getter::operator()(const Process_Info &pi) const
{
  if (pi.m_fi.m_nloewtype==nlo_type::real && pi.m_fi.m_nloqcdtype==nlo_type::lo
      || pi.m_fi.m_nloqcdtype==nlo_type::lo && pi.m_fi.m_nloewtype==nlo_type::lo) {
    Flavour_Vector fl=pi.ExtractFlavours();
    if (fl.size() != 5) return 0.;
    if (fl[0].IsLepton() && fl[1]==fl[0].Bar() &&
        fl[2].IsPhoton() && fl[3].IsFermion() && fl[4]==fl[3].Bar()) {
      if ((pi.m_oqcd==0 || pi.m_oqcd==99) && (pi.m_oew==2 || pi.m_oew==99)) {
        Particle_Vector particles;
        Vec4D p_i(0.,0.,0.,0.);
        for(size_t i=0;i<5;i++){
          Particle part(i, fl[i], p_i,'a');
          Particle *p_part;
          p_part = &part;
          particles.push_back(p_part);
        }
        return new nunuee_EW_Real_ME(pi, fl, particles);
      }
    }
  }
  return NULL;
}



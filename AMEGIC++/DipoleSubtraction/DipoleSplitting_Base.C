#include "AMEGIC++/DipoleSubtraction/DipoleSplitting_Base.H"
#include "AMEGIC++/Main/ColorSC.H"
#include "MODEL/Main/Model_Base.H"
#include "ATOOLS/Org/Run_Parameter.H"
#include "ATOOLS/Org/Exception.H"
#include "ATOOLS/Org/Message.H"
#include "ATOOLS/Org/Shell_Tools.H"
#include "ATOOLS/Org/MyStrStream.H"
#include "ATOOLS/Org/Scoped_Settings.H"

#include <typeinfo>

using namespace ATOOLS;
using namespace AMEGIC;
using namespace MODEL;
using namespace std;

#define SQRT_05 0.70710678118654757

std::ostream& operator<<(std::ostream& str, const DipoleCase& dc)
{
  if     (dc == DipoleCase::CS)  return str << "DipoleCase = CS\n";
  else if(dc == DipoleCase::IDa)  return str << "DipoleCase = IDa\n";
  else if(dc == DipoleCase::IDb)  return str << "DipoleCase = IDb\n";
  else if(dc == DipoleCase::IDin) return str << "DipoleCase = IDin\n";
  else                            return str << "DipoleCase = non\n";
  return str;
}

DipoleSplitting_Base::DipoleSplitting_Base(sbt::subtype st,
                                           spt::splittingtype ft,
                                           dpt::dipoletype dt,
                                           int m,int i,int j,int k) :
  m_name(""), m_alpha(1.), m_kt2max(std::numeric_limits<double>::max()),
  m_amin(max(ATOOLS::Accu(),1.e-8)), m_kappa(2./3.),
  m_Q2(0.), m_kt2(0.), m_a(-1.),
  m_sff(std::nan), m_av(std::nan), m_fac(1.),
  m_mcmode(-1), m_mcsign(-1),
  m_stype(st), m_dtype(dt), m_ftype(ft),
  m_i(i), m_j(j), m_k(k), m_tij(-1), m_tk(-1),
  m_m(m), m_es(-1), m_k0sqf(-1.), m_k0sqi(-1.),
  m_spfdef(0.), m_collVFF(true), m_Vsubmode(1),
  p_cpl(NULL)
{
  p_nlomc=NULL;
  Settings& s = Settings::GetMainSettings();
  m_subtype = s["NLO_SUBTRACTION_SCHEME"].Get<int>();
  m_name=ToString(m_dtype)+"["+ToString(m_ftype)
         +"("+ToString(m_i)+","+ToString(m_j)+")("+ToString(m_k)+")]"
         +"("+ToString(m_stype)+")";
  DEBUG_FUNC(m_name);

  m_amin=s["DIPOLES"]["AMIN"].Get<double>();
  m_kappa=s["DIPOLES"]["KAPPA"].Get<double>();

  if (m_stype==sbt::qed) m_collVFF=false;
  m_collVFF=s["DIPOLES"]["COLLINEAR_VFF_SPLITTINGS"].Get<int>();
  m_Vsubmode=s["DIPOLES"]["V_SUBTRACTION_MODE"].Get<int>();
  if (m_Vsubmode<0 || m_Vsubmode>2) THROW(fatal_error,"Unknown mode.");
  std::string vsm("");
  if      (m_Vsubmode==0) vsm="scalar";
  else if (m_Vsubmode==1) vsm="fermionic";
  else if (m_Vsubmode==2) vsm="eikonal";
  msg_Tracking()<<"Use "<<vsm<<" V->VP subtraction term."<<std::endl;

  if (st==sbt::qcd) {
    if      (m_ftype==spt::g2qq) m_fac = CSC.TR/CSC.CA;
    else if (m_ftype==spt::g2gg) m_fac = 2.;
  }

  m_pfactors.clear(); m_pfactors.push_back(1.);
  m_k0sqf = s["CSS_FS_PT2MIN"].Get<double>();
  m_k0sqi = s["CSS_IS_PT2MIN"].Get<double>();
  m_es = s["CSS_EVOLUTION_SCHEME"].Get<int>();
  if (m_subtype==1) m_kappa=1.0;

  std::string dipole_string = s["DIPOLES"]["CASE"].Get<std::string>();
  if      (dipole_string == "CS")   m_dipolecase=DipoleCase::CS;
  else if (dipole_string == "IDa")  m_dipolecase=DipoleCase::IDa;
  else if (dipole_string == "IDb")  m_dipolecase=DipoleCase::IDb;
  else if (dipole_string == "IDin") m_dipolecase=DipoleCase::IDin;
  else                              m_dipolecase=DipoleCase::CS;
  switch(m_dipolecase){
    case IDa: m_pIDspec.resize(1); break; // # of kinematic spectators
    case IDb: m_pIDspec.resize(2); break; // # of kinematic spectators
    default: break;
  }
}

void DipoleSplitting_Base::SetCoupling(const MODEL::Coupling_Map *cpls)
{
  std::string cplname("");
  if      (m_stype==sbt::qcd) cplname="Alpha_QCD";
  else if (m_stype==sbt::qed) cplname="Alpha_QED";
  else THROW(fatal_error,"Cannot set coupling for subtraction type"
                         +ToString(m_stype));
  msg_Debugging()<<Name()<<" : "<<cplname<<std::endl;

  if (cpls->find(cplname)!=cpls->end()) p_cpl=cpls->find(cplname)->second;
  else THROW(fatal_error,"Coupling not found");
  msg_Tracking()<<METHOD<<"(): "<<cplname<<" = "<<*p_cpl<<std::endl;
  m_spfdef = -8.*M_PI*p_cpl->Default();
}

void DipoleSplitting_Base::CalcVectors(Vec4D& p1, Vec4D& p2, double B)
{
  m_dpollist.clear();
  m_pfactors.clear();

  Vec3D pv(p2);
  Vec3D ptp=Vec3D(p1)-(p1[0]/p2[0])*pv;
  Vec3D ptt=cross(ptp,pv);

  m_dpollist.push_back(Vec4D(0.,ptt/ptt.Abs()));
  m_pfactors.push_back(1.);

  Vec4D vh(0.,ptp/ptp.Abs());
  m_dpollist.push_back(vh);
  m_pfactors.push_back((B-1.)/B);
}

double DipoleSplitting_Base::GetR(const Vec4D* mom,const Vec4D* LOmom)
{
  double ptijk=2.*(m_ptij*m_pj)*(m_ptk*m_pj)/(m_ptij*m_ptk);
  double spt=0.;
  for(int a=2;a<m_m;a++) {
    for(int b=a+1;b<m_m;b++) {
      for(int c=2;c<m_m;c++)
	if(c!=a&&c!=b) spt+=0.5*sqr((LOmom[a]+LOmom[b])*LOmom[c])/
			    ((LOmom[a]*LOmom[b])*(LOmom[b]*LOmom[c])
			     *(LOmom[c]*LOmom[a]));
    }
  }
  return 1./(1+ptijk*spt);
}

bool DipoleSplitting_Base::Reject(const double &alpha)
{
  if (IsBad(m_av))
    msg_Error()<<METHOD<<"(): Average is "<<m_av<<" in "
	       <<Demangle(typeid(*this).name())
	       <<"[type="<<m_ftype<<"]"<<std::endl;
  if (m_mcmode==1) {
    // discard m_av>0 check, to improve integration in resonance-aware matching
    int da(/*m_av>0.0 &&*/ (m_kt2<m_kt2max || IsEqual(m_kt2,m_kt2max,1.0e-6))),
        ds(alpha<=m_alpha);
    msg_Debugging()<<"kt = "<<sqrt(m_kt2)<<", ktmax = "<<sqrt(m_kt2max)
		   <<" -> DA = "<<da<<", DS = "<<ds<<" -> DA-DS = "<<da-ds<<"\n";
    m_mcsign=ds-da;
    return m_mcsign==0;
  }
  if (m_mcmode==2) {
    // discard m_av>0 check, to improve integration in resonance-aware matching
    m_mcsign=/*m_av>0.0 &&*/ (m_kt2<m_kt2max || IsEqual(m_kt2,m_kt2max,1.0e-6));
    msg_Debugging()<<"kt = "<<sqrt(m_kt2)<<", ktmax = "<<sqrt(m_kt2max)
		   <<" -> DA = "<<m_mcsign<<"\n";
    return m_mcsign==0;
  }
  return alpha>m_alpha || m_kt2>m_kt2max;
}

double DipoleSplitting_Base::GetF()
{
  DEBUG_FUNC("a="<<m_a<<", alpha="<<m_alpha<<", sf="<<m_sff<<", av="<<m_av);
  if (Reject(m_a)) return 0.;
  else return GetValue();
}

double DipoleSplitting_Base::GetValue()
{
  THROW(fatal_error, "Virtual function not reimplemented.");
  return 0.0;
}

bool DipoleSplitting_Base::KinCheck() const
{
  if (m_amin>0.0) return m_a>m_amin;
  return m_kt2>-m_amin;
}

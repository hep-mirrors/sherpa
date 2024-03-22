#include "PHASIC++/Process/ME_Generator_Base.H"

#include "PHASIC++/Process/ME_Generators.H"
#include "ATOOLS/Org/Run_Parameter.H"
#include "ATOOLS/Org/Scoped_Settings.H"
#include "ATOOLS/Math/Function_Base.H"
#include "ATOOLS/Math/Poincare.H"
#include "ATOOLS/Math/BreitBoost.H"
#include "ATOOLS/Phys/Flavour.H"
#include "MODEL/Main/Model_Base.H"
#include "REMNANTS/Main/Remnant_Handler.H"
#include <algorithm>

#define COMPILE__Getter_Function
#define OBJECT_TYPE PHASIC::ME_Generator_Base
#define PARAMETER_TYPE PHASIC::ME_Generator_Key
#define EXACTMATCH false
#include "ATOOLS/Org/Getter_Function.C"

using namespace PHASIC;
using namespace ATOOLS;

ME_Generator_Base::ME_Generator_Base(const std::string &name):
  m_name(name), m_massmode(0), p_gens(NULL), p_remnant(NULL)
{
  RegisterDefaults();
}

ME_Generator_Base::~ME_Generator_Base()
{
}

void ME_Generator_Base::RegisterDefaults()
{
  RegisterDipoleParameters();
  RegisterNLOParameters();
}

Process_Base *ME_Generator_Base::InitializeProcess
(Cluster_Amplitude *const ampl,const int mode,
 const std::string &gen,const std::string &addname)
{
  Process_Info pi;
  pi.m_addname=addname;
  pi.m_megenerator=gen.length()?gen:m_name;
  for (size_t i(0);i<ampl->NIn();++i) {
    Flavour fl(ampl->Leg(i)->Flav().Bar());
    if (Flavour(kf_jet).Includes(fl)) fl=Flavour(kf_jet);
    pi.m_ii.m_ps.push_back(Subprocess_Info(fl,"",""));
  }
  for (size_t i(ampl->NIn());i<ampl->Legs().size();++i) {
    Flavour fl(ampl->Leg(i)->Flav());
    if (Flavour(kf_jet).Includes(fl)) fl=Flavour(kf_jet);
    pi.m_fi.m_ps.push_back(Subprocess_Info(fl,"",""));
  }
  if (mode&8) {
    pi.m_maxcpl[0]=pi.m_mincpl[0]=ampl->OrderQCD();
    pi.m_maxcpl[1]=pi.m_mincpl[1]=ampl->OrderEW();
  }
  PHASIC::Process_Base *proc=p_gens->InitializeProcess(pi,mode&1);
  if (proc==NULL) return proc;
  proc->SetSelector(Selector_Key{});
  std::string stag("VAR{"+ToString(sqr(rpa->gen.Ecms()))+"}");
  proc->SetScale(Scale_Setter_Arguments(MODEL::s_model,stag,"Alpha_QCD 1"));
  proc->SetKFactor(KFactor_Setter_Arguments("None"));
  proc->PerformTests();
  return proc;
}

void ME_Generator_Base::SetPSMasses()
{
  Settings& s = Settings::GetMainSettings();
  ATOOLS::Flavour_Vector allflavs(MODEL::s_model->IncludedFlavours());
  std::vector<size_t> defpsmassive,defpsmassless;
  const std::vector<size_t> psmassive  { s["MASSIVE_PS"].GetVector<size_t>()  };
  const std::vector<size_t> psmassless { s["MASSLESS_PS"].GetVector<size_t>() };
  const bool respect{ s["RESPECT_MASSIVE_FLAG"].Get<bool>() };
  // check consistency
  for (size_t i(0);i<psmassive.size();++i)
    if (std::find(psmassless.begin(),psmassless.end(),psmassive[i])!=
        psmassless.end()) THROW(fatal_error,"Inconsistent input.");
  for (size_t i(0);i<psmassless.size();++i)
    if (Flavour(psmassless[i]).IsMassive())
      THROW(fatal_error,"Cannot shower massive particle massless.");
  // set defaults
  // respect=false -> def: dusgy massless, rest massive
  // respect=true  -> def: only massive massive, rest massless
  // TODO: need to fill in those that are massive already?
  if (!respect) {
    defpsmassless.push_back(kf_d);
    defpsmassless.push_back(kf_u);
    defpsmassless.push_back(kf_s);
    defpsmassless.push_back(kf_gluon);
    defpsmassless.push_back(kf_photon);
    // in a DIS-like setup, always respect massive flag of lepton
    if(rpa->gen.Beam2().IsLepton() && !rpa->gen.Beam2().IsMassive()
       && rpa->gen.Beam1().IsHadron()) {
        defpsmassless.push_back(rpa->gen.Beam2().Kfcode());
    }
    else if(rpa->gen.Beam1().IsLepton() && !rpa->gen.Beam1().IsMassive()
            && rpa->gen.Beam2().IsHadron()) {
        defpsmassless.push_back(rpa->gen.Beam1().Kfcode());
    }
    for (size_t i(0);i<allflavs.size();++i) {
      if (allflavs[i].IsDummy()) continue;
      size_t kf(allflavs[i].Kfcode());
      bool add(true);
      for (size_t j(0);j<defpsmassive.size();++j)
        if (kf==defpsmassive[j]) { add=false; break; }
      for (size_t j(0);j<defpsmassless.size();++j)
        if (kf==defpsmassless[j]) { add=false; break; }
      if (add)  defpsmassive.push_back(kf);
    }
  }
  else {
    for (size_t i(0);i<allflavs.size();++i) {
      if (allflavs[i].IsDummy()) continue;
      size_t kf(allflavs[i].Kfcode());
      bool add(true);
      for (size_t j(0);j<defpsmassive.size();++j)
        if (kf==defpsmassive[j]) { add=false; break; }
      for (size_t j(0);j<defpsmassless.size();++j)
        if (kf==defpsmassless[j]) { add=false; break; }
      if (add && allflavs[i].IsMassive())  defpsmassive.push_back(kf);
      if (add && !allflavs[i].IsMassive()) defpsmassless.push_back(kf);
    }
  }
  // then remove and add those specified manually
  for (size_t i(0);i<psmassive.size();++i) {
    defpsmassless.erase(std::remove(defpsmassless.begin(),defpsmassless.end(),
                                    psmassive[i]),defpsmassless.end());
    if (std::find(defpsmassive.begin(),defpsmassive.end(),psmassive[i])==
        defpsmassive.end()) defpsmassive.push_back(psmassive[i]);
  }
  for (size_t i(0);i<psmassless.size();++i) {
    defpsmassive.erase(std::remove(defpsmassive.begin(),defpsmassive.end(),
                                   psmassless[i]),defpsmassive.end());
    if (std::find(defpsmassless.begin(),defpsmassless.end(),psmassless[i])==
        defpsmassless.end()) defpsmassless.push_back(psmassless[i]);
  }
  // fill massive ones into m_psmass
  for (size_t i(0);i<defpsmassive.size();++i) {
    Flavour fl(defpsmassive[i],0);
    m_psmass.insert(fl);
    m_psmass.insert(fl.Bar());
    msg_Tracking()<<METHOD<<"(): "<<m_name<<": Using massive PS for "<<fl<<".\n";
  }
  Flavour_Vector mf;
  for (Flavour_Set::iterator fit(m_psmass.begin());fit!=m_psmass.end();++fit)
    if (fit->Mass(true)!=fit->Mass(false)) mf.push_back(*fit);
  msg_Info()<<"Massive PS flavours for "<<m_name<<": "
                    <<mf<<std::endl;
}

void ME_Generator_Base::RegisterDipoleParameters()
{
  // most (but not all) are used by both COMIX and AMEGIC
  // some are also used by PHASIC (e.g. KP_Terms)

  Scoped_Settings s{ Settings::GetMainSettings()["DIPOLES"] };
  s["AMIN"].SetDefault(Max(rpa->gen.Accu(), 1.0e-8));
  const auto& amax = s["ALPHA"].SetDefault(1.0).Get<double>();
  s["ALPHA_FF"].SetDefault(amax);
  s["ALPHA_FI"].SetDefault(amax);
  s["ALPHA_IF"].SetDefault(amax);
  s["ALPHA_II"].SetDefault(amax);
  s["NF_GSPLIT"].SetDefault(Flavour(kf_jet).Size() / 2);
  s["KT2MAX"].SetDefault(sqr(rpa->gen.Ecms()));
  s["COLLINEAR_VFF_SPLITTINGS"].SetDefault(1);
  s["V_SUBTRACTION_MODE"].SetDefault(1);  // 0: scalar, 1: fermionic
  s["PFF_IS_SPLIT_SCHEME"].SetDefault(1);
  s["PFF_FS_SPLIT_SCHEME"].SetDefault(0);
  s["PFF_IS_RECOIL_SCHEME"].SetDefault(0);
  s["PFF_FS_RECOIL_SCHEME"].SetDefault(0);
  s["IS_CLUSTER_TO_LEPTONS"].SetDefault(0);
  s["LIST"].SetDefault(0);
  s.DeclareVectorSettingsWithEmptyDefault({ "BORN_FLAVOUR_RESTRICTIONS" });
  s["ONSHELL_SUBTRACTION_WINDOW"].SetDefault(5.0);
}

void ME_Generator_Base::RegisterNLOParameters()
{
  SetParameter("NLO_SMEAR_THRESHOLD", 0.0);
  SetParameter("NLO_SMEAR_POWER", 0.5);
}

template <typename T>
void ME_Generator_Base::SetParameter(const std::string& param,
                                     const T& def)
{
  Scoped_Settings s{ Settings::GetMainSettings()[param] };
  s.SetDefault(def);
  rpa->gen.SetVariable(param, ToString(s.Get<T>()));
}

namespace PHASIC {

  class ShiftMasses_Energy: public Function_Base {
  private:
    std::vector<double>::size_type m_nentries;
    std::vector<double> m_m2, m_p2;
  public:
    ShiftMasses_Energy(Mass_Selector *const ms,
		    Cluster_Amplitude *const ampl,int mode)
    {
      const auto nin = ampl->NIn();
      auto offset = 0;
      if (mode < 0) {
        m_nentries = nin;
      } else {
        offset = nin;
        m_nentries = ampl->Legs().size() - nin;
      }
      m_p2.reserve(m_nentries);
      m_m2.reserve(m_nentries);
      const auto end = offset + m_nentries;
      for (int i {offset}; i < end; ++i) {
        m_p2.push_back(ampl->Leg(i)->Mom().PSpat2());
        m_m2.push_back(ms->Mass2(ampl->Leg(i)->Flav()));
      }
    }
    virtual double operator()(double x)
    {
      const auto x2=x*x;
      auto E=0.0;
      for (size_t i {0}; i < m_nentries; ++i)
	E+=sqrt(m_m2[i]+x2*m_p2[i]);
      return E;
    }
  };// end of class ShiftMasses_Energy

  class ShiftMasses_DIS: public Function_Base {
  private:
    std::vector<double>::size_type m_nentries;
    std::vector<double> m_m2, m_xy2, m_z;
    double m_pZplus, m_pZminus, m_pZreserve;
  public:
    ShiftMasses_DIS(Mass_Selector *const ms,
                         Cluster_Amplitude *const ampl)
    {
      const int nin = ampl->NIn();
      m_nentries = ampl->Legs().size()-nin;
      m_xy2.reserve(m_nentries);
      m_m2.reserve(m_nentries);
      m_z.reserve(m_nentries);
      m_pZplus = 0;
      m_pZminus = 0;
      for(int i {nin}; i<ampl->Legs().size(); i++) {
        if(ampl->Leg(i)->Flav().Strong()) {
          const Vec4D mom = ampl->Leg(i)->Mom();
          m_xy2.push_back(sqr(mom[1])+sqr(mom[2]));
          m_z.push_back(mom[3]);
          m_m2.push_back(ms->Mass2(ampl->Leg(i)->Flav()));
          if(mom[3] > 0) m_pZplus += abs(mom[3]);
          else           m_pZminus += abs(mom[3]);
        }
      }
      m_nentries = m_m2.size();
      m_pZreserve = std::min(m_pZplus,m_pZminus);
    }

    double scaledZ(double z, double xi) {
      /// scale the z-component such that the overall  p_z of the final state
      /// is conserved (given by DIS kinematics), and preserve the relative
      /// p_z ordering in both hemispheres
      double totZ = z > 0 ? m_pZplus : m_pZminus;
      double frac = z/totZ;
      return (1-xi)*frac*(totZ-m_pZreserve) + xi*z;
    }

    virtual double operator()(double xi)
    {
      const double xi2=xi*xi;
      double E=0.0;
      for (size_t i {0}; i < m_nentries; ++i) {
        E+=sqrt(xi2*m_xy2[i]+sqr(scaledZ(m_z[i],xi))+m_m2[i]);
      }
      return E;
    }
  };// end of class ShiftMasses_DIS
}// end of namespace PHASIC

int ME_Generator_Base::ShiftMasses(Cluster_Amplitude *const ampl)
{
  /// first figure out if ampl has a flavour we want to shift on-shell
  if (m_psmass.empty()) return 0;
  bool run=false;
  Vec4D cms;
  for (size_t i(0);i<ampl->Legs().size();++i) {
    if (i<ampl->NIn()) cms-=ampl->Leg(i)->Mom();
    if (m_psmass.find(ampl->Leg(i)->Flav())!=
	m_psmass.end()) run=true;
  }
  if (!run) return 1;
  /// if so treat DIS as special case
  if(ampl->NIn() <= 1 ||
     (!(ampl->Leg(0)->Flav().IsLepton() && ampl->Leg(1)->Flav().Strong()) &&
      !(ampl->Leg(0)->Flav().Strong() && ampl->Leg(1)->Flav().IsLepton()) ) ) {
    return ShiftMassesDefault(ampl, cms);
  }
  else {
    return ShiftMassesDIS(ampl, cms);
  }
}

int ME_Generator_Base::ShiftMassesDefault(Cluster_Amplitude *const ampl, Vec4D cms)
{
  DEBUG_FUNC(m_name);
  msg_Debugging()<<"Before shift: "<<*ampl<<"\n";
  Poincare boost(cms);
  boost.Boost(cms);
  for (size_t i(0);i<ampl->Legs().size();++i)
    ampl->Leg(i)->SetMom(boost*ampl->Leg(i)->Mom());
  boost.Invert();
  if (ampl->NIn()>1) {
    ShiftMasses_Energy etot(this,ampl,-1);
    double xi(etot.WDBSolve(cms[0],0.0,1.0));
    if (!IsEqual(etot(xi),cms[0],rpa->gen.Accu())) {
      if (m_massmode==0) xi=etot.WDBSolve(cms[0],1.0,2.0);
      if (!IsEqual(etot(xi),cms[0],rpa->gen.Accu())) return -1;
    }
    for (size_t i(0);i<ampl->NIn();++i) {
      Vec4D p(xi*ampl->Leg(i)->Mom());
      p[0]=-sqrt(Mass2(ampl->Leg(i)->Flav())+p.PSpat2());
      ampl->Leg(i)->SetMom(boost*p);
    }
  }
  ShiftMasses_Energy etot(this,ampl,1);
  double xi(etot.WDBSolve(cms[0],0.0,1.0));
  if (!IsEqual(etot(xi),cms[0],rpa->gen.Accu())) {
    if (m_massmode==0) xi=etot.WDBSolve(cms[0],1.0,2.0);
    if (!IsEqual(etot(xi),cms[0],rpa->gen.Accu())) return -1;
  }
  for (size_t i(ampl->NIn());i<ampl->Legs().size();++i) {
    Vec4D p(xi*ampl->Leg(i)->Mom());
    p[0]=sqrt(Mass2(ampl->Leg(i)->Flav())+p.PSpat2());
    ampl->Leg(i)->SetMom(boost*p);
  }
  msg_Debugging()<<"After shift: "<<*ampl<<"\n";
  return 1;
}


/// convenience function
Vec4D MomSum(Cluster_Amplitude *const ampl) {
  Vec4D ret(0,0,0,0);
  for(auto l: ampl->Legs()) ret += l->Mom();
  return ret;
}

int ME_Generator_Base::ShiftMassesDIS(Cluster_Amplitude *const ampl, Vec4D cms) {
  DEBUG_FUNC(m_name);
  msg_Debugging()<<"Before shift: "<<*ampl<<"\n";

  /// currently assume leg 0 is the electron/QCD siglet
  /// is this ever wrong?
  const Vec4D pLepIn = ampl->Leg(0)->Mom();

  BreitBoost breit(ampl);
  breit.Apply(ampl);
  msg_Debugging()<<"In Breit frame: "<<*ampl<<"\n";
  msg_Debugging()<<"Momentum conservation: "<<MomSum(ampl)<<".\n";
  breit.Invert();

  double Ein = 0;
  for(size_t i=0; i<ampl->NIn(); i++) {
    Vec4D p = ampl->Leg(i)->Mom();
    if(ampl->Leg(i)->Flav().Strong()) {
      p[0]=-sqrt(Mass2(ampl->Leg(i)->Flav())+p.PSpat2());
      Ein = -p[0];
    }
    ampl->Leg(i)->SetMom(p);
  }

  ShiftMasses_DIS etot(this,ampl);
  // need at least the energy to produce all masses
  // while preserving pZ
  double EoutMin = 0;
  for(size_t i=ampl->NIn(); i<ampl->Legs().size(); i++) {
    if(ampl->Leg(i)->Flav().Strong()) {
      EoutMin += sqrt(Mass2(ampl->Leg(i)->Flav()) +
                      sqr(etot.scaledZ(ampl->Leg(i)->Mom()[3],0)));
    }
  }
  if(Ein < EoutMin) {
    msg_Debugging()<<"Not enough energy, Ein = "<<Ein
                   <<" vs "<<EoutMin<<".\n";
    return -1;
  }
  double xi(etot.WDBSolve(Ein,0.0,1.0));
  for (size_t i(ampl->NIn());i<ampl->Legs().size();++i) {
    Vec4D p = ampl->Leg(i)->Mom();
    if(ampl->Leg(i)->Flav().Strong()) {
      p[1] *= xi;
      p[2] *= xi;
      p[3] = etot.scaledZ(p[3],xi);
      p[0] = sqrt(Mass2(ampl->Leg(i)->Flav())+p.PSpat2());
    }
    ampl->Leg(i)->SetMom(p);
  }
  msg_Debugging()<<"After shift (xi = "<<xi<<") in Breit frame: "<<*ampl<<"\n";
  msg_Debugging()<<"Momentum conservation: "<<MomSum(ampl)<<".\n";
  msg_Debugging()<<"DIS variables: Q2 = "<<breit.Q2()<<" vs "<<BreitBoost(ampl).Q2()
                 <<" and x = "<<breit.x()<<" vs "<<BreitBoost(ampl).x()<<".\n";
  breit.Apply(ampl);
  msg_Debugging()<<"After boost back: "<<*ampl<<"\n";
  msg_Debugging()<<"Momentum conservation: "<<MomSum(ampl)<<".\n";
  msg_Debugging()<<"DIS variables: Q2 = "<<breit.Q2()<<" vs "<<BreitBoost(ampl).Q2()
                 <<" and x = "<<breit.x()<<" vs "<<BreitBoost(ampl).x()<<".\n";
  /// shifted in momentum will be a bit of the z axis, so construct a new one that is
  /// aligned and push recoil to hadronic final state
  Vec4D pStrongOut(0,0,0,0);
  for (size_t i(ampl->NIn());i<ampl->Legs().size();++i) {
    if(ampl->Leg(i)->Flav().Strong()) pStrongOut += ampl->Leg(i)->Mom();
  }
  Vec4D pStrongIn;
  for(size_t i=0; i<ampl->NIn(); i++) {
    if(ampl->Leg(i)->Flav().Strong()) pStrongIn = -ampl->Leg(i)->Mom();
  }
  /// hadronic final state should have finite mass, otherwise we had a
  /// we should be in a 2->2 case with massless partons, and never arrive here
  if(!IsZero(pStrongOut.Abs2())) {
    /// might assume positive z for hadronic IS?
    const double m2 = pStrongIn.Abs2();
    const double inOutProd = pStrongIn*pStrongOut;
    const double eDiff = pStrongOut[0]-pStrongIn[0];
    const double zDiff = pStrongOut[3]-pStrongIn[3];
    double newPz = (inOutProd-m2)/(eDiff+zDiff);
    newPz *= zDiff - eDiff * sqrt(1+m2*(sqr(zDiff)-sqr(eDiff))/sqr(inOutProd-m2));
    newPz /= eDiff - zDiff;
    const Vec4D newInMom(sqrt(sqr(newPz)+m2),0,0,newPz);

    Poincare oldHCM(pStrongOut);
    Poincare newHCM(pStrongOut-pStrongIn+newInMom);

    for(size_t i=0; i<ampl->Legs().size(); i++) {
      if(i < ampl->NIn()) {
        if(ampl->Leg(i)->Flav().Strong()) ampl->Leg(i)->SetMom(-newInMom);
        else                              ampl->Leg(i)->SetMom(pLepIn);
      }
      else if(ampl->Leg(i)->Flav().Strong()) {
        Vec4D p = ampl->Leg(i)->Mom();
        oldHCM.Boost(p);
        newHCM.BoostBack(p);
        ampl->Leg(i)->SetMom(p);
      }
    }
  }
  msg_Debugging()<<"After full shift: "<<*ampl<<"\n";
  msg_Debugging()<<"Momentum conservation: "<<MomSum(ampl)<<".\n";
  if(!IsZero(MomSum(ampl).PSpat(),1e-6)) {
    msg_Error()<<"Mass shift could not be completed. "
               <<"Momentum conseravtion fails by "<<MomSum(ampl)<<"\n";
    return -1;
  }
  msg_Debugging()<<"Finished DIS mass shift with Q2 = "
                 <<breit.Q2()<<" vs "<<BreitBoost(ampl).Q2()<<" and x = "
                 <<breit.x()<<" vs "<<BreitBoost(ampl).x()<<".\n";
  return 1;
}


double ME_Generator_Base::Mass(const ATOOLS::Flavour &fl) const
{
  if (m_massmode==0) return fl.Mass();
  if (m_psmass.find(fl)!=m_psmass.end()) return fl.Mass(true);
  return fl.Mass();
}

void ME_Generator_Base::ShowSyntax(const int mode)
{
  if (!msg_LevelIsInfo() || mode==0) return;
  msg_Out()<<METHOD<<"(): {\n\n";
  ME_Generator_Getter::PrintGetterInfo(msg->Out(),15);
  msg_Out()<<"\n}"<<std::endl;
}

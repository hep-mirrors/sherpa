#include "PHASIC++/Process/ME_Generator_Base.H"

#include "PHASIC++/Process/ME_Generators.H"
#include "ATOOLS/Org/Run_Parameter.H"
#include "ATOOLS/Org/Scoped_Settings.H"
#include "ATOOLS/Math/Function_Base.H"
#include "ATOOLS/Math/Poincare.H"
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
  msg_Info()<<METHOD<<"(): Massive PS flavours for "<<m_name<<": "
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

}// end of namespace PHASIC

int ME_Generator_Base::ShiftMasses(Cluster_Amplitude *const ampl)
{
  if (m_psmass.empty()) return 0;
  bool run=false;
  Vec4D cms;
  for (size_t i(0);i<ampl->Legs().size();++i) {
    if (i<ampl->NIn()) cms-=ampl->Leg(i)->Mom();
    if (m_psmass.find(ampl->Leg(i)->Flav())!=
	m_psmass.end()) run=true;
  }
  if (!run) return 1;

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
  for (int i = 0; i < 2; i++) {
    if (p_remnant != NULL && !(p_remnant->GetRemnant(ampl->Leg(i)->Mom()[3] < 0.0 ? 0 : 1)
              ->TestExtract(ampl->Leg(i)->Flav().Bar(), -ampl->Leg(i)->Mom())))
      return -1;
  }
  msg_Debugging() << "After shift: " << *ampl << "\n";
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

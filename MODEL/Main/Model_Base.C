#include "MODEL/Main/Model_Base.H"

#define COMPILE__Getter_Function
#define OBJECT_TYPE MODEL::Model_Base
#define PARAMETER_TYPE MODEL::Model_Arguments
#include "ATOOLS/Org/Getter_Function.C"

#include "MODEL/Main/Single_Vertex.H"
#include "MODEL/Main/Running_AlphaS.H"
#include "MODEL/Main/Strong_Coupling.H"
#include "MODEL/Main/Running_Fermion_Mass.H"
#include "MODEL/Main/Running_AlphaQED.H"
#include "ATOOLS/Org/Run_Parameter.H"
#include "ATOOLS/Org/Message.H"
#include "ATOOLS/Org/Exception.H"
#include "ATOOLS/Org/Scoped_Settings.H"
#include "ATOOLS/Phys/KF_Table.H"

#include <algorithm>

using namespace MODEL;
using namespace ATOOLS;
using std::string;

namespace MODEL 
{
  Model_Base *s_model;
}

std::ostream &MODEL::operator<<(std::ostream &str,const ew_scheme::code &c)
{
  if      (c==ew_scheme::UserDefined) return str<<"UserDefined";
  else if (c==ew_scheme::alpha0)      return str<<"alpha0";
  else if (c==ew_scheme::alphamZ)     return str<<"alphamZ";
  else if (c==ew_scheme::Gmu)         return str<<"Gmu";
  else if (c==ew_scheme::alphamZsW)   return str<<"alphamZsW";
  else if (c==ew_scheme::alphamWsW)   return str<<"alphamWsW";
  else if (c==ew_scheme::GmumZsW)     return str<<"GmumZsW";
  else if (c==ew_scheme::GmumWsW)     return str<<"GmumWsW";
  else if (c==ew_scheme::FeynRules)   return str<<"FeynRules";
  return str<<"undefined";
}

std::istream &MODEL::operator>>(std::istream &str,ew_scheme::code &c)
{
  std::string tag;
  str>>tag;
  c=ew_scheme::Undefined;
  if      (tag=="UserDefined") c=ew_scheme::UserDefined;
  else if (tag=="0")           c=ew_scheme::UserDefined;
  else if (tag=="alpha0")      c=ew_scheme::alpha0;
  else if (tag=="1")           c=ew_scheme::alpha0;
  else if (tag=="alphamZ")     c=ew_scheme::alphamZ;
  else if (tag=="2")           c=ew_scheme::alphamZ;
  else if (tag=="Gmu")         c=ew_scheme::Gmu;
  else if (tag=="3")           c=ew_scheme::Gmu;
  else if (tag=="alphamZsW")   c=ew_scheme::alphamZsW;
  else if (tag=="4")           c=ew_scheme::alphamZsW;
  else if (tag=="alphamWsW")   c=ew_scheme::alphamWsW;
  else if (tag=="5")           c=ew_scheme::alphamWsW;
  else if (tag=="GmumZsW")     c=ew_scheme::GmumZsW;
  else if (tag=="6")           c=ew_scheme::GmumZsW;
  else if (tag=="GmumWsW")     c=ew_scheme::GmumWsW;
  else if (tag=="7")           c=ew_scheme::GmumWsW;
  else if (tag=="FeynRules")   c=ew_scheme::FeynRules;
  else if (tag=="10")          c=ew_scheme::FeynRules;
  else                         c=ew_scheme::Undefined;
  return str;
}

Model_Base::Model_Base(bool _elementary) :
  m_elementary(_elementary),
  m_hasnegativecouplingorders(false),
  p_numbers(NULL), p_constants(NULL), p_complexconstants(NULL),
  p_functions(NULL)
{
  p_numbers          = new MODEL::ScalarNumbersMap();
  p_constants        = new MODEL::ScalarConstantsMap();
  p_complexconstants = new MODEL::ComplexConstantsMap();
  p_functions        = new MODEL::ScalarFunctionsMap();
}

Model_Base::~Model_Base() 
{
  if (p_numbers!=NULL) delete p_numbers;
  if (p_functions!=NULL) {
    while (!p_functions->empty()) {
      delete p_functions->begin()->second;
      p_functions->erase(p_functions->begin());
    }
    delete p_functions;
  }
  if (p_constants!=NULL)         delete p_constants;
  if (p_complexconstants!=NULL)  delete p_complexconstants;
}

void Model_Base::RegisterDefaults() const
{
  Settings& s = Settings::GetMainSettings();
  s["1/ALPHAQED(MZ)"].SetDefault(128.802);
  s["1/ALPHAQED(MW)"].SetDefault(132.17);
  s["1/ALPHAQED(0)"].SetDefault(137.03599976);
  s["ALPHAS(MZ)"].SetDefault(0.118);
  s["ORDER_ALPHAS"].SetDefault(2);
  s["THRESHOLD_ALPHAS"].SetDefault(1);
  s["Q2_AS"].SetDefault(1.0);
  s["AS_FORM"].SetDefault("Smooth");
  s["CKM"]["Output"].SetDefault(false);
  s["JET_MASS_THRESHOLD"].SetDefault(10.0);
  s["YUKAWA_MASSES"].SetDefault("Running");
  s["WIDTH_SCHEME"].SetDefault("CMS");
  s["SIN2THETAW"].SetDefault(0.23155);
  s["VEV"].SetDefault(246.0);
  s["GF"].SetDefault(1.16639e-5);
  s["GMU_CMS_AQED_CONVENTION"].SetDefault(0);
  s["CKM"]["Order"].SetDefault(0);
  s["CKM"]["Cabibbo"].SetDefault(0.22537);
  s["CKM"]["A"].SetDefault(0.814);
  s["CKM"]["Eta"].SetDefault(0.353);
  s["CKM"]["Rho"].SetDefault(0.117);
  s["DECOMPOSE_4G_VERTEX"].SetDefault(1);
  s["EW_SCHEME"].SetDefault(ew_scheme::Gmu);
  const ew_scheme::code ewscheme( s["EW_SCHEME"].Get<ew_scheme::code>() );
  s["EW_REN_SCHEME"].SetDefault(ewscheme);
  if (ewscheme == ew_scheme::UserDefined)
    s["ALPHAQED_DEFAULT_SCALE"].SetDefault(sqr(Flavour(kf_Z).Mass()));
  else
    s["ALPHAQED_DEFAULT_SCALE"].SetDefault(0.0);
}

void Model_Base::RotateVertices()
{
  int nv=m_v.size(); 
  for (int i=0;i<nv;++i) {
    std::vector<size_t> id(m_v[i].id);
    if (m_v[i].dec>=0) {
      for (size_t k=0;k<m_v[i].in.size();++k) {
	Flavour fl(m_v[i].in[k]);
	if (m_maxlegs.find(fl)==m_maxlegs.end()) m_maxlegs[fl]=0;
	if (id.size()>m_maxlegs[fl]) m_maxlegs[fl]=id.size();
      }
    }
    if (m_v[i].dec&2) continue;
    for (size_t k=0;k<id.size()-1;++k) {
      for (int lid=id[id.size()-1], l=id.size()-1;l>=0;--l) id[l]=l?id[l-1]:lid;
      Single_Vertex v(m_v[i]);
      for (int j=0;j<v.in.size();++j) v.in[j]=m_v[i].in[v.id[j]=id[j]];
      if(find(m_v.begin(),m_v.end(),v)==m_v.end()) m_v.push_back(v);
    }
  }
}

void Model_Base::GetCouplings(Coupling_Map &cpls)
{
  DEBUG_FUNC(&cpls);
  for (ScalarFunctionsMap::const_iterator
	 cit(p_functions->begin());cit!=p_functions->end();++cit) {
    std::string tag(cit->second->Name());
    cpls.insert(std::make_pair(tag,new Coupling_Data
      (cit->second,tag,ScalarConstant(cit->first))));
    msg_Debugging()<<"  '"<<tag<<"' -> ("<<cpls.lower_bound(tag)->second<<")"
		   <<*cpls.lower_bound(tag)->second<<"\n";
  }
}

// To be called in ModelInit, default value will be set to aqed_def argument
void Model_Base::SetAlphaQED(const double& aqed_def)
{
  Settings& s = Settings::GetMainSettings();
  double alphaQED0=1./s["1/ALPHAQED(0)"].Get<double>();
  aqed=new Running_AlphaQED(alphaQED0);
  aqed->SetDefault(aqed_def);
  p_functions->insert(make_pair(std::string("alpha_QED"),aqed));
  p_constants->insert(make_pair(std::string("alpha_QED"),aqed_def));
}

// To be called in ModelInit, default will be set to AlphaQED at scale2
void Model_Base::SetAlphaQEDByScale(const double& scale2)
{
  Settings& s = Settings::GetMainSettings();
  double alphaQED0=1./s["1/ALPHAQED(0)"].Get<double>();
  aqed=new Running_AlphaQED(alphaQED0);
  aqed->SetDefault((*aqed)(scale2));
  p_functions->insert(make_pair(std::string("alpha_QED"),aqed));
  p_constants->insert(make_pair(std::string("alpha_QED"),aqed->Default()));
}

// To be called in ModelInit, default will be set to AlphaQED as input
void Model_Base::SetAlphaQEDByInput(const std::string& tag)
{
  Settings& s = Settings::GetMainSettings();
  double alphaQED0=1./s[tag].Get<double>();
  aqed=new Running_AlphaQED(alphaQED0);
  p_functions->insert(make_pair(std::string("alpha_QED"),aqed));
  p_constants->insert(make_pair(std::string("alpha_QED"),aqed->Default()));
}

// To be called in ModelInit, alphaS argument is alphaS input at MZ
void Model_Base::SetAlphaQCD(const PDF::ISR_Handler_Map& isr, const double& alphaS)
{
  Settings& s = Settings::GetMainSettings();
  int    order_alphaS   = s["ORDER_ALPHAS"].Get<int>();
  int    th_alphaS      = s["THRESHOLD_ALPHAS"].Get<int>();
  double MZ2            = sqr(Flavour(kf_Z).Mass());
  as = new Running_AlphaS(alphaS,MZ2,order_alphaS,th_alphaS,isr);
  p_constants->insert(make_pair(string("alpha_S"),alphaS));
  p_functions->insert(make_pair(string("alpha_S"),as));
  double Q2aS = s["Q2_AS"].Get<int>();
  string asf  = s["AS_FORM"].Get<string>();
  asform::code as_form(asform::smooth);
  if (asf==string("Constant"))    as_form = asform::constant;
  else if (asf==string("Frozen")) as_form = asform::frozen;
  else if (asf==string("Smooth")) as_form = asform::smooth;
  else if (asf==string("IR0"))    as_form = asform::IR0;
  else if (asf==string("GDH"))    as_form = asform::GDH_inspired;
  Strong_Coupling * strong_cpl(new Strong_Coupling(as,as_form,Q2aS));
  p_functions->insert(make_pair(string("strong_cpl"),strong_cpl));
  p_constants->insert(make_pair(string("strong_cpl"),alphaS));
}

// To be called in ModelInit 
void Model_Base::SetRunningFermionMasses()
{
  for (size_t i=0;i<17; ++i) {
    if (i==7) i=11;
    Flavour yfl((kf_code)i);
    if (yfl.Yuk()==0.0) continue;
    Running_Fermion_Mass *rfm(new Running_Fermion_Mass(yfl, yfl.Yuk(), as));
    p_functions->insert(make_pair("m"+yfl.IDName(),rfm));
    p_constants->insert(make_pair("m"+yfl.IDName(),yfl.Yuk()));
  }
}

void Model_Base::SetRunningBosonMasses()
{
  for (size_t i=23;i<26; ++i) {
    Flavour yfl((kf_code)i);
    if (yfl.Yuk()==0.0) continue;
    p_constants->insert(make_pair("m"+yfl.IDName(),yfl.Yuk()));
  }
}

void Model_Base::ReadExplicitCKM(CMatrix& CKM)
{
  Settings& s = Settings::GetMainSettings();
  for (const auto& key : s["CKM"]["Matrix_Elements"].GetKeys()) {
    const auto indizes = ToVector<int>(key, ',');
    if (indizes.size() != 2) {
      THROW(fatal_error, "Please give two comma-separated indizes for each"
                         + std::string(" CKM:Matrix_Elements definition."));
    }
    if (indizes[0] > CKM.Rank() || indizes[1] > CKM.Rank())
      THROW(fatal_error, "Trying to read in CKM element beyond range.");
    const auto values = s["CKM"]["Matrix_Elements"][key]
      .SetDefault({0.0, 0.0})
      .GetVector<double>();
    if (values.size() == 1) {
      CKM[indizes[0]][indizes[1]] = Complex(values[0], 0.0);
    } else if (values.size() == 2) {
      CKM[indizes[0]][indizes[1]] = Complex(values[0], values[1]);
    } else {
      THROW(fatal_error, "Please provide either one or two values for each."
                         + std::string(" CKM:Matrix_Elements definition."));
    }
  }
}

void Model_Base::OutputCKM()
{
  Settings& s = Settings::GetMainSettings();
  if (!s["CKM"]["Output"].Get<bool>())
    return;
  msg_Info()<<" CKM Matrix:\n";
  msg_Info()<<std::setw(8)<<"V_ij";
  size_t rank(ScalarConstant("CKM_DIMENSION"));
  size_t width(96/rank);
  // flavour output assumes quarks are all consecutive from 1...n
  for (size_t i(0);i<rank;++i)
    msg_Info()<<std::setw(width)<<Flavour(2*i+1);
  msg_Info()<<std::endl;
  for (size_t i(0);i<rank;++i) {
    msg_Info()<<std::setw(8)<<Flavour(2*i+2);
    for (size_t j(0);j<rank;++j) {
      msg_Info()<<std::setw(width)
                <<ComplexConstant("CKM_"+ToString(i)+"_"+ToString(j));
    }
    msg_Info()<<std::endl;
  }
  msg_Info()<<std::endl;
}

void Model_Base::ShowSyntax(const size_t i)
{
  if (!msg_LevelIsInfo() || i==0) return;
  msg_Out()<<METHOD<<"(): {\n\n"
	   <<"   // available model implementations (specified by MODEL: <value>)\n\n";
  Model_Getter_Function::PrintGetterInfo(msg->Out(),25);
  msg_Out()<<"\n}"<<std::endl;
}

void Model_Base::ReadParticleData()
{
  auto pdata = Settings::GetMainSettings()["PARTICLE_DATA"];
  for (const auto& ptclname : pdata.GetKeys()) {
    kf_code kf = ToType<kf_code>(ptclname);
    const auto it = s_kftable.find(kf);
    if (it != s_kftable.end()) {
      for (const auto& propertyname : pdata[ptclname].GetKeys()) {
        auto s = pdata[ptclname][propertyname];
        if (propertyname == "Mass") {
          const auto mass = s.SetDefault(it->second->m_mass).Get<double>();
          it->second->m_mass = mass;
          it->second->m_hmass = mass;
        } else if (propertyname == "Width") {
          it->second->m_width = s.SetDefault(it->second->m_width).Get<double>();
        } else if (propertyname == "Active") {
          it->second->m_on = s.SetDefault(it->second->m_on).Get<bool>();
        } else if (propertyname == "Stable") {
          it->second->m_stable = s.SetDefault(it->second->m_stable).Get<int>();
        } else if (propertyname == "Massive") {
          it->second->m_massive = s.SetDefault(it->second->m_massive).Get<bool>();
        } else if (propertyname == "IntCharge") {
          it->second->m_icharge = s.SetDefault(it->second->m_icharge).Get<int>();
        } else if (propertyname == "StrongCharge") {
          it->second->m_strong = s.SetDefault(it->second->m_strong).Get<int>();
        } else if (propertyname == "Yukawa") {
          it->second->m_yuk = s.SetDefault(it->second->m_yuk).Get<double>();
        } else if (propertyname == "Priority") {
          it->second->m_priority = s.SetDefault(it->second->m_priority).Get<int>();
        } else {
          THROW(fatal_error,
                "Unknown PARTICLE_DATA property '" + propertyname + "'");
        }
      }
    }
  }
}

void Model_Base::AddStandardContainers()
{
  if (s_kftable.find(kf_jet) != s_kftable.end()) {
    return;
  }
  // kf,mass,width,icharge,strong,spin,majo,on,stable,massive,
  //   idname,antiidname,texname,antitexname,dummy,group
  s_kftable[kf_jet] = new
    Particle_Info(kf_jet,0.,0.,0.,0,1, 2,1,1,1,0,"j","j","j","j",1,1);
  s_kftable[kf_ewjet] = new
    Particle_Info(kf_ewjet,0.,0.,0.,0,1, 2,1,1,1,0,"ewj","ewj","ewj","ewj",1,1);
  s_kftable[kf_quark] = new
    Particle_Info(kf_quark,0.,0.,0.,0,1,1,0,1,1,0,"Q","Q","Q","Q",1,1);
  s_kftable[kf_lepton] = new
    Particle_Info(kf_lepton,0.,0.,0.,-3,0,1,0,1,1,0,"l","l","\\ell","\\ell",1,1);
  s_kftable[kf_neutrino] = new
    Particle_Info(kf_neutrino,0.,0.,0.,0,0,1,0,1,1,0,"v","v","\\nu","\\nu",1,1);
  s_kftable[kf_fermion] = new
    Particle_Info(kf_fermion,0.,0.,0.,0,0,1,0,1,1,0,"f","f","f","f",1,1);
  s_kftable[kf_lepton]->m_priority=2;
  s_kftable[kf_neutrino]->m_priority=1;
  s_kftable[kf_jet]->Clear();
  s_kftable[kf_ewjet]->Clear();
  s_kftable[kf_quark]->Clear();
  s_kftable[kf_lepton]->Clear();
  s_kftable[kf_neutrino]->Clear();
  s_kftable[kf_fermion]->Clear();
  Settings& s = Settings::GetMainSettings();
  const double jet_mass_threshold{ s["JET_MASS_THRESHOLD"].Get<double>() };
  for (int i=1;i<7;i++) {
    Flavour addit((kf_code)i);
    if ((addit.Mass()==0.0 || !addit.IsMassive()) && addit.IsOn()) {
      if (addit.Mass(true)<=jet_mass_threshold) {
        s_kftable[kf_jet]->Add(addit);
        s_kftable[kf_jet]->Add(addit.Bar());
        s_kftable[kf_ewjet]->Add(addit);
        s_kftable[kf_ewjet]->Add(addit.Bar());
        s_kftable[kf_quark]->Add(addit);
        s_kftable[kf_quark]->Add(addit.Bar());
        s_kftable[kf_fermion]->Add(addit);
        s_kftable[kf_fermion]->Add(addit.Bar());
      }
      else {
        msg_Info()<<"Ignoring "<<addit<<" due to JET_MASS_THRESHOLD.\n";
      }
    }
  }
  s_kftable[kf_jet]->Add(Flavour(kf_gluon));
  s_kftable[kf_jet]->SetResummed();
  s_kftable[kf_ewjet]->Add(Flavour(kf_gluon));
  s_kftable[kf_ewjet]->Add(Flavour(kf_photon));
  s_kftable[kf_ewjet]->SetResummed();
  for (int i=11;i<17;i+=2) {
    Flavour addit((kf_code)i);
    if ((addit.Mass()==0.0 || !addit.IsMassive()) && addit.IsOn()) {
      s_kftable[kf_lepton]->Add(addit);
      s_kftable[kf_lepton]->Add(addit.Bar());
      s_kftable[kf_fermion]->Add(addit);
      s_kftable[kf_fermion]->Add(addit.Bar());
      if (s_kftable[i]->m_priority)
	msg_Error()<<METHOD<<"(): Changing "<<addit<<" sort priority: "
		   <<s_kftable[i]->m_priority<<" -> "
		   <<s_kftable[kf_lepton]->m_priority<<std::endl;
      s_kftable[i]->m_priority=s_kftable[kf_lepton]->m_priority;
    }
  }
  for (int i=12;i<17;i+=2) {
    Flavour addit((kf_code)i);
    if ((addit.Mass()==0.0) && addit.IsOn()) {
      s_kftable[kf_neutrino]->Add(addit);
      s_kftable[kf_neutrino]->Add(addit.Bar());
      s_kftable[kf_fermion]->Add(addit);
      s_kftable[kf_fermion]->Add(addit.Bar());
      if (s_kftable[i]->m_priority)
	msg_Error()<<METHOD<<"(): Changing "<<addit<<" sort priority: "
		   <<s_kftable[i]->m_priority<<" -> "
		   <<s_kftable[kf_neutrino]->m_priority<<std::endl;
      s_kftable[i]->m_priority=s_kftable[kf_neutrino]->m_priority;
    }
  }
}

void Model_Base::CustomContainerInit()
{
  DEBUG_FUNC("");
  auto s = Settings::GetMainSettings()["PARTICLE_CONTAINERS"];
  for (const auto& containeridstring : s.GetKeys()) {
    auto cs = s[containeridstring];
    const auto containerid = ToType<long int>(containeridstring);
    if (s_kftable.find(containerid) != s_kftable.end()) {
      THROW(critical_error,
            "Particle ID " + containeridstring + " already exists.");
    }
    auto name = cs["Name"].SetDefault(containeridstring).Get<std::string>();
    auto barname = name;
    if (name.find('|')!=std::string::npos) {
      name=name.substr(0,name.find('|'));
      barname=barname.substr(barname.find('|')+1);
    }
    const auto majorana = cs["Majorana"].SetDefault(false).Get<bool>();
    if (barname == name && !majorana)
      barname+="b";
    s_kftable[containerid]
      = new Particle_Info(containerid,
                          cs["Mass"].SetDefault(0.0).Get<double>(),
                          cs["Radius"].SetDefault(0.0).Get<double>(),
                          cs["Width"].SetDefault(0.0).Get<double>(),
                          cs["ICharge"].SetDefault(0).Get<int>(),
                          cs["Strong"].SetDefault(0).Get<int>(),
                          cs["Spin"].SetDefault(0).Get<int>(),
                          majorana,
                          1, 1, 0, name, barname, name, barname);
    s_kftable[containerid]->m_priority
      = cs["Priority"].SetDefault(0).Get<int>();
    s_kftable[containerid]->Clear();
    const auto flavs
      = cs["Flavs"].SetSynonyms({"Flavours","Flavors"}).SetDefault<long int>({}).GetVector<long int>();
    for (const auto flav : flavs) {
      s_kftable[containerid]->Add(Flavour((kf_code)std::abs(flav), flav < 0));
      if (s_kftable[std::abs(flav)]->m_priority) {
        msg_Error()<<METHOD<<"(): Changing "<<Flavour(flav)<<" sort priority: "
                   <<s_kftable[std::abs(flav)]->m_priority<<" -> "
                   <<s_kftable[containerid]->m_priority<<std::endl;
        s_kftable[std::abs(flav)]->m_priority
          = s_kftable[containerid]->m_priority;
      }
    }
    s_kftable[containerid]->SetIsGroup(true);
  }
}

void Model_Base::InitializeInteractionModel()
{
  InitVertices();
  for (std::vector<Single_Vertex>::iterator
	 vit(m_v.begin());vit!=m_v.end();) {
    for (size_t i(0);i<vit->cpl.size();)
      if (vit->cpl[i].Value().real()==0.0 &&
	  vit->cpl[i].Value().imag()==0.0) {
	vit->cpl.erase(vit->cpl.begin()+i);
	vit->Color.erase(vit->Color.begin()+i);
	vit->Lorentz.erase(vit->Lorentz.begin()+i);
      }
      else { ++i; }
    if (vit->cpl.empty()) vit=m_v.erase(vit);
    else ++vit;
  }
  CheckForNegativeCouplingOrders();
  m_ov=m_v;
  RotateVertices();
  InitMEInfo();
}

int Model_Base::ScalarNumber(const std::string _name) const {
  if (p_numbers->count(_name)>0) return (*p_numbers)[_name];
  THROW(fatal_error, "Key "+_name+" not found");
}


double Model_Base::ScalarConstant(const std::string _name) const {
  if (p_constants->count(_name)>0) return (*p_constants)[_name];
  THROW(fatal_error, "Key "+_name+" not found");
}


Complex Model_Base::ComplexConstant(const std::string _name) const {
  if (p_complexconstants->count(_name)>0) return (*p_complexconstants)[_name];
  THROW(fatal_error, "Key "+_name+" not found");
}


Function_Base * Model_Base::GetScalarFunction(const std::string _name) {
  if (p_functions->count(_name)>0) return (*p_functions)[_name];
  THROW(fatal_error, "Key "+_name+" not found");
}


double Model_Base::ScalarFunction(const std::string _name,double _t) {
  if (p_functions->count(_name)>0) return (*(*p_functions)[_name])(_t);
  THROW(fatal_error, "Key "+_name+" not found");
}


bool Model_Base::CheckFlavours(int nin, int nout, Flavour* flavs)
{
  return true;
}

void Model_Base::InitMEInfo()
{
  msg_Debugging()<<METHOD<<"(): {\n";
  m_fls.clear();
  std::set<Flavour> fls;
  msg_Debugging()<<"\n  add vertices\n\n";
  std::vector<Single_Vertex> &all(m_v);
  for (size_t i=0;i<all.size();++i) {
    m_vmap.insert(VMap_Key(all[i].PID(),&all[i]));
    m_vtable[all[i].in[0].Bar()].push_back(&all[i]);
    for (int j(0);j<all[i].NLegs();++j) fls.insert(all[i].in[j]);
    if (msg_LevelIsDebugging()) {
      msg_Debugging()
        <<"  "<<all[i].PID()<<" ["<<all[i].id[0];
      for (size_t j(1);j<all[i].id.size();++j) msg_Out()<<","<<all[i].id[j];
      msg_Out()<<"] "<<all[i].order<<" "<<(all[i].dec>0?'{':(all[i].dec<0?'(':'['))
               <<all[i].Lorentz.front()<<","<<all[i].Color[0].PID();
      for (size_t j(1);j<all[i].Lorentz.size();++j)
        msg_Out()<<"|"<<all[i].Lorentz[j]<<","<<all[i].Color[j].PID();
      msg_Out()<<(all[i].dec>0?'}':(all[i].dec<0?')':']'));
      for (size_t l(0);l<all[i].cpl.size();++l)
          msg_Out()<<", C"<<l<<" = "<<all[i].Coupling(l);
      msg_Out()<<"\n";
    }
  }
  msg_Debugging()<<"\n  add particles\n\n";
  for (std::set<Flavour>::const_iterator 
	 fit(fls.begin());fit!=fls.end();++fit) {
      m_fls.push_back(*fit);
      msg_Debugging()<<"  "<<*fit<<"\n";
  }
  msg_Debugging()<<"\n}\n";
}

int Model_Base::MaxNumber() const
{
  return m_v.size();
}

const std::vector<Single_Vertex> &Model_Base::Vertices() const
{
  return m_v;
}

const std::vector<Single_Vertex> &Model_Base::OriginalVertices() const
{
  return m_ov;
}

size_t Model_Base::IndexOfOrderKey(const std::string& key) const
{
  if (key == "QCD")
    return 0;
  if (key == "EW")
    return 1;
  THROW(fatal_error, "Unknown Orders key '" + key + "'.");
}

void Model_Base::CheckForNegativeCouplingOrders()
{
  for (std::vector<MODEL::Single_Vertex>::const_iterator vit(m_v.begin());
       vit!=m_v.end();++vit) {
    for (std::vector<int>::const_iterator oit(vit->order.begin());
         oit!=vit->order.end();++oit) {
      if (*oit<0) {
        m_hasnegativecouplingorders=true;
        return;
      }
    }
  }
}

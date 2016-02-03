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
#include "ATOOLS/Org/Data_Reader.H"
#include "ATOOLS/Org/Run_Parameter.H"
#include "ATOOLS/Org/Message.H"
#include "ATOOLS/Org/Exception.H"

#include <algorithm>

using namespace MODEL;
using namespace ATOOLS;
using std::string;

namespace MODEL 
{
  Model_Base *s_model;
}

Model_Base::Model_Base(std::string _dir,std::string _file,bool _elementary) :
  m_dir(_dir), m_file(_file), m_elementary(_elementary), 
  p_dataread(NULL), p_numbers(NULL), p_constants(NULL), p_complexconstants(NULL), 
  p_functions(NULL)
{
  p_dataread = new Data_Reader(" ",";","!","=");
  p_dataread->AddComment("#");
  p_dataread->AddWordSeparator("\t");
  p_dataread->SetInputPath(m_dir);
  p_dataread->SetInputFile(m_file);

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
  if (p_dataread!=NULL)          delete p_dataread;
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
void Model_Base::SetAlphaQED(const double& aqed_def){
  double alphaQED0=1./p_dataread->GetValue<double>("1/ALPHAQED(0)",137.03599976);
  aqed=new Running_AlphaQED(alphaQED0);
  aqed->SetDefault(aqed_def);
  p_functions->insert(make_pair(std::string("alpha_QED"),aqed));
  p_constants->insert(make_pair(std::string("alpha_QED"),aqed_def));
}

// To be called in ModelInit, default will be set to AlphaQED at scale2
void Model_Base::SetAlphaQEDByScale(const double& scale2){
  double alphaQED0=1./p_dataread->GetValue<double>("1/ALPHAQED(0)",137.03599976);;
  aqed=new Running_AlphaQED(alphaQED0);
  aqed->SetDefault((*aqed)(scale2));
  p_functions->insert(make_pair(std::string("alpha_QED"),aqed));
  p_constants->insert(make_pair(std::string("alpha_QED"),aqed->Default()));
}

// To be called in ModelInit, alphaS argument is alphaS input at MZ
void Model_Base::SetAlphaQCD(const PDF::ISR_Handler_Map& isr, const double& alphaS)
{
  int    order_alphaS	= p_dataread->GetValue<int>("ORDER_ALPHAS",1);
  int    th_alphaS	= p_dataread->GetValue<int>("THRESHOLD_ALPHAS",1);
  double MZ2            = sqr(Flavour(kf_Z).Mass());
  as = new Running_AlphaS(alphaS,MZ2,order_alphaS,th_alphaS,isr);
  p_constants->insert(make_pair(string("alpha_S"),alphaS));
  p_functions->insert(make_pair(string("alpha_S"),as));
  double Q2aS      = p_dataread->GetValue<double>("Q2_AS",1.);
  string asf  = p_dataread->GetValue<string>("As_Form",string("smooth"));
  asform::code as_form(asform::smooth);
  if (asf==string("constant"))    as_form = asform::constant;
  else if (asf==string("frozen")) as_form = asform::frozen;
  else if (asf==string("smooth")) as_form = asform::smooth;
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

void Model_Base::ReadExplicitCKM(CMatrix& CKM)
{
  Data_Reader dr(" ",";","!","=");
  dr.AddComment("#");
  dr.AddWordSeparator("\t");
  dr.AddIgnore("[");
  dr.AddIgnore("]");
  dr.AddIgnore("(");
  dr.AddIgnore(")");
  dr.AddIgnore(",");
  dr.SetAddCommandLine(true);
  dr.SetInputPath(m_dir);
  dr.SetInputFile(m_file);
  std::map<size_t,std::map<size_t,Complex> > ckms;
  std::vector<std::vector<double> > helpdvv;
  if (dr.MatrixFromFile(helpdvv,"CKM_ELEMENT")) {
    for (size_t i(0);i<helpdvv.size();++i) {
      if (helpdvv[i][0]>CKM.Rank() || helpdvv[i][1]>CKM.Rank())
        THROW(fatal_error,"Trying to read in CKM element beyond range.");
      if (helpdvv[i].size()==3)
        CKM[size_t(helpdvv[i][0])][size_t(helpdvv[i][1])]
          =Complex(helpdvv[i][2],0.);
      else if (helpdvv[i].size()==4) {
        CKM[size_t(helpdvv[i][0])][size_t(helpdvv[i][1])]
          =Complex(helpdvv[i][2],helpdvv[i][3]);
      }
      else
        THROW(fatal_error,"Unrecognised CKM_ELEMENT input.");
    }
  }
}

void Model_Base::OutputCKM()
{
  if (!p_dataread->GetValue<int>("CKM_OUTPUT",0)) return;
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
	   <<"   // available model implementations (specified by MODEL=<value>)\n\n";
  Model_Getter_Function::PrintGetterInfo(msg->Out(),25);
  msg_Out()<<"\n}"<<std::endl;
}

void Model_Base::ReadParticleData() {
  std::map<int,double> cdm, cdw, cdy;
  std::map<int,int> cia, cis, cim, cic, csc, cip;
  Data_Reader dr(" ",";","!","=");
  dr.AddComment("#");
  dr.AddWordSeparator("\t");
  dr.AddIgnore("[");
  dr.AddIgnore("]");
  dr.SetAddCommandLine(true);
  dr.SetInputPath(m_dir);
  dr.SetInputFile(m_file);
  std::vector<std::vector<double> > helpdvv;
  if (dr.MatrixFromFile(helpdvv,"MASS"))
    for (size_t i(0);i<helpdvv.size();++i)
      if (helpdvv[i].size()==2) cdm[int(helpdvv[i][0])]=helpdvv[i][1];
  if (dr.MatrixFromFile(helpdvv,"WIDTH")) 
    for (size_t i(0);i<helpdvv.size();++i) 
      if (helpdvv[i].size()==2) cdw[int(helpdvv[i][0])]=helpdvv[i][1];
  if (dr.MatrixFromFile(helpdvv,"ACTIVE"))
    for (size_t i(0);i<helpdvv.size();++i)
      if (helpdvv[i].size()==2) cia[int(helpdvv[i][0])]=int(helpdvv[i][1]);
  if (dr.MatrixFromFile(helpdvv,"STABLE"))
    for (size_t i(0);i<helpdvv.size();++i)
      if (helpdvv[i].size()==2) cis[int(helpdvv[i][0])]=int(helpdvv[i][1]);
  if (dr.MatrixFromFile(helpdvv,"MASSIVE"))
    for (size_t i(0);i<helpdvv.size();++i)
      if (helpdvv[i].size()==2) cim[int(helpdvv[i][0])]=int(helpdvv[i][1]);
  if (dr.MatrixFromFile(helpdvv,"INTCHARGE"))
    for (size_t i(0);i<helpdvv.size();++i)
      if (helpdvv[i].size()==2) cic[int(helpdvv[i][0])]=int(helpdvv[i][1]);
  if (dr.MatrixFromFile(helpdvv,"STRONGCHARGE"))
    for (size_t i(0);i<helpdvv.size();++i)
      if (helpdvv[i].size()==2) csc[int(helpdvv[i][0])]=int(helpdvv[i][1]);
  if (dr.MatrixFromFile(helpdvv,"YUKAWA"))
    for (size_t i(0);i<helpdvv.size();++i)
      if (helpdvv[i].size()==2) cdy[int(helpdvv[i][0])]=helpdvv[i][1];
  if (dr.MatrixFromFile(helpdvv,"PRIORITY"))
    for (size_t i(0);i<helpdvv.size();++i)
      if (helpdvv[i].size()==2) cip[int(helpdvv[i][0])]=int(helpdvv[i][1]);

  //set masses
  std::map<int,double>::const_iterator dit=cdm.begin();
  for (;dit!=cdm.end();dit++) {
    if (s_kftable.find(dit->first)!=s_kftable.end()) {
      s_kftable[dit->first]->m_mass = dit->second;
      msg_Tracking()<<" set mass of "<<Flavour(dit->first)<<" to "<<dit->second<<" GeV"<<std::endl; 
    }
  }
  //set widths
  dit=cdw.begin();
  for (;dit!=cdw.end();dit++) {
    if (s_kftable.find(dit->first)!=s_kftable.end()) {
      s_kftable[dit->first]->m_width = dit->second;
      msg_Tracking()<<" set width of "<<Flavour(dit->first)<<" to "<<dit->second<<" GeV"<<std::endl; 
    }
  }
  //set (in)active
  std::map<int,int>::const_iterator iit=cia.begin();
  for (;iit!=cia.end();iit++) {
    if (s_kftable.find(iit->first)!=s_kftable.end()) {
      s_kftable[iit->first]->m_on = iit->second;
      if (iit->second==0) {
	msg_Tracking()<<" set flavour "<<Flavour(iit->first)<<" inactive "<<std::endl; 
      }
      else {
	msg_Tracking()<<" set flavour "<<Flavour(iit->first)<<" active "<<std::endl; 
      }
    }
  }
  //set (un)stable
  iit=cis.begin();
  for (;iit!=cis.end();iit++) {
    if (s_kftable.find(iit->first)!=s_kftable.end()) {
      s_kftable[iit->first]->m_stable = iit->second;
      if (iit->second==0) {
	msg_Tracking()<<" set flavour "<<Flavour(iit->first)<<" unstable "<<std::endl; 
      }
      else {
	msg_Tracking()<<" set flavour "<<Flavour(iit->first)<<" stable "<<std::endl; 
      }
    }
  }
  //set massive/massless
  iit=cim.begin();
  for (;iit!=cim.end();iit++) {
    if (s_kftable.find(iit->first)!=s_kftable.end()) {
      s_kftable[iit->first]->m_massive = iit->second;
      if (iit->second==0) {
	msg_Tracking()<<" set flavour "<<Flavour(iit->first)<<" massless "<<std::endl; 
      }
      else {
	msg_Tracking()<<" set flavour "<<Flavour(iit->first)<<" massive "<<std::endl; 
      }
    }
  }
  //set electrical charges
  iit=cic.begin();
  for (;iit!=cic.end();iit++) {
    if (s_kftable.find(iit->first)!=s_kftable.end()) {
      s_kftable[iit->first]->m_icharge = iit->second;
      msg_Tracking()<<" set charge of "<<Flavour(iit->first)<<" to "
		    <<Flavour(iit->first).Charge()<<std::endl; 
    }
  }
  //set strong charges
  iit=csc.begin();
  for (;iit!=csc.end();iit++) {
    if (s_kftable.find(iit->first)!=s_kftable.end()) {
      s_kftable[iit->first]->m_strong = iit->second;
      msg_Tracking()<<" set strong charge of "<<Flavour(iit->first)<<" to "
		    <<Flavour(iit->first).StrongCharge()<<std::endl; 
    }
  }
  // set yukawas (i.e. overwrite the one from mass, if requested by user)
  dit=cdy.begin();
  for (;dit!=cdy.end();dit++) {
    if (s_kftable.find(dit->first)!=s_kftable.end()) {
      s_kftable[dit->first]->m_yuk = dit->second;
      msg_Tracking()<<" set yukawa of "<<Flavour(dit->first)<<" to "<<dit->second<<" GeV"<<std::endl; 
    }
  }
  // set sorting priority
  iit=cip.begin();
  for (;iit!=cip.end();iit++) {
    if (s_kftable.find(iit->first)!=s_kftable.end()) {
      s_kftable[iit->first]->m_priority = iit->second;
      msg_Tracking()<<" set priority of "<<Flavour(iit->first)<<" to "
		    <<Flavour(iit->first).Priority()<<std::endl; 
    }
  }
}

void Model_Base::AddStandardContainers()
{
  s_kftable[kf_resummed] = new
    Particle_Info(kf_resummed,0.,0.,0,1,2,1,1,1,0,"r","r","r","r",0,1);
  s_kftable[kf_jet] = new
    Particle_Info(kf_jet,0.,0.,0,1, 2,1,1,1,0,"j","j","j","j",1,1);
  s_kftable[kf_quark] = new
    Particle_Info(kf_quark,0.,0.,0,1,1,0,1,1,0,"Q","Q","Q","Q",1,1);
  s_kftable[kf_lepton] = new
    Particle_Info(kf_lepton,0.,0.,-3,0,1,0,1,1,0,"l","l","l","l",1,1);
  s_kftable[kf_neutrino] = new
    Particle_Info(kf_neutrino,0.,0.,0,0, 1,0,1,1,0,"v","v","v","v",1,1);
  s_kftable[kf_lepton]->m_priority=2;
  s_kftable[kf_neutrino]->m_priority=1;
  s_kftable[kf_resummed]->Clear();
  s_kftable[kf_jet]->Clear();
  s_kftable[kf_quark]->Clear();
  s_kftable[kf_lepton]->Clear();
  s_kftable[kf_neutrino]->Clear();
  double jet_mass_threshold=p_dataread->GetValue<double>("JET_MASS_THRESHOLD", 10.0);
  for (int i=1;i<7;i++) {
    Flavour addit((kf_code)i);
    if ((addit.Mass()==0.0 || !addit.IsMassive()) && addit.IsOn()) {
      if (addit.Mass(true)<=jet_mass_threshold) {
        s_kftable[kf_jet]->Add(addit);
        s_kftable[kf_jet]->Add(addit.Bar());
        s_kftable[kf_quark]->Add(addit);
        s_kftable[kf_quark]->Add(addit.Bar());
      }
      else {
        msg_Info()<<"Ignoring "<<addit<<" due to JET_MASS_THRESHOLD.\n";
      }
    }
  }
  s_kftable[kf_jet]->Add(Flavour(kf_gluon));
  s_kftable[kf_jet]->SetResummed();
  for (int i=11;i<17;i+=2) {
    Flavour addit((kf_code)i);
    if ((addit.Mass()==0.0 || !addit.IsMassive()) && addit.IsOn()) {
      s_kftable[kf_lepton]->Add(addit);
      s_kftable[kf_lepton]->Add(addit.Bar());
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
  std::vector<std::vector<std::string> > helpsvv;
  if (!p_dataread->MatrixFromFile(helpsvv,"PARTICLE_CONTAINER")) return;
  for (size_t i(0);i<helpsvv.size();++i) {
    if (helpsvv[i].size()<3) continue;
    std::string props, kfs(helpsvv[i][0]);
    size_t opos(kfs.find('['));
    if (opos<std::string::npos) {
      props=kfs.substr(opos+1,kfs.find(']')-opos-1);
      kfs=kfs.substr(0,opos);
    }
    long int nkf(ToType<long int>(kfs));
    if (s_kftable.find(nkf)!=s_kftable.end())
      THROW(critical_error,"Particle ID "+helpsvv[i][0]+" already exists.");
    msg_Debugging()<<helpsvv[i][1]<<" ("<<nkf<<") ["<<props<<"] = {";
    Data_Reader ppread("|",";","#",":");
    ppread.SetAddCommandLine(false);
    ppread.SetString(props);
    std::string name(helpsvv[i][1]), barname(name);
    if (name.find('|')!=std::string::npos) {
      name=name.substr(0,name.find('|'));
      barname=barname.substr(barname.find('|')+1);
    }
    if (barname==name && ppread.StringValue<int>("M",0)==0) barname+="b";
    s_kftable[nkf] = new Particle_Info
      (nkf,ppread.StringValue<double>("m",0.0),//Mass
       ppread.StringValue<double>("W",0.0),//Width
       ppread.StringValue<int>("C",0),//ICharge
       ppread.StringValue<int>("Q",0),//Strong
       ppread.StringValue<int>("S",0),//Spin
       ppread.StringValue<int>("M",0),//Majorana
       1,1,0,name,barname,name,barname);
    s_kftable[nkf]->m_priority=ppread.StringValue<int>("P",0);
    s_kftable[nkf]->Clear();
    for (size_t j(2);j<helpsvv[i].size();++j) {
      msg_Debugging()<<" "<<helpsvv[i][j];
      long int kfc(ToType<long int>(helpsvv[i][j]));
      s_kftable[nkf]->Add(Flavour((kf_code)std::abs(kfc),kfc<0));
      if (s_kftable[std::abs(kfc)]->m_priority)
	msg_Error()<<METHOD<<"(): Changing "<<Flavour(kfc)<<" sort priority: "
		   <<s_kftable[std::abs(kfc)]->m_priority<<" -> "
		   <<s_kftable[nkf]->m_priority<<std::endl;
      s_kftable[std::abs(kfc)]->m_priority=s_kftable[std::abs(nkf)]->m_priority;
    }
    s_kftable[nkf]->SetIsGroup(true);
    msg_Debugging()<<" }\n";
  }
}

void Model_Base::InitializeInteractionModel()
{
  InitVertices();
  for (std::vector<Single_Vertex>::iterator
	 vit(m_v.begin());vit!=m_v.end();) {
    for (size_t i(0);i<vit->cpl.size();++i)
      if (vit->cpl[i].Value().real()==0.0 &&
	  vit->cpl[i].Value().imag()==0.0) {
	vit->cpl.erase(vit->cpl.begin()+i);
	vit->Color.erase(vit->Color.begin()+i);
	vit->Lorentz.erase(vit->Lorentz.begin()+i);
      }
    if (vit->cpl.empty()) vit=m_v.erase(vit);
    else ++vit;
  }
  m_ov=m_v;
  RotateVertices();
  InitMEInfo();
}

int Model_Base::ScalarNumber(const std::string _name) {
  if (p_numbers->count(_name)>0) return (*p_numbers)[_name];
  THROW(fatal_error, "Key "+_name+" not found");
}


double Model_Base::ScalarConstant(const std::string _name) {
  if (p_constants->count(_name)>0) return (*p_constants)[_name];
  THROW(fatal_error, "Key "+_name+" not found");
}


Complex Model_Base::ComplexConstant(const std::string _name) {
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

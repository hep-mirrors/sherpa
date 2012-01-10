#include "MODEL/Main/Model_Base.H"

#include "MODEL/Main/Spectrum_Generator_Base.H"
#include "MODEL/Interaction_Models/Interaction_Model_Base.H"
#include "ATOOLS/Org/Data_Reader.H"
#include "ATOOLS/Org/Run_Parameter.H"
#include "ATOOLS/Org/Message.H"
#include "MODEL/Interaction_Models/Vertex.H"
#include "ATOOLS/Org/Exception.H"

using namespace MODEL;
using namespace ATOOLS;

namespace MODEL 
{
  Model_Base *s_model;
}

Model_Base::Model_Base(std::string _dir,std::string _file,bool _elementary) :
  p_model(NULL), m_dir(_dir), m_file(_file), m_elementary(_elementary), 
  p_dataread(NULL), p_numbers(NULL), p_constants(NULL), p_complexconstants(NULL), 
  p_functions(NULL), p_matrices(NULL), p_spectrumgenerator(NULL), p_vertex(NULL), 
  p_vertextable(NULL), m_vinfo(0)
{
  p_dataread = new Data_Reader(" ",";","!","=");
  p_dataread->AddComment("#");
  p_dataread->AddWordSeparator("\t");
  p_dataread->SetInputPath(m_dir);
  p_dataread->SetInputFile(m_file);
}

Model_Base::~Model_Base() 
{
  delete p_model;
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
  if (p_matrices!=NULL)          delete p_matrices;
  if (p_dataread!=NULL)          delete p_dataread;
  if (p_spectrumgenerator!=NULL) delete p_spectrumgenerator;
  if (p_vertex!=NULL)            delete p_vertex;
  if (p_vertextable!=NULL)       delete p_vertextable;
}

void Model_Base::GetCouplings(Coupling_Map &cpls) const
{
  DEBUG_FUNC(&cpls);
  for (ScalarFunctionsMap::const_iterator
	 cit(p_functions->begin());cit!=p_functions->end();++cit) {
    std::string tag(cit->second->Name());
    cpls.insert(std::make_pair(tag,new Coupling_Data
      (cit->second,tag,p_model->ScalarFunction(cit->first,rpa->gen.CplScale()))));
    msg_Debugging()<<"  '"<<tag<<"' -> ("<<cpls.lower_bound(tag)->second<<")"
		   <<*cpls.lower_bound(tag)->second<<"\n";
  }
}

void Model_Base::ShowSyntax(const size_t i)
{
  if (!msg_LevelIsInfo() || i==0) return;
  msg_Out()<<METHOD<<"(): {\n\n"
	   <<"   // available model implementations (specified through MODEL=<value>)\n\n";
  Model_Getter_Function::PrintGetterInfo(msg->Out(),25);
  msg_Out()<<"\n   // available sets of interaction vertices (specified through SIGNAL_MODEL=<value> in ME.dat)\n"
	   <<"   // default given by MODEL switch\n\n";
  Interaction_Model_Base::Interaction_Model_Getter_Function::
    PrintGetterInfo(msg->Out(),25);
  msg_Out()<<"\n}"<<std::endl;
}

void Model_Base::ReadParticleData() {
  std::map<int,double> cdm, cdw;
  std::map<int,int> cia, cis, cim, cic;
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

  //set masses
  std::map<int,double>::const_iterator dit=cdm.begin();
  for (;dit!=cdm.end();dit++) {
    if (s_kftable.find(dit->first)!=s_kftable.end()) {
      s_kftable[dit->first]->m_mass = dit->second;
      s_kftable[dit->first]->m_yuk = dit->second;
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
    s_kftable[nkf] = new Particle_Info
      (nkf,ppread.StringValue<double>("m",0.0),//Mass
       ppread.StringValue<double>("W",0.0),//Width
       ppread.StringValue<int>("C",0),//ICharge
       ppread.StringValue<int>("I",0),//Isoweak
       ppread.StringValue<int>("Q",0),//Strong
       ppread.StringValue<int>("S",0),//Spin
       ppread.StringValue<int>("M",0),//Majorana
       1,1,0,helpsvv[i][1],helpsvv[i][1]);
    s_kftable[nkf]->Clear();
    for (size_t j(2);j<helpsvv[i].size();++j) {
      msg_Debugging()<<" "<<helpsvv[i][j];
      long int kfc(ToType<long int>(helpsvv[i][j]));
      s_kftable[nkf]->Add(Flavour((kf_code)abs(kfc),kfc<0));
    }
    msg_Debugging()<<" }\n";
  }
}

void Model_Base::InitializeInteractionModel()
{
  Data_Reader read(" ",";","!","=");
  read.AddComment("#");
  read.AddWordSeparator("\t");
  read.SetInputPath(m_dir);
  read.SetInputFile(rpa->gen.Variable("ME_DATA_FILE"));
  std::string modeltype   = read.GetValue<std::string>("SIGNAL_MODEL",m_name);
  std::string cplscheme   = read.GetValue<std::string>("COUPLING_SCHEME","Running_alpha_S");
  std::string massscheme  = read.GetValue<std::string>("YUKAWA_MASSES","Running");
  std::string widthscheme = read.GetValue<std::string>("WIDTH_SCHEME","Fixed");
  
  p_model = Interaction_Model_Base::Interaction_Model_Getter_Function::GetObject
    (modeltype,Interaction_Model_Arguments(this,cplscheme,massscheme));
  
  if (p_model==NULL) THROW(not_implemented,"Interaction model not implemented");

  p_vertex        = new Vertex(p_model);
  p_vertextable   = new Vertex_Table;
  for (int i=0;i<p_vertex->MaxNumber();++i) {
    if ((*p_vertex)[i]->on) {
      (*p_vertextable)[(*p_vertex)[i]->in[0]].push_back((*p_vertex)[i]);
    }
  }
  InitMEInfo();
}

int Model_Base::ScalarNumber(const std::string _name) {
  if (p_numbers->empty()) {
    msg_Error()<<"Error in Model_Base::ScalarNumber("<<_name<<") : "<<std::endl
	       <<"   No numbers stored in model "<<m_name<<". Return 0."<<std::endl;
    return 0;
  }
  if (p_numbers->count(_name)>0) return (*p_numbers)[_name];

  msg_Error()<<"Error in Model_Base::ScalarNumber("<<_name<<") : "<<std::endl
	     <<"   Key not found in model "<<m_name<<". Return 0."<<std::endl;
  return 0;
}


double Model_Base::ScalarConstant(const std::string _name) {
  if (p_constants->empty()) {
    msg_Error()<<"Error in Model_Base::ScalarConstant("<<_name<<") : "<<std::endl
	       <<"   No constants stored in model "<<m_name<<". Return 0."<<std::endl;
    return 0.;
  }
  if (p_constants->count(_name)>0) return (*p_constants)[_name];

  msg_Error()<<"Error in Model_Base::ScalarConstant("<<_name<<") : "<<std::endl
	     <<"   Key not found in model "<<m_name<<". Return 0."<<std::endl;
  return 0.;
}


Complex Model_Base::ComplexConstant(const std::string _name) {
  if (p_complexconstants->empty()) {
    msg_Error()<<"Error in Model_Base::ComplexConstant("<<_name<<") : "<<std::endl
	       <<"   No constants stored in model "<<m_name<<". Return 0."<<std::endl;
    return 0.;
  }
  if (p_complexconstants->count(_name)>0) return (*p_complexconstants)[_name];

  msg_Error()<<"Error in Model_Base::ComplexConstant("<<_name<<") : "<<std::endl
	     <<"   Key not found in model "<<m_name<<". Return 0."<<std::endl;
  return Complex(0.,0.);
}


Function_Base * Model_Base::GetScalarFunction(const std::string _name) {
  if (p_functions->empty()) {
    msg_Error()<<"Error in Model_Base::ScalarFunction("<<_name<<") : "<<std::endl
	       <<"   No functions stored in model "<<m_name<<". Return 0."<<std::endl;
    return NULL;
  }
  if (p_functions->count(_name)>0) return (*p_functions)[_name];

  msg_Error()<<"Error in Model_Base::ScalarFunction("<<_name<<") : "<<std::endl
	     <<"   Key not found in model "<<m_name<<". Return 0."<<std::endl;
  return NULL;
}


double Model_Base::ScalarFunction(const std::string _name,double _t) {
  if (p_functions->empty()) {
    msg_Error()<<"Error in Model_Base::ScalarNumber("<<_name<<") : "<<std::endl
	       <<"   No functions stored in model "<<m_name<<". Return 0."<<std::endl;
    return 0.;
  }
  if (p_functions->count(_name)>0) return (*(*p_functions)[_name])(_t);

  msg_Error()<<"Error in Model_Base::ScalarNumber("<<_name<<") : "<<std::endl
	     <<"   Key not found in model "<<m_name<<". Return 0."<<std::endl;
  return 0.;
}


double Model_Base::ScalarFunction(const std::string _name) {
  if (p_functions->empty()) {
    msg_Error()<<"Error in Model_Base::ScalarNumber("<<_name<<") : "<<std::endl
	       <<"   No functions stored in model "<<m_name<<". Return 0."<<std::endl;
    return 0.;
  }
  if (p_functions->count(_name)>0) return (*(*p_functions)[_name])();

  msg_Error()<<"Error in Model_Base::ScalarNumber("<<_name<<") : "<<std::endl
	     <<"   Key not found in model "<<m_name<<". Return 0."<<std::endl;
  return 0.;
}


CMatrix Model_Base::ComplexMatrix(const std::string _name) {
  if (p_matrices->empty()) {
    msg_Error()<<"Error in Model_Base::ComplexMatrix("<<_name<<") : "<<std::endl
	       <<"   No matrices stored in model "<<m_name<<". Return 0."<<std::endl;
    return CMatrix(1);
  }
  if (p_matrices->count(_name)>0) return (*p_matrices)[_name];

  msg_Error()<<"Error in Model_Base::ComplexMatrix("<<_name<<") : "<<std::endl
	     <<"   Key not found in model "<<m_name<<". Return 0."<<std::endl;
  return CMatrix(1);
}


Complex Model_Base::ComplexMatrixElement(const std::string _name,const int _i,const int _j) {
  if (p_matrices->empty()) {
    msg_Error()<<"Error in Model_Base::ComplexMatrixElement("<<_name<<")("<<_i<<","<<_j<<") : "<<std::endl
	       <<"   No matrices stored in model "<<m_name<<". Return 0."<<std::endl;
    return 0;
  }
  if (p_matrices->count(_name)>0) {
    int rank = (*p_matrices)[_name].Rank();
    if (_i<rank && _j<rank && 0<=_i && 0<=_j) return (*p_matrices)[_name][_i][_j];
  }

  msg_Error()<<"Error in Model_Base::ComplexMatrixElement("<<_name<<")("<<_i<<","<<_j<<") : "<<std::endl
	     <<"   Key not found in model "<<m_name<<". Return 0."<<std::endl;
  return 0;
}

bool Model_Base::CheckFlavours(int nin, int nout, Flavour* flavs)
{
  return true;
}

void Model_Base::InitMEInfo()
{
  msg_Debugging()<<METHOD<<"(): {\n";
  m_fls.clear();
  bool hasndec(false);
  std::set<Flavour> fls;
  msg_Debugging()<<"\n  add vertices\n\n";
  std::vector<Single_Vertex> &all(p_vertex->Vertices());
  for (size_t i=0;i<all.size();++i) {
    if (all[i].on) {
      m_vmap.insert(VMap_Key(all[i].PID(),&all[i]));
      m_vtable[all[i].in[0]].push_back(&all[i]);
      if (all[i].nleg>3) {
	if (all[i].dec<0) m_vinfo|=2;
	else hasndec=true;
      }
      for (int j(0);j<all[i].nleg;++j) fls.insert(all[i].in[j]);
      if (msg_LevelIsDebugging()) {
	msg_Debugging()
	  <<"  "<<all[i].PID()<<" ("<<all[i].oew
	  <<","<<all[i].oqcd<<") "<<(all[i].dec>0?'{':(all[i].dec<0?'(':'['))
	  <<all[i].Lorentz.front()->Type()<<","<<all[i].Color[0].PID();
	for (size_t j(1);j<all[i].Lorentz.size();++j)
	  msg_Out()<<"|"<<all[i].Lorentz[j]->Type()
		   <<","<<all[i].Color[j].PID();
	msg_Out()<<(all[i].dec>0?'}':(all[i].dec<0?')':']'))<<"\n";
      }
    }
  }
  if (hasndec) m_vinfo|=1;
  msg_Debugging()<<"\n  add particles\n\n";
  for (std::set<Flavour>::const_iterator 
	 fit(fls.begin());fit!=fls.end();++fit) {
      m_fls.push_back(*fit);
      msg_Debugging()<<"  "<<*fit<<"\n";
  }
  msg_Debugging()<<"\n}\n";
}

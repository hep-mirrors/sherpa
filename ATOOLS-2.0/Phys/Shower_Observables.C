#include "Shower_Observables.H"
#include "MyStrStream.H"
#include "MathTools.H"
#include "Run_Parameter.H"
#include "Vector.H"
#include "Primitive_Analysis.H"

#include <algorithm>

#ifdef PROFILE__Analysis_Phase
#include "prof.hh"
#else 
#define PROFILE_HERE {}
#define PROFILE_LOCAL(LOCALNAME) {}
#endif

using namespace ATOOLS;
using namespace std;


bool ParticleIsInList(const Particle * const p,  const Particle_List & pl) 
{
  bool hit=0;
  for (size_t j=0;j<pl.size();++j) {
    if (p==pl[j]) {
      hit=1;
      break;
    }
  }
  return hit;
}




void Shower_Observables::InitObservables() {
  all_obs.flav     = Flavour(kf::none);
  all_obs.jet_ini  = 0;
  all_obs.jetrates = new Jetrates(11,1.e-6,1.,80,0);
  all_obs.multi    = new Multiplicity(00,-0.5,50.5,51,0,"FinalState");
  all_obs.wz_pt    = new PT_Distribution(00,0.,200.,100,1,Flavour(kf::W));
  all_obs.jet_pt   = new PT_Distribution(00,0.,200.,100,6,Flavour(kf::jet));
  all_obs.sum      =0.;
}


Shower_Observables::Shower_Observables(int _type,double _xmin,double _xmax,int _nbins,
      Selector_Base * _sel,int _njet,int _nflavs,int _dohad) 
{
  p_histo=0;
  m_splitt_flag=false;
  m_type = _type; m_xmin = _xmin; m_xmax = _xmax; m_nbins = _nbins; p_sel = _sel;
  m_name  = std::string("shower_obs");
  InitObservables();
}

void Shower_Observables::SetAnalysis(Primitive_Analysis* ana)
{
  p_ana=ana;
  all_obs.jetrates->SetAnalysis(ana); 
  all_obs.multi->SetAnalysis(ana); 
  all_obs.wz_pt->SetAnalysis(ana); 
  all_obs.jet_pt->SetAnalysis(ana);
}

void Shower_Observables::Evaluate(const Blob_List & blobs ,double value, int ncount) {
  bool do_fl=0,do_jet=1;


  Particle_List pl;

  int njet_ini=0;
  Flavour lfl;
  for (Blob_Const_Iterator blit=blobs.begin();blit!=blobs.end();++blit) {
    if ((*blit)->Type()[0]=='S') {
      njet_ini=(*blit)->NOutP();
      lfl =(*blit)->OutParticle(0)->Flav();
      for (int i=0;i<(*blit)->NInP();++i) {
	Particle * p =(*blit)->InParticle(i);
	if (!ParticleIsInList(p,pl)) pl.push_back(p);
      }
    }
  }

  
  Particle_List * pl_fs = p_ana->GetParticleList("FinalState");
  if (pl_fs!=0) {
    copy(pl_fs->begin(),pl_fs->end(),back_inserter(pl));
  }
  else {
    std::cout<<" ERROR in Shower_Observables::Evaluate : no particle list recieved "<<std::endl;
  }

  all_obs.jetrates->Evaluate(pl,value,ncount);
  all_obs.multi->Evaluate(pl,value,ncount);
  all_obs.wz_pt->Evaluate(pl,value,ncount);
  all_obs.jet_pt->Evaluate(pl,value,ncount);
  all_obs.sum+=value;

  size_t nc=0;
  if (do_fl) {
  nc=fl_obs.size();
  for (size_t i=0;i<nc;++i) {
    if (fl_obs[i].flav==lfl) 
      nc=i;
  }
  if (nc==fl_obs.size()) {
    Event_Obi obs;
    obs.flav     = lfl;  
    obs.jet_ini  = 0;
    obs.jetrates = new Jetrates(all_obs.jetrates,lfl.Name()+std::string("_"));
    obs.multi    = new Multiplicity(all_obs.multi,lfl.Name()+std::string("_"));
    obs.wz_pt    = new PT_Distribution(all_obs.wz_pt,lfl.Name()+std::string("_"));
    obs.jet_pt   = new PT_Distribution(all_obs.jet_pt,lfl.Name()+std::string("_"));
    obs.jetrates->SetAnalysis(p_ana); 
    obs.multi->SetAnalysis(p_ana); 
    obs.wz_pt->SetAnalysis(p_ana); 
    obs.jet_pt->SetAnalysis(p_ana);
    obs.sum      = 0.;
    fl_obs.push_back(obs);
  }
  fl_obs[nc].jetrates->Evaluate(pl,value,ncount);
  fl_obs[nc].multi->Evaluate(pl,value,ncount);
  fl_obs[nc].wz_pt->Evaluate(pl,value,ncount);
  fl_obs[nc].jet_pt->Evaluate(pl,value,ncount);
  fl_obs[nc].sum+=value;
  }


  if (do_jet) {
  // fill histograms (for each number of seed)
  nc=jet_obs.size();
  for (size_t i=0;i<nc;++i) {
    if (jet_obs[i].jet_ini==njet_ini) 
      nc=i;
  }
  if (nc==jet_obs.size()) {
    std::string jname;
    MyStrStream sstr;
    sstr<<"j";
    sstr<<njet_ini;
    sstr<<"_";
    sstr>>jname;

    Event_Obi obs;
    obs.flav     = Flavour(kf::none);  
    obs.jet_ini  = njet_ini;
    obs.jetrates = new Jetrates(all_obs.jetrates,jname);
    obs.multi    = new Multiplicity(all_obs.multi,jname);
    obs.wz_pt    = new PT_Distribution(all_obs.wz_pt,jname);
    obs.jet_pt   = new PT_Distribution(all_obs.jet_pt,jname);
    obs.jetrates->SetAnalysis(p_ana); 
    obs.multi->SetAnalysis(p_ana); 
    obs.wz_pt->SetAnalysis(p_ana); 
    obs.jet_pt->SetAnalysis(p_ana);
    obs.sum      = 0.;
    jet_obs.push_back(obs);
  }
  jet_obs[nc].jetrates->Evaluate(pl,value,ncount);
  jet_obs[nc].multi->Evaluate(pl,value,ncount);
  jet_obs[nc].wz_pt->Evaluate(pl,value,ncount);
  jet_obs[nc].jet_pt->Evaluate(pl,value,ncount);
  jet_obs[nc].sum+=value;
  }

  if (do_fl & do_jet) {
  // fill histograms (for each number of seeds AND each initial flavour)
  nc=fl_jet_obs.size();
  for (size_t i=0;i<nc;++i) {
    if ((fl_jet_obs[i].jet_ini==njet_ini) && (fl_jet_obs[i].flav==lfl)) 
      nc=i;
  }
  if (nc==fl_jet_obs.size()) {
    std::string jname;
    MyStrStream sstr;
    sstr<<"j"<<njet_ini<<"_"<<lfl.Name()<<"_";
    sstr>>jname;

    Event_Obi obs;
    obs.flav     = lfl;  
    obs.jet_ini  = njet_ini;
    obs.jetrates = new Jetrates(all_obs.jetrates,jname);
    obs.multi    = new Multiplicity(all_obs.multi,jname);
    obs.wz_pt    = new PT_Distribution(all_obs.wz_pt,lfl.Name()+std::string("_"));
    obs.jet_pt   = new PT_Distribution(all_obs.jet_pt,lfl.Name()+std::string("_"));
    obs.jetrates->SetAnalysis(p_ana); 
    obs.multi->SetAnalysis(p_ana); 
    obs.wz_pt->SetAnalysis(p_ana); 
    obs.jet_pt->SetAnalysis(p_ana);
    obs.sum      = 0.;
    fl_jet_obs.push_back(obs);
  }
  fl_jet_obs[nc].jetrates->Evaluate(pl,value,ncount);
  fl_jet_obs[nc].multi->Evaluate(pl,value,ncount);
  fl_jet_obs[nc].wz_pt->Evaluate(pl,value,ncount);
  fl_jet_obs[nc].jet_pt->Evaluate(pl,value,ncount);
  fl_jet_obs[nc].sum+=value;
  }
}


void Shower_Observables::EndEvaluation() {
  all_obs.jetrates->EndEvaluation(); 
  all_obs.multi->EndEvaluation();
  all_obs.wz_pt->EndEvaluation(); 
  all_obs.jet_pt->EndEvaluation(); 
  for (size_t i=0;i<fl_obs.size();++i) {
    double scale =fl_obs[i].sum/all_obs.sum;
    fl_obs[i].jetrates->EndEvaluation(scale); 
    fl_obs[i].multi->EndEvaluation(scale);
    fl_obs[i].wz_pt->EndEvaluation(scale); 
    fl_obs[i].jet_pt->EndEvaluation(scale); 
  }
  for (size_t i=0;i<jet_obs.size();++i) {
    double scale =jet_obs[i].sum/all_obs.sum;
    jet_obs[i].jetrates->EndEvaluation(scale); 
    jet_obs[i].multi->EndEvaluation(scale);
    jet_obs[i].wz_pt->EndEvaluation(scale); 
    jet_obs[i].jet_pt->EndEvaluation(scale); 
  }
  for (size_t i=0;i<fl_jet_obs.size();++i) {
    double scale =fl_jet_obs[i].sum/all_obs.sum;
    fl_jet_obs[i].jetrates->EndEvaluation(scale); 
    fl_jet_obs[i].multi->EndEvaluation(scale);
    fl_jet_obs[i].wz_pt->EndEvaluation(scale); 
    fl_jet_obs[i].jet_pt->EndEvaluation(scale); 
  }
}

void Shower_Observables::Output(const std::string & pname) {
  all_obs.jetrates->Output(pname); 
  all_obs.multi->Output(pname);
  all_obs.wz_pt->Output(pname); 
  all_obs.jet_pt->Output(pname); 
  for (size_t i=0;i<fl_obs.size();++i) {
    fl_obs[i].jetrates->Output(pname); 
    fl_obs[i].multi->Output(pname);
    fl_obs[i].wz_pt->Output(pname); 
    fl_obs[i].jet_pt->Output(pname); 
  }
  for (size_t i=0;i<jet_obs.size();++i) {
    jet_obs[i].jetrates->Output(pname); 
    jet_obs[i].multi->Output(pname);
    jet_obs[i].wz_pt->Output(pname); 
    jet_obs[i].jet_pt->Output(pname); 
  }
  for (size_t i=0;i<fl_jet_obs.size();++i) {
    fl_jet_obs[i].jetrates->Output(pname); 
    fl_jet_obs[i].multi->Output(pname);
    fl_jet_obs[i].wz_pt->Output(pname); 
    fl_jet_obs[i].jet_pt->Output(pname); 
  }
}


//======================================================================

std::ostream & ATOOLS::operator<<(std::ostream & s, const Jetrate_Data * jd)
{
  for (size_t i=0;i<jd->jets.size();++i) s<<jd->jets[i]<<"\t"; s<<std::endl;
  for (size_t i=0;i<jd->ys.size();++i) s<<jd->ys[i]<<"\t"; s<<std::endl;
  return s;
}

template <> Jetrate_Data * Blob_Data_Base::Get<Jetrate_Data*>()
{
  return ((Blob_Data<Jetrate_Data*>*)this)->Get();
}

template <> Blob_Data<Jetrate_Data*>::~Blob_Data() {
  delete m_data;
}

template class Blob_Data<Jetrate_Data*>;

template Jetrate_Data * Blob_Data_Base::Get<Jetrate_Data*>();
//----------------------------------------------------------------------

Jetrates::Jetrates(int _type,double _xmin,double _xmax,int _nbins,
		   Selector_Base * _sel, const std::string & lname)
{
  p_partner = 0;
  m_type = _type; m_xmin = _xmin; m_xmax = _xmax; m_nbins = _nbins; p_sel = _sel; m_listname=lname;
  m_name  = lname+std::string("jetrate");
  //    if (m_listname!="") m_name=m_listname+std::string("_")+m_name;

  p_histo = 0;
  p_jfind = 0; //new Jet_Finder(rpa.gen.Ycut(),1);
  p_ktalg = new Kt_Algorithm();

  int had=0;
  double xmin=m_xmin;
  double xmax=m_xmax;
  if (rpa.gen.Beam1()==Flavour(kf::p_plus) || rpa.gen.Beam1()==Flavour(kf::p_plus).Bar()) {
    had=1;
    xmin=0.01*xmin;
    xmax=0.01*xmax;
  }

  m_jets.push_back(6);
  m_histos.push_back(new Histogram(11,xmin,xmax,m_nbins));
  m_rates.push_back(new Histogram(11,xmin,xmax,m_nbins));
  m_jets.push_back(5);
  m_histos.push_back(new Histogram(11,xmin,xmax,m_nbins));
  m_rates.push_back(new Histogram(11,xmin,xmax,m_nbins));
  m_jets.push_back(4);
  m_histos.push_back(new Histogram(11,xmin,xmax,m_nbins));
  m_rates.push_back(new Histogram(11,xmin,xmax,m_nbins));
  m_jets.push_back(3);
  m_histos.push_back(new Histogram(11,xmin,xmax,m_nbins));
  m_rates.push_back(new Histogram(11,xmin,xmax,m_nbins));
  m_jets.push_back(2);
  m_histos.push_back(new Histogram(11,xmin,xmax,m_nbins));
  m_rates.push_back(new Histogram(11,xmin,xmax,m_nbins));
  if (had) {
    m_jets.push_back(1);
    m_histos.push_back(new Histogram(11,xmin,xmax,m_nbins));
    m_rates.push_back(new Histogram(11,xmin,xmax,m_nbins));
    m_jets.push_back(0);
    m_histos.push_back(new Histogram(11,xmin,xmax,m_nbins));
    m_rates.push_back(new Histogram(11,xmin,xmax,m_nbins));
  }
  m_ymax=-1;
  m_ymin=2;
};


Jetrates::Jetrates(Jetrates * _partner, std::string _prefix)
{
  p_partner = _partner;
  m_type = p_partner->m_type; m_xmin = p_partner->m_xmin; m_xmax = p_partner->m_xmax; m_nbins = p_partner->m_nbins; 
  p_sel = p_partner->p_sel; m_listname = p_partner->m_listname;
  m_name  = _prefix; 
  p_histo = 0;
  p_jfind = 0;

  for (size_t i=0;i<p_partner->m_jets.size();++i) {
    m_jets.push_back(p_partner->m_jets[i]);
    m_histos.push_back(new Histogram(p_partner->m_histos[i]));
    m_histos.back()->Reset();  
    m_rates.push_back(new Histogram(p_partner->m_rates[i]));
    m_rates.back()->Reset();
  }

  m_ymax=-1;
  m_ymin=2;
};


Primitive_Observable_Base * Jetrates::Copy() const {
  return new Jetrates(m_type,m_xmin,m_xmax,m_nbins,p_sel,m_listname);
}

void Jetrates::Evaluate(const Blob_List &,double value, int ncount) 
{
  PROFILE_HERE;
  Particle_List * pl=0;
  if (m_listname=="") pl=p_ana->GetParticleList("FinalState");
  else pl=p_ana->GetParticleList(m_listname);
  Evaluate(*pl,value, ncount);
}


void Jetrates::Evaluate(const Particle_List & pl,double value, int ncount) {
  PROFILE_HERE;
  m_ys.clear();

  if (p_partner==0) {
    std::string key=std::string("KtJetrates")+m_listname;

    if (!(*p_ana)[key]) {
      if (value==0.) {
	for (size_t i=0;i<m_jets.size();++i) m_ys.push_back(-1.);
      }
      else {
	PROFILE_LOCAL("Jetrates::Evaluate-Construct Jets");
	if (p_jfind) p_jfind->ConstructJets(&pl,m_jets,m_ys,true);
	else m_ys.resize(m_jets.size(), -1.);
      }
      Jetrate_Data * jd = new Jetrate_Data;
      jd->jets=m_jets;
      jd->ys  =m_ys;
      p_ana->AddData(key,new Blob_Data<Jetrate_Data*>(jd));
    }
    else {
      Jetrate_Data * jd=(*p_ana)[key]->Get<Jetrate_Data*>();
      m_ys=jd->ys;
    }
    // create jetlist
    if (m_listname=="") {
      Particle_List * pl_jets = p_ana->GetParticleList("KtJets");
      if (pl_jets==0) {
	PROFILE_LOCAL("Jetrates::Evaluate-CreateJetlist");
	if (p_jfind) {
	  pl_jets = new Particle_List();
	  copy(pl.begin(),pl.end(),back_inserter(*pl_jets));
	  p_jfind->ConstructJets(pl_jets,rpa.gen.Ycut(),true);
	  p_ana->AddParticleList("KtJets",pl_jets);
	  
	}
	else if (p_ktalg) {
	  pl_jets = new Particle_List();
	  p_ktalg->ConstructJets(&pl,pl_jets,0,0.49);
	  p_ana->AddParticleList("KtJets",pl_jets);
	}
      }
    }
  }
  else m_ys=p_partner->m_ys;


  for (size_t k=0; k<m_jets.size(); ++k ){
    m_histos[k]->Insert(m_ys[k], value, ncount);
  }

  for (size_t k=0; k<m_jets.size()-1; ++k ){
    m_rates[k]->InsertRange(m_ys[k],m_ys[k+1], value);
  }
  m_rates.back()->InsertRange(m_ys.back(),1., value);

  m_ymax=Max(m_ymax,m_ys.back());
  m_ymin=m_ys.back();
}

void  Jetrates::EndEvaluation() {
  PROFILE_HERE;
  for (size_t i=0; i<m_histos.size();++i) {
    m_histos[i]->Finalize();
    m_histos[i]->Output();
    m_rates[i]->Finalize();
    m_rates[i]->Output();
  }
  msg.Debugging()<<std::endl<<std::endl;
}

void  Jetrates::EndEvaluation(double scale) {
  PROFILE_HERE;
  for (size_t i=0; i<m_histos.size();++i) {
    m_histos[i]->Finalize();
    m_histos[i]->Scale(scale);
    m_histos[i]->Output();
    m_rates[i]->Finalize();
    m_rates[i]->Scale(scale);
    m_rates[i]->Output();
  }
  msg.Debugging()<<std::endl<<std::endl;
}

void Jetrates::Output(const std::string & pname) {
  PROFILE_HERE;
  int  mode_dir = 448;
  mkdir((pname).c_str(),mode_dir); 
  for (size_t i=0; i<m_histos.size();++i) {
    std::string fname;
    MyStrStream s;
    s<<m_jets[i];
    s<<".dat";
    s>>fname;
    //    if (m_listname!="") m_name=m_listname+std::string("_")+m_name;
    if (m_name==("jetrate")) {
      m_histos[i]->Output((pname+std::string("/diffrate")+fname).c_str());
      m_rates[i]->Output((pname+std::string("/totrate")+fname).c_str());
    } 
    else{
      m_histos[i]->Output((pname+std::string("/")+m_listname+std::string("_diffrate")+fname).c_str());
      m_rates[i]->Output((pname+std::string("/")+m_listname+std::string("_totrate")+fname).c_str());
    }
  }
}

Primitive_Observable_Base & Jetrates::operator+=(const Primitive_Observable_Base & ob) {
  PROFILE_HERE;
  if (m_xmin!=ob.Xmin() || m_xmax!=ob.Xmax() || m_nbins!=ob.Nbins()) {
    std::cout<<m_xmin<<"\t"<<m_xmax<<"\t"<<m_nbins<<std::endl;
    std::cout<<ob.Xmin()<<"\t"<<ob.Xmax()<<"\t"<<ob.Nbins()<<std::endl;
    std::cout<<" ERROR: in Jetrates::operator+="<<std::endl;
    return *this;
  }
  Jetrates * jr = ((Jetrates*)(&ob));

  if (m_histos.size()==jr->m_histos.size()) {
    for (size_t i=0; i<m_histos.size();++i) {
      (*m_histos[i])+=(*jr->m_histos[i]);
      (*m_rates[i]) +=(*jr->m_rates[i]);
    }
  }
  else {
    std::cout<<" ERROR: in Jetrates::operator+="<<std::endl;
    std::cout<<"        histo size not equal!"<<std::endl;
  }
  return *this;
}

void Jetrates::Reset()
{
  for (size_t i=0; i<m_histos.size();++i) {
    m_histos[i]->Reset();
    m_rates[i]->Reset();
  }
}

Jetrates::~Jetrates() {
  if (p_jfind) delete p_jfind;
  if (p_ktalg) delete p_ktalg;
}

//----------------------------------------------------------------------
Multiplicity::Multiplicity(int _type,double _xmin,double _xmax,int _nbins, int mode, const std::string & lname)
{
  m_type = _type; m_xmin = _xmin; m_xmax = _xmax; m_nbins = _nbins; p_sel = 0; m_mode=mode; m_listname=lname;
  switch (m_mode) {
  case 0:  
    if (lname=="FinalState")
      m_name = std::string("multi.dat");
    else
      m_name  = lname+std::string("_multi.dat");
    break;
  case 2:
    if (lname=="FinalState")
      m_name = std::string("multi_charged_hadrons.dat");
    else
    m_name  = lname+std::string("_multi_charged_hadrons.dat");
    break;
  default:
    if (lname=="FinalState")
      m_name = std::string("multi_total.dat");
    else
      m_name  = lname+std::string("_multi_total.dat");
  }
  p_histo = new Histogram(m_type,m_xmin,m_xmax,m_nbins);
};

Multiplicity::Multiplicity(Multiplicity * old) {
  m_type  = old->Type();
  p_sel   = old->Sel();
  m_xmin  = old->Xmin();
  m_xmax  = old->Xmax();
  m_nbins = old->Nbins();
  m_name  = old->Name();
  m_listname = old->m_listname;
  p_histo = new Histogram(old->Histo());
  p_histo -> Reset();
}


Multiplicity::Multiplicity(Multiplicity * _partner, std::string _prefix)
{
  p_partner = _partner;
  m_type = p_partner->m_type; m_xmin = p_partner->m_xmin; m_xmax = p_partner->m_xmax; m_nbins = p_partner->m_nbins; 
  p_sel = p_partner->p_sel; m_listname = p_partner->m_listname;
  m_name  = _prefix+std::string("multi.dat");
  p_histo = new Histogram(m_type,m_xmin,m_xmax,m_nbins);
}


void Multiplicity::Evaluate(const Particle_List & pl,double value, int ncount) {
  PROFILE_HERE;
  int multi=0;
  if (value!=0.) {
    multi=pl.size();
    std::string key=m_listname+std::string("Multiplicity");
    if (!(*p_ana)[key]) {
      p_ana->AddData(key,new Blob_Data<int>(multi));
    }
  }
  if (m_mode==0) p_histo->Insert(multi,value,ncount);
  else p_histo->Insert(pl.size()-2,value,ncount);
}

void  Multiplicity::Evaluate(const Blob_List & blobs,double value, int ncount) {
  PROFILE_HERE;
  Particle_List * pl=p_ana->GetParticleList(m_listname);
  Evaluate(*pl,value, ncount);
}


Primitive_Observable_Base * Multiplicity::Copy() const{
  return new Multiplicity(m_type,m_xmin,m_xmax,m_nbins,m_mode,m_listname);
}

//----------------------------------------------------------------------




ME_Rate::ME_Rate(int _type,double _xmin,double _xmax,int _nbins,std::string prefix)
{
  m_type = _type; m_xmin = _xmin; m_xmax = _xmax; m_nbins = _nbins; p_sel = 0;
  m_name  = prefix+std::string("_rates.dat");
  p_histo = new Histogram(m_type,m_xmin,m_xmax,m_nbins);
  m_sum =0.;
}

ME_Rate::ME_Rate(ME_Rate * old) {
  m_type  = old->Type();
  p_sel   = old->Sel();
  m_xmin  = old->Xmin();
  m_xmax  = old->Xmax();
  m_nbins = old->Nbins();
  m_name  = old->Name();
  p_histo = new Histogram(old->Histo());
  p_histo -> Reset();
}


void ME_Rate::Evaluate(const Blob_List & blobs,double value, int ncount) {
  for (Blob_Const_Iterator bit=blobs.begin();bit!=blobs.end();++bit) {
    if ((*bit)->Type().find("Signal") !=std::string::npos ) {
      p_histo->Insert((*bit)->NOutP(),value);
      m_sum+=value;

      int jets=(*bit)->NOutP();
      Flavour flavs[3];
      flavs[0]=(*bit)->OutParticle(0)->Flav();
      if (jets>=3)       flavs[1]=(*bit)->OutParticle(2)->Flav();
      if (jets>=5)       flavs[2]=(*bit)->OutParticle(4)->Flav();
      size_t nc=m_all_rates.size();
      for (size_t i=0;i<nc;++i) {
	if (m_all_rates[i].flavs[0]==flavs[0]
	    && m_all_rates[i].flavs[1]==flavs[1] && m_all_rates[i].flavs[2]==flavs[2] ) {
	  nc=i;
	}
      }

      if (nc==m_all_rates.size()) {
	ME_Data  me_data;
	me_data.jets     = jets;
	me_data.sum      = 0.;
	me_data.flavs[0] = flavs[0];
	me_data.flavs[1] = flavs[1];
	me_data.flavs[2] = flavs[2];
	if (jets<=2) {
	  me_data.name     = flavs[0].Name()+m_name;
	}
	else if (jets<=4) {
	  me_data.name     = string(flavs[0].Name())
	                       +string(flavs[1].Name())+m_name;
	}
	else {
	  me_data.name     = string(flavs[0].Name())
	    +string(flavs[1].Name())+ string(flavs[2].Name())+m_name;
	}
	me_data.histo    = new Histogram(m_type,m_xmin,m_xmax,m_nbins);
	m_all_rates.push_back(me_data);
      }
      m_all_rates[nc].histo->Insert((*bit)->NOutP(),value);
      m_all_rates[nc].sum+=value;
    }
  }
}


void  ME_Rate::EndEvaluation() {
  Primitive_Observable_Base::EndEvaluation();
  for (size_t i=0; i<m_all_rates.size();++i) {
    m_all_rates[i].histo->Finalize();
    m_all_rates[i].histo->Scale(m_all_rates[i].sum/m_sum);
    m_all_rates[i].histo->Output();
  }
}

void  ME_Rate::EndEvaluation(double scale) {
  Primitive_Observable_Base::EndEvaluation(scale);
  for (size_t i=0; i<m_all_rates.size();++i) { 
    m_all_rates[i].histo->Finalize();
    m_all_rates[i].histo->Scale(scale);
    m_all_rates[i].histo->Output();
  }
}


void ME_Rate::Output(const std::string & pname) {
  Primitive_Observable_Base::Output(pname);

  for (size_t i=0;i<m_all_rates.size();++i) {
    m_all_rates[i].histo->Output((pname+std::string("/")+m_all_rates[i].name).c_str());
  }
}

Primitive_Observable_Base & ME_Rate::operator+=(const Primitive_Observable_Base & ob) {
  Primitive_Observable_Base::operator+=(ob);
  if (m_xmin!=ob.Xmin() || m_xmax!=ob.Xmax() || m_nbins!=ob.Nbins()) {
    return *this;
  }

  ME_Rate * mr = ((ME_Rate*)(&ob));
  if (m_all_rates.size()==mr->m_all_rates.size()) {
    for (size_t i=0; i<m_all_rates.size();++i) {
      (*m_all_rates[i].histo) +=(*mr->m_all_rates[i].histo);
    }
  }
  return *this;
}

void ME_Rate::Reset()
{
  Primitive_Observable_Base::Reset();
  for (size_t i=0; i<m_all_rates.size();++i) {
    m_all_rates[i].histo->Reset();
  }
}

//======================================================================


PT_Distribution::PT_Distribution(int _type,double _xmin,double _xmax,int _nbins,
				 int _maxn, Flavour _fl,std::string lname)
{
  m_checkfl=_fl;
  m_type = _type; m_xmin = _xmin; m_xmax = _xmax; m_nbins = _nbins; m_minn=0; m_maxn = _maxn; m_listname=lname;
  m_name  = std::string("pt_dist_")+string(m_checkfl.Name())+string("_");
  if (lname!="FinalState") m_name=lname+std::string("_")+m_name;
  p_histo =  0;
  for (int i=0;i<m_maxn+1;++i)
    m_histos.push_back(new Histogram(m_type,m_xmin,m_xmax,m_nbins));
}

PT_Distribution::PT_Distribution(int _type,double _xmin,double _xmax,int _nbins,
				 int _minn, int _maxn, Flavour _fl,std::string lname)
{
  m_checkfl=_fl;
  m_type = _type; m_xmin = _xmin; m_xmax = _xmax; m_nbins = _nbins; m_minn = _minn; m_maxn = _maxn; m_listname=lname;
  m_name  = std::string("pt_dist_")+string(m_checkfl.Name())+string("_");
  if (lname!="FinalState") m_name=lname+std::string("_")+m_name;
  if (m_minn!=0) {
    MyStrStream str;
    str<<m_name<<m_minn<<"_";
    str>>m_name;
  }

  p_histo =  0;
  for (int i=0;i<m_maxn+1;++i)
    m_histos.push_back(new Histogram(m_type,m_xmin,m_xmax,m_nbins));
}

PT_Distribution::PT_Distribution(PT_Distribution * partner, std::string _prefix)
{
  m_checkfl=partner->m_checkfl;
  m_type = partner->m_type; 
  m_xmin = partner->m_xmin; m_xmax = partner->m_xmax; m_nbins = partner->m_nbins; 
  m_maxn = partner->m_maxn; m_minn = partner->m_minn;
  m_listname = partner->m_listname;
  m_name  = std::string("pt_dist_");
  if (m_listname!="FinalState") m_name=m_listname+std::string("_")+m_name;
  m_name  = _prefix+m_name+string(m_checkfl.Name())+string("_");
  p_histo =  0;
  for (size_t i=0;i<partner->m_histos.size();++i) {
    m_histos.push_back(new Histogram(partner->m_histos[i]));
    m_histos.back()->Reset();
  }
}



void PT_Distribution::Evaluate(const Particle_List & pl,double weight, int ncount) 
{
  PROFILE_HERE;
  if (m_checkfl==Flavour(kf::Z) || m_checkfl==Flavour(kf::W) || m_checkfl==Flavour(kf::W).Bar()) {
    Vec4D mom;
    Vec4D mom_w;
    int count=0;
    for (size_t i=0;i<pl.size();i++) {
      if (pl[i]->Flav().IsLepton()) {
	mom+=pl[i]->Momentum();
	++count;
      }
      if (pl[i]->Flav()==m_checkfl || pl[i]->Flav().Bar()==m_checkfl)
	mom_w=pl[i]->Momentum();
    }
    if (count==2) {
      double pt=sqrt(sqr(mom[1]) + sqr(mom[2]));
      m_histos[0]->Insert(pt,weight,ncount);
    }
    else {
      if (mom_w[0]!=0) {
	double pt=sqrt(sqr(mom_w[1]) + sqr(mom_w[2]));
	m_histos[0]->Insert(pt,weight,ncount);
      }
    }
    return;
  }
    
  vector<double> pts;

  if (pts.size()==0) {
    for (size_t i=0;i<pl.size();i++) {
      if (m_checkfl.Includes(pl[i]->Flav()) || m_checkfl==pl[i]->Flav()) {
	Vec4D mom=pl[i]->Momentum();
	double pt=sqrt(sqr(mom[1]) + sqr(mom[2]));
	pts.push_back(pt);
      }
    }
  }

  std::sort(pts.begin(),pts.end());

  if ((int)pts.size()>=m_minn) {
    for (size_t i=1;i<=pts.size();++i) {
      m_histos[0]->Insert(pts[pts.size()-i],weight,ncount);
      if (i<m_histos.size()) m_histos[i]->Insert(pts[pts.size()-i],weight,ncount);
    }
    for (size_t i=pts.size()+1;i<m_histos.size();++i) m_histos[i]->Insert(0.,0.);
  }
  else {
    for (size_t i=0;i<m_histos.size();++i) m_histos[i]->Insert(0.,0.);
  }
}

void  PT_Distribution::Evaluate(const Blob_List & blobs,double value, int ncount) 
{
  PROFILE_HERE;
  Particle_List * pl=p_ana->GetParticleList(m_listname);
  Evaluate(*pl,value, ncount);
}


void  PT_Distribution::EndEvaluation() 
{
  PROFILE_HERE;
  msg.Debugging()<<"  "<<m_name<<" : "<<std::endl;
  for (size_t i=0; i<m_histos.size();++i) {
    m_histos[i]->Finalize();
    m_histos[i]->Output();
  }
}

void  PT_Distribution::EndEvaluation(double scale) 
{
  PROFILE_HERE;
  msg.Debugging()<<"  "<<m_name<<" : "<<std::endl;
  for (size_t i=0; i<m_histos.size();++i) {
    m_histos[i]->Finalize();
    m_histos[i]->Scale(scale);
    m_histos[i]->Output();
  }
}

void PT_Distribution::Output(const std::string & pname) 
{
  PROFILE_HERE;
  int  mode_dir = 448;
  mkdir((pname).c_str(),mode_dir); 
  for (size_t i=0; i<m_histos.size();++i) {
    std::string fname;
    MyStrStream s;
    s<<i;
    s<<".dat"; 
    s>>fname;
    if (m_name==std::string("pt_dist_")) {
      m_histos[i]->Output((pname+std::string("/pt_dist_")+fname).c_str());
    } 
    else {
      m_histos[i]->Output((pname+std::string("/")+m_name+fname).c_str());
    }
  }
}


Primitive_Observable_Base & PT_Distribution::operator+=(const Primitive_Observable_Base & ob) 
{
  PROFILE_HERE;
  if (m_xmin!=ob.Xmin() || m_xmax!=ob.Xmax() || m_nbins!=ob.Nbins()) {
    std::cout<<" ERROR: in PT_Distribution::operator+="<<std::endl;
    return *this;
  }
  PT_Distribution * ptd = ((PT_Distribution*)(&ob));

  if (m_histos.size()==ptd->m_histos.size()) {
    for (size_t i=0; i<m_histos.size();++i) {
      (*m_histos[i])+=(*ptd->m_histos[i]);
    }
  }
  return *this;
}

void PT_Distribution::Reset() {
  for (size_t i=0; i<m_histos.size();++i) {
    m_histos[i]->Reset();
  }  
}

Primitive_Observable_Base * PT_Distribution::Copy() const {
  return new PT_Distribution(m_type,m_xmin,m_xmax,m_nbins,m_minn,m_maxn,m_checkfl,m_listname);
}




#include "Shower_Observables.H"
#include "MyStrStream.H"
#include "MathTools.H"
#include "Run_Parameter.H"
#include "Vector.H"

#include <algorithm>

using namespace ATOOLS;
using namespace std;


bool ParticleIsInList(const Particle * const p,  const Particle_List & pl) 
{
  bool hit=0;
  for (int j=0;j<pl.size();++j) {
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
  all_obs.multi    = new Multiplicity(00,-0.5,50.5,51,0);
  all_obs.wz_pt    = new PT_Distribution(00,0.,200.,100,1,Flavour(kf::W));
  all_obs.jet_pt   = new PT_Distribution(00,0.,200.,100,6,Flavour(kf::jet));
  all_obs.sum      =0.;
}


Shower_Observables::Shower_Observables(int _type,double _xmin,double _xmax,int _nbins,
      Selector_Base * _sel,int _njet,int _nflavs,int _dohad) 
{
  histo=0;
  type = _type; xmin = _xmin; xmax = _xmax; nbins = _nbins; sel = _sel;
  name  = std::string("shower_obs");
  InitObservables();
}

void Shower_Observables::Evaluate(const Blob_List & blobs ,double value) {
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

  for (Blob_Const_Iterator blit=blobs.begin();blit!=blobs.end();++blit) {
    for (int i=0;i<(*blit)->NOutP();++i) {
      Particle * p = (*blit)->OutParticle(i);

      if (!ParticleIsInList(p,pl) && (*blit)->Status()==1 && 
	  (((p->Info()=='F') || (p->Info()==' ') || p->Info()=='H') && p->Status()!=2 )) {
	pl.push_back(p);
      }
    }
  }

  all_obs.jetrates->Evaluate(pl,value);
  all_obs.multi->Evaluate(pl,value);
  all_obs.wz_pt->Evaluate(pl,value);
  all_obs.jet_pt->Evaluate(pl,value);
  all_obs.sum+=value;

  int nc=0;
  if (do_fl) {
  nc=fl_obs.size();
  for (int i=0;i<nc;++i) {
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
    obs.sum      = 0.;
    fl_obs.push_back(obs);
  }
  fl_obs[nc].jetrates->Evaluate(pl,value);
  fl_obs[nc].multi->Evaluate(pl,value);
  fl_obs[nc].wz_pt->Evaluate(pl,value);
  fl_obs[nc].jet_pt->Evaluate(pl,value);
  fl_obs[nc].sum+=value;
  }


  if (do_jet) {
  // fill histograms (for each number of seed)
  nc=jet_obs.size();
  for (int i=0;i<nc;++i) {
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
    obs.sum      = 0.;
    jet_obs.push_back(obs);
  }
  jet_obs[nc].jetrates->Evaluate(pl,value);
  jet_obs[nc].multi->Evaluate(pl,value);
  jet_obs[nc].wz_pt->Evaluate(pl,value);
  jet_obs[nc].jet_pt->Evaluate(pl,value);
  jet_obs[nc].sum+=value;
  }

  if (do_fl & do_jet) {
  // fill histograms (for each number of seeds AND each initial flavour)
  nc=fl_jet_obs.size();
  for (int i=0;i<nc;++i) {
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
    obs.sum      = 0.;
    fl_jet_obs.push_back(obs);
  }
  fl_jet_obs[nc].jetrates->Evaluate(pl,value);
  fl_jet_obs[nc].multi->Evaluate(pl,value);
  fl_jet_obs[nc].wz_pt->Evaluate(pl,value);
  fl_jet_obs[nc].jet_pt->Evaluate(pl,value);
  fl_jet_obs[nc].sum+=value;
  }
}


void Shower_Observables::EndEvaluation() {
  all_obs.jetrates->EndEvaluation(); 
  all_obs.multi->EndEvaluation();
  all_obs.wz_pt->EndEvaluation(); 
  all_obs.jet_pt->EndEvaluation(); 
  for (int i=0;i<fl_obs.size();++i) {
    double scale =fl_obs[i].sum/all_obs.sum;
    fl_obs[i].jetrates->EndEvaluation(scale); 
    fl_obs[i].multi->EndEvaluation(scale);
    fl_obs[i].wz_pt->EndEvaluation(scale); 
    fl_obs[i].jet_pt->EndEvaluation(scale); 
  }
  for (int i=0;i<jet_obs.size();++i) {
    double scale =jet_obs[i].sum/all_obs.sum;
    jet_obs[i].jetrates->EndEvaluation(scale); 
    jet_obs[i].multi->EndEvaluation(scale);
    jet_obs[i].wz_pt->EndEvaluation(scale); 
    jet_obs[i].jet_pt->EndEvaluation(scale); 
  }
  for (int i=0;i<fl_jet_obs.size();++i) {
    double scale =fl_jet_obs[i].sum/all_obs.sum;
    fl_jet_obs[i].jetrates->EndEvaluation(scale); 
    fl_jet_obs[i].multi->EndEvaluation(scale);
    fl_jet_obs[i].wz_pt->EndEvaluation(scale); 
    fl_jet_obs[i].jet_pt->EndEvaluation(scale); 
  }
}

void Shower_Observables::Output(std::string pname) {
  all_obs.jetrates->Output(pname); 
  all_obs.multi->Output(pname);
  all_obs.wz_pt->Output(pname); 
  all_obs.jet_pt->Output(pname); 
  for (int i=0;i<fl_obs.size();++i) {
    fl_obs[i].jetrates->Output(pname); 
    fl_obs[i].multi->Output(pname);
    fl_obs[i].wz_pt->Output(pname); 
    fl_obs[i].jet_pt->Output(pname); 
  }
  for (int i=0;i<jet_obs.size();++i) {
    jet_obs[i].jetrates->Output(pname); 
    jet_obs[i].multi->Output(pname);
    jet_obs[i].wz_pt->Output(pname); 
    jet_obs[i].jet_pt->Output(pname); 
  }
  for (int i=0;i<fl_jet_obs.size();++i) {
    fl_jet_obs[i].jetrates->Output(pname); 
    fl_jet_obs[i].multi->Output(pname);
    fl_jet_obs[i].wz_pt->Output(pname); 
    fl_jet_obs[i].jet_pt->Output(pname); 
  }
}


//----------------------------------------------------------------------

Jetrates::Jetrates(int _type,double _xmin,double _xmax,int _nbins,
		   Selector_Base * _sel)
{
  partner = 0;
  type = _type; xmin = _xmin; xmax = _xmax; nbins = _nbins; sel = _sel;
  name  = std::string("jetrate");
  histo = 0;
  jfind = new Jet_Finder(rpa.gen.Ycut(),1);

  int had=0;
  if (rpa.gen.Beam1()==Flavour(kf::p_plus) || rpa.gen.Beam1()==Flavour(kf::p_plus).Bar()) {
    had=1;
    xmin=0.01*xmin;
    xmax=0.01*xmax;
  }

  jets.push_back(6);
  histos.push_back(new Histogram(11,xmin,xmax,nbins));
  rates.push_back(new Histogram(11,xmin,xmax,nbins));
  jets.push_back(5);
  histos.push_back(new Histogram(11,xmin,xmax,nbins));
  rates.push_back(new Histogram(11,xmin,xmax,nbins));
  jets.push_back(4);
  histos.push_back(new Histogram(11,xmin,xmax,nbins));
  rates.push_back(new Histogram(11,xmin,xmax,nbins));
  jets.push_back(3);
  histos.push_back(new Histogram(11,xmin,xmax,nbins));
  rates.push_back(new Histogram(11,xmin,xmax,nbins));
  jets.push_back(2);
  histos.push_back(new Histogram(11,xmin,xmax,nbins));
  rates.push_back(new Histogram(11,xmin,xmax,nbins));
  if (had) {
    jets.push_back(1);
    histos.push_back(new Histogram(11,xmin,xmax,nbins));
    rates.push_back(new Histogram(11,xmin,xmax,nbins));
    jets.push_back(0);
    histos.push_back(new Histogram(11,xmin,xmax,nbins));
    rates.push_back(new Histogram(11,xmin,xmax,nbins));
  }
  ymax=-1;
  ymin=2;
};


Jetrates::Jetrates(Jetrates * _partner, std::string _prefix)
{
  partner = _partner;
  type = partner->type; xmin = partner->xmin; xmax = partner->xmax; nbins = partner->nbins; sel = partner->sel;
  name  = _prefix; 
  histo = 0;
  jfind = 0;

  for (int i=0;i<partner->jets.size();++i) {
    jets.push_back(partner->jets[i]);
    histos.push_back(new Histogram(partner->histos[i]));
    histos.back()->Reset();  
    rates.push_back(new Histogram(partner->rates[i]));
    rates.back()->Reset();
  }

  ymax=-1;
  ymin=2;
};



void Jetrates::Evaluate(const Particle_List & pl,double value) {
  ys.clear();
  if (partner==0)  jfind->ConstructJets(&pl,jets,ys);
  else ys=partner->ys;


  for (int k=0; k<jets.size(); ++k ){
    histos[k]->Insert(ys[k], value);
  }

  for (int k=0; k<jets.size()-1; ++k ){
    rates[k]->InsertRange(ys[k],ys[k+1], value);
  }
  rates.back()->InsertRange(ys.back(),1., value);

  ymax=Max(ymax,ys.back());
  ymin=ys.back();
}

void  Jetrates::EndEvaluation() {
  msg.Debugging()<<"  "<<name<<" : "<<std::endl;
  for (int i=0; i<histos.size();++i) {
    histos[i]->Finalize();
    histos[i]->Output();
    rates[i]->Finalize();
    rates[i]->Output();
  }
  msg.Debugging()<<std::endl<<std::endl;
}

void  Jetrates::EndEvaluation(double scale) {
  msg.Debugging()<<"  "<<name<<" : "<<std::endl;
  for (int i=0; i<histos.size();++i) {
    histos[i]->Finalize();
    histos[i]->Scale(scale);
    histos[i]->Output();
    rates[i]->Finalize();
    rates[i]->Scale(scale);
    rates[i]->Output();
  }
  msg.Debugging()<<std::endl<<std::endl;
}

void Jetrates::Output(std::string pname) {
  int  mode_dir = 448;
  mkdir((pname).c_str(),mode_dir); 
  for (int i=0; i<histos.size();++i) {
    std::string fname;
    MyStrStream s;
    s<<jets[i];
    s<<".dat";
    s>>fname;
    if (name==std::string("jetrate")) {
      histos[i]->Output((pname+std::string("/diffrate")+fname).c_str());
      rates[i]->Output((pname+std::string("/totrate")+fname).c_str());
    } 
    else{
      histos[i]->Output((pname+std::string("/")+name+std::string("diffrate")+fname).c_str());
      rates[i]->Output((pname+std::string("/")+name+std::string("totrate")+fname).c_str());
    }
  }
}

//----------------------------------------------------------------------
Multiplicity::Multiplicity(int _type,double _xmin,double _xmax,int _nbins, int mode)
{
  type = _type; xmin = _xmin; xmax = _xmax; nbins = _nbins; sel = 0; m_mode=mode;
  switch (m_mode) {
  case 0:  
    name  = std::string("multi.dat");
    break;
  case 2:
    name  = std::string("multi_charged_hadrons.dat");
    break;
  default:
    name  = std::string("multi_total.dat");
  }
  histo = new Histogram(type,xmin,xmax,nbins);
};

Multiplicity::Multiplicity(Multiplicity * old) {
  type  = old->Type();
  sel   = old->Sel();
  xmin  = old->Xmin();
  xmax  = old->Xmax();
  nbins = old->Nbins();
  name  = old->Name();
  histo = new Histogram(old->Histo());
  histo -> Reset();
}


Multiplicity::Multiplicity(Multiplicity * _partner, std::string _prefix)
{
  partner = _partner;
  type = partner->type; xmin = partner->xmin; xmax = partner->xmax; nbins = partner->nbins; sel = partner->sel;
  name  = _prefix+std::string("multi.dat");
  histo = new Histogram(type,xmin,xmax,nbins);
}


void Multiplicity::Evaluate(const Particle_List & pl,double value) {
  if (m_mode==0) histo->Insert(pl.size()-2,value);
  else histo->Insert(pl.size()-2,value);
}

void  Multiplicity::Evaluate(Vec4D *,Flavour *,double value) {
}

void  Multiplicity::Evaluate(const Blob_List & blobs,double value) {
  Particle_List pl;

  if (m_mode==0) { // compatibity mode
    for (Blob_Const_Iterator blit=blobs.begin();blit!=blobs.end();++blit) {
      if ((*blit)->Type()[0]=='S') {
	for (int i=0;i<(*blit)->NInP();++i) {
	  Particle * p = (*blit)->InParticle(i);
	  if (!ParticleIsInList(p,pl)) pl.push_back(p);
	}
      }
    }
    
    for (Blob_Const_Iterator blit=blobs.begin();blit!=blobs.end();++blit) {
      for (int i=0;i<(*blit)->NOutP();++i) {
	Particle * p = (*blit)->OutParticle(i);
	if (!ParticleIsInList(p,pl) && (*blit)->Status()==1 && 
	    ((p->Info()=='F') || ((*blit)->Type()[0]!='S' &&  p->Info()=='H')))
	  pl.push_back(p);
      }
    }
  }
  else {
    for (Blob_Const_Iterator blit=blobs.begin();blit!=blobs.end();++blit) {
      for (int i=0;i<(*blit)->NOutP();++i) {
	Particle * p = (*blit)->OutParticle(i);
	if (!ParticleIsInList(p,pl) && (*blit)->Status()==1 && p->DecayBlob()==0) {
	  if (m_mode==2 && p->Flav().IsHadron() && p->Flav().Charge()!=0.) {
	    pl.push_back(p);
	  }
	  else if (m_mode!=2) {
	    pl.push_back(p);
	  }
	}
      }
    }
  }


  Evaluate(pl,value);

}
//----------------------------------------------------------------------




ME_Rate::ME_Rate(int _type,double _xmin,double _xmax,int _nbins,std::string prefix)
{
  type = _type; xmin = _xmin; xmax = _xmax; nbins = _nbins; sel = 0;
  name  = prefix+std::string("_rates.dat");
  histo = new Histogram(type,xmin,xmax,nbins);
  sum =0.;
}

ME_Rate::ME_Rate(ME_Rate * old) {
  type  = old->Type();
  sel   = old->Sel();
  xmin  = old->Xmin();
  xmax  = old->Xmax();
  nbins = old->Nbins();
  name  = old->Name();
  histo = new Histogram(old->Histo());
  histo -> Reset();
}


void ME_Rate::Evaluate(const Blob_List & blobs,double value) {
  for (Blob_Const_Iterator bit=blobs.begin();bit!=blobs.end();++bit) {
    if ((*bit)->Type().find("Signal") !=-1 ) {
      histo->Insert((*bit)->NOutP(),value);
      sum+=value;

      int jets=(*bit)->NOutP();
      Flavour flavs[3];
      flavs[0]=(*bit)->OutParticle(0)->Flav();
      if (jets>=3)       flavs[1]=(*bit)->OutParticle(2)->Flav();
      if (jets>=5)       flavs[2]=(*bit)->OutParticle(4)->Flav();
      int nc=all_rates.size();
      for (int i=0;i<nc;++i) {
	if (all_rates[i].flavs[0]==flavs[0]
	    && all_rates[i].flavs[1]==flavs[1] && all_rates[i].flavs[2]==flavs[2] ) {
	  nc=i;
	}
      }

      if (nc==all_rates.size()) {
	ME_Data  me_data;
	me_data.jets     = jets;
	me_data.sum      = 0.;
	me_data.flavs[0] = flavs[0];
	me_data.flavs[1] = flavs[1];
	me_data.flavs[2] = flavs[2];
	if (jets<=2) {
	  me_data.name     = flavs[0].Name()+name;
	}
	else if (jets<=4) {
	  me_data.name     = string(flavs[0].Name())
	                       +string(flavs[1].Name())+name;
	}
	else {
	  me_data.name     = string(flavs[0].Name())
	    +string(flavs[1].Name())+ string(flavs[2].Name())+name;
	}
	me_data.histo    = new Histogram(type,xmin,xmax,nbins);
	all_rates.push_back(me_data);
      }
      all_rates[nc].histo->Insert((*bit)->NOutP(),value);
      all_rates[nc].sum+=value;
    }
  }
}


void  ME_Rate::EndEvaluation() {
  Primitive_Observable_Base::EndEvaluation();
  for (int i=0; i<all_rates.size();++i) {
    all_rates[i].histo->Finalize();
    all_rates[i].histo->Scale(all_rates[i].sum/sum);
    all_rates[i].histo->Output();
  }
}

void  ME_Rate::EndEvaluation(double scale) {
  Primitive_Observable_Base::EndEvaluation(scale);
  for (int i=0; i<all_rates.size();++i) { 
    all_rates[i].histo->Finalize();
    all_rates[i].histo->Scale(scale);
    all_rates[i].histo->Output();
  }
}


void ME_Rate::Output(std::string pname) {
  Primitive_Observable_Base::Output(pname);

  for (int i=0;i<all_rates.size();++i) {
    all_rates[i].histo->Output((pname+std::string("/")+all_rates[i].name).c_str());
  }
}


PT_Distribution::PT_Distribution(int _type,double _xmin,double _xmax,int _nbins,
				 int _maxn, Flavour _fl)
{
  checkfl=_fl;
  type = _type; xmin = _xmin; xmax = _xmax; nbins = _nbins; maxn = _maxn;
  name  = std::string("pt_dist_")+string(checkfl.Name())+string("_");
  histo =  0;
  for (int i=0;i<maxn+1;++i)
    histos.push_back(new Histogram(type,xmin,xmax,nbins));
}

PT_Distribution::PT_Distribution(PT_Distribution * partner, std::string _prefix)
{
  checkfl=partner->checkfl;
  type = partner->type; 
  xmin = partner->xmin; xmax = partner->xmax; nbins = partner->nbins; 
  maxn = partner->maxn;
  name  = _prefix+std::string("pt_dist_")+string(checkfl.Name())+string("_");
  histo =  0;
  for (int i=0;i<partner->histos.size();++i) {
    histos.push_back(new Histogram(partner->histos[i]));
    histos.back()->Reset();
  }
}



void PT_Distribution::Evaluate(const Particle_List & pl,double weight) {

  if (checkfl==Flavour(kf::Z) || checkfl==Flavour(kf::W) || checkfl==Flavour(kf::W).Bar()) {
    Vec4D mom;
    Vec4D mom_w;
    int count=0;
    for (int i=2;i<pl.size();i++) {
      if (pl[i]->Flav().IsLepton()) {
	mom+=pl[i]->Momentum();
	++count;
      }
      if (pl[i]->Flav()==checkfl || pl[i]->Flav().Bar()==checkfl)
	mom_w=pl[i]->Momentum();
    }
    if (count==2) {
      double pt=sqrt(sqr(mom[1]) + sqr(mom[2]));
      histos[0]->Insert(pt,weight);
    }
    else {
      if (mom_w[0]!=0) {
	double pt=sqrt(sqr(mom_w[1]) + sqr(mom_w[2]));
	histos[0]->Insert(pt,weight);
      }
    }
    return;
  }
    
  vector<double> pts;

  for (int i=2;i<pl.size();i++) {
    if (checkfl.Includes(pl[i]->Flav())) {
      Vec4D mom=pl[i]->Momentum();
      double pt=sqrt(sqr(mom[1]) + sqr(mom[2]));
      pts.push_back(pt);
    }
  }

  std::sort(pts.begin(),pts.end());

  for (int i=1;i<=pts.size();++i) {
    histos[0]->Insert(pts[pts.size()-i],weight);
    if (i<histos.size()) histos[i]->Insert(pts[pts.size()-i],weight);
  }
}

void  PT_Distribution::Evaluate(const Blob_List & blobs,double value) {
  Particle_List pl;
  for (Blob_Const_Iterator blit=blobs.begin();blit!=blobs.end();++blit) {
    if ((*blit)->Type()[0]=='S') {
      for (int i=0;i<(*blit)->NInP();++i) {
	Particle * p = (*blit)->InParticle(i);
	if (!ParticleIsInList(p,pl)) pl.push_back(p);
      }
    }
  }

  // looking for Final State Shower blob
  for (Blob_Const_Iterator blit=blobs.begin();blit!=blobs.end();++blit) {
    for (int i=0;i<(*blit)->NOutP();++i) {
      Particle * p = (*blit)->OutParticle(i);
      if ((p->Info()=='F') ||  
	  (p->Info()=='H' && p->Status()!=2 ) && (*blit)->Status()==1 ){
	if (!ParticleIsInList(p,pl)) pl.push_back(p);
      }
    }
  }
  Evaluate(pl,value);

}


void  PT_Distribution::EndEvaluation() {
  msg.Debugging()<<"  "<<name<<" : "<<std::endl;
  for (int i=0; i<histos.size();++i) {
    histos[i]->Finalize();
    histos[i]->Output();
  }
}

void  PT_Distribution::EndEvaluation(double scale) {
  msg.Debugging()<<"  "<<name<<" : "<<std::endl;
  for (int i=0; i<histos.size();++i) {
    histos[i]->Finalize();
    histos[i]->Scale(scale);
    histos[i]->Output();
  }
}

void PT_Distribution::Output(std::string pname) {
  int  mode_dir = 448;
  mkdir((pname).c_str(),mode_dir); 
  for (int i=0; i<histos.size();++i) {
    std::string fname;
    MyStrStream s;
    s<<i;
    s<<".dat"; 
    s>>fname;
    if (name==std::string("pt_dist_")) {
      histos[i]->Output((pname+std::string("/pt_dist_")+fname).c_str());
    } 
    else{
      histos[i]->Output((pname+std::string("/")+name+fname).c_str());
    }
  }
}

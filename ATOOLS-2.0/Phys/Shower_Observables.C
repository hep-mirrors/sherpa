#include "Shower_Observables.H"
#include "MyStrStream.H"
#include "MathTools.H"
#include "Run_Parameter.H"

using namespace APHYTOOLS;
using namespace AORGTOOLS;
using namespace std;

void Shower_Observables::InitObservables() {
  all_obs.flav     = Flavour(kf::none);
  all_obs.jet_ini  = 0;
  all_obs.jetrates = new Jetrates(11,1.e-6,1.,180,0);
  all_obs.multi    = new Multiplicity(00,-0.5,50.5,51,0);
  all_obs.sum      =0.;
}


Shower_Observables::Shower_Observables(int _type,double _xmin,double _xmax,int _nbins,
      APHYTOOLS::Selector_Base * _sel,int _njet,int _nflavs,int _dohad) 
{
  histo=0;
  type = _type; xmin = _xmin; xmax = _xmax; nbins = _nbins; sel = _sel;
  name  = std::string("shower_obs");
  InitObservables();
}

void Shower_Observables::Evaluate(const APHYTOOLS::Blob_List & blobs ,double value) {
  Parton_List pl;

  // determin "leading flavour and njet_ini"
  int njet_ini=0;
  Flavour lfl;
  // looking for a hard blob
  for (Blob_Const_Iterator blit=blobs.begin();blit!=blobs.end();++blit) {
    if ((*blit)->Type()[0]=='H') {
      njet_ini=(*blit)->NOutP();
      lfl =(*blit)->OutParton(0)->Flav();
      for (int i=0;i<(*blit)->NInP();++i) {
	pl.push_back((*blit)->InParton(i));
      }
    }
  }

  // looging for Final State Shower blob
  /*
    // old scheme:
  for (Blob_Const_Iterator blit=blobs.begin();blit!=blobs.end();++blit) {
    if ((*blit)->Type()[0]=='F') {
      for (int i=0;i<(*blit)->NOutP();++i) {
	pl.push_back((*blit)->OutParton(i));
      }
    }
  }
  */
  for (Blob_Const_Iterator blit=blobs.begin();blit!=blobs.end();++blit) {
    for (int i=0;i<(*blit)->NOutP();++i) {
      Parton * p = (*blit)->OutParton(i);
      if (p->Info()=='F')
	pl.push_back(p);
    }
  }


  //  fill histograms (for all events)
  all_obs.jetrates->Evaluate(pl,value);
  all_obs.multi->Evaluate(pl,value);
  all_obs.sum+=value;


  // fill histograms (for each flavour)
  //  cout<<" look for flavour channel :"<<lfl<<endl;
  int nc=fl_obs.size();
  for (int i=0;i<nc;++i) {
    if (fl_obs[i].flav==lfl) 
      nc=i;
  }
  if (nc==fl_obs.size()) {
    cout<<" create flavour channel :"<<lfl<<endl;
    Event_Obi obs;
    obs.flav     = lfl;  
    obs.jet_ini  = 0;
    obs.jetrates = new Jetrates(all_obs.jetrates,lfl.Name()+std::string("_"));
    obs.multi    = new Multiplicity(all_obs.multi,lfl.Name()+std::string("_"));
    obs.sum      = 0.;
    fl_obs.push_back(obs);
  }
  fl_obs[nc].jetrates->Evaluate(pl,value);
  fl_obs[nc].multi->Evaluate(pl,value);
  fl_obs[nc].sum+=value;


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

    cout<<" create "<<njet_ini<<"-jet channel : "<<jname<<endl;

    Event_Obi obs;
    obs.flav     = Flavour(kf::none);  
    obs.jet_ini  = njet_ini;
    obs.jetrates = new Jetrates(all_obs.jetrates,jname);
    obs.multi    = new Multiplicity(all_obs.multi,jname);
    obs.sum      = 0.;
    jet_obs.push_back(obs);
  }
  jet_obs[nc].jetrates->Evaluate(pl,value);
  jet_obs[nc].multi->Evaluate(pl,value);
  jet_obs[nc].sum+=value;

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
    cout<<" create "<<njet_ini<<"-jet, "<<lfl.Name()<<" channel : "<<jname<<endl;

    Event_Obi obs;
    obs.flav     = lfl;  
    obs.jet_ini  = njet_ini;
    obs.jetrates = new Jetrates(all_obs.jetrates,jname);
    obs.multi    = new Multiplicity(all_obs.multi,jname);
    obs.sum      = 0.;
    fl_jet_obs.push_back(obs);
  }
  fl_jet_obs[nc].jetrates->Evaluate(pl,value);
  fl_jet_obs[nc].multi->Evaluate(pl,value);
  fl_jet_obs[nc].sum+=value;



}


void Shower_Observables::EndEvaluation() {
  cout<<" EndEvaluation "<<endl;
  all_obs.jetrates->EndEvaluation(); 
  all_obs.multi->EndEvaluation();
  for (int i=0;i<fl_obs.size();++i) {
    double scale =fl_obs[i].sum/all_obs.sum;
    fl_obs[i].jetrates->EndEvaluation(scale); 
    fl_obs[i].multi->EndEvaluation(scale);
  }
  for (int i=0;i<jet_obs.size();++i) {
    double scale =jet_obs[i].sum/all_obs.sum;
    jet_obs[i].jetrates->EndEvaluation(scale); 
    jet_obs[i].multi->EndEvaluation(scale);
  }
  for (int i=0;i<fl_jet_obs.size();++i) {
    double scale =fl_jet_obs[i].sum/all_obs.sum;
    fl_jet_obs[i].jetrates->EndEvaluation(scale); 
    fl_jet_obs[i].multi->EndEvaluation(scale);
  }
}

void Shower_Observables::Output(std::string pname) {
  cout<<" Output "<<pname<<endl;
  all_obs.jetrates->Output(pname); 
  all_obs.multi->Output(pname);
  for (int i=0;i<fl_obs.size();++i) {
    fl_obs[i].jetrates->Output(pname); 
    fl_obs[i].multi->Output(pname);
  }
  for (int i=0;i<jet_obs.size();++i) {
    jet_obs[i].jetrates->Output(pname); 
    jet_obs[i].multi->Output(pname);
  }
  for (int i=0;i<fl_jet_obs.size();++i) {
    fl_jet_obs[i].jetrates->Output(pname); 
    fl_jet_obs[i].multi->Output(pname);
  }
}


//----------------------------------------------------------------------

Jetrates::Jetrates(int _type,double _xmin,double _xmax,int _nbins,
				     APHYTOOLS::Selector_Base * _sel)
{
  partner = 0;
  type = _type; xmin = _xmin; xmax = _xmax; nbins = _nbins; sel = _sel;
  name  = std::string("jetrate");
  histo = 0;
  jfind = new APHYTOOLS::Jet_Finder(2.,1);

  int had=0;
  if (rpa.gen.Beam1()==Flavour(kf::p_plus) || rpa.gen.Beam1()==Flavour(kf::p_plus).Bar()) {
    had=1;
    xmin=0.01*xmin;
    xmax=0.01*xmax;
  }

  jets.push_back(6);
  histos.push_back(new AMATOOLS::Histogram(11,xmin,xmax,nbins));
  rates.push_back(new AMATOOLS::Histogram(11,xmin,xmax,nbins));
  jets.push_back(5);
  histos.push_back(new AMATOOLS::Histogram(11,xmin,xmax,nbins));
  rates.push_back(new AMATOOLS::Histogram(11,xmin,xmax,nbins));
  jets.push_back(4);
  histos.push_back(new AMATOOLS::Histogram(11,xmin,xmax,nbins));
  rates.push_back(new AMATOOLS::Histogram(11,xmin,xmax,nbins));
  jets.push_back(3);
  histos.push_back(new AMATOOLS::Histogram(11,xmin,xmax,nbins));
  rates.push_back(new AMATOOLS::Histogram(11,xmin,xmax,nbins));
  jets.push_back(2);
  histos.push_back(new AMATOOLS::Histogram(11,xmin,xmax,nbins));
  rates.push_back(new AMATOOLS::Histogram(11,xmin,xmax,nbins));
  if (had) {
    jets.push_back(1);
    histos.push_back(new AMATOOLS::Histogram(11,xmin,xmax,nbins));
    rates.push_back(new AMATOOLS::Histogram(11,xmin,xmax,nbins));
    jets.push_back(0);
    histos.push_back(new AMATOOLS::Histogram(11,xmin,xmax,nbins));
    rates.push_back(new AMATOOLS::Histogram(11,xmin,xmax,nbins));
  }
  ymax=-1;
  ymin=2;
};


Jetrates::Jetrates(Jetrates * _partner, std::string _prefix)
{
  partner = _partner;
  type = partner->type; xmin = partner->xmin; xmax = partner->xmax; nbins = partner->nbins; sel = partner->sel;
  name  = _prefix; //+std::string("jetrate");
  histo = 0;
  jfind = 0;

  for (int i=0;i<partner->jets.size();++i) {
    jets.push_back(partner->jets[i]);
    histos.push_back(new AMATOOLS::Histogram(partner->histos[i]));
    histos.back()->Reset();  // clear histogram 
    rates.push_back(new AMATOOLS::Histogram(partner->rates[i]));
    rates.back()->Reset();
  }

  ymax=-1;
  ymin=2;
};



void Jetrates::Evaluate(const APHYTOOLS::Parton_List & pl,double value) {
  ys.clear();
  if (partner==0)  jfind->ConstructJets(&pl,jets,ys);
  else ys=partner->ys;


  // differential n+1 -> n Jetrates 
  for (int k=0; k<jets.size(); ++k ){
    histos[k]->Insert(ys[k], value);
  }

  // total n Jetrates(ycut)
  for (int k=0; k<jets.size()-1; ++k ){
    rates[k]->InsertRange(ys[k],ys[k+1], value);
  }
  rates.back()->InsertRange(ys.back(),1., value);

  // some extra statistics
  ymax=AMATOOLS::Max(ymax,ys.back());
  ymin=ys.back();
}

void  Jetrates::EndEvaluation() {
  AORGTOOLS::msg.Tracking()<<"  "<<name<<" : "<<std::endl;
  for (int i=0; i<histos.size();++i) {
    histos[i]->Finalize();
    histos[i]->Output();
    rates[i]->Finalize();
    rates[i]->Output();
  }
  AORGTOOLS::msg.Tracking()<<std::endl<<std::endl;

  std::cout<<name<<" ymax ="<<ymax<<std::endl;
}

void  Jetrates::EndEvaluation(double scale) {
  AORGTOOLS::msg.Tracking()<<"  "<<name<<" : "<<std::endl;
  for (int i=0; i<histos.size();++i) {
    histos[i]->Finalize();
    histos[i]->Scale(scale);
    histos[i]->Output();
    rates[i]->Finalize();
    rates[i]->Scale(scale);
    rates[i]->Output();
  }
  AORGTOOLS::msg.Tracking()<<std::endl<<std::endl;

  std::cout<<" =================== "<<std::endl;
  std::cout<<" ymax ="<<ymax<<std::endl;
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

Multiplicity::Multiplicity(Multiplicity * _partner, std::string _prefix)
{
  partner = _partner;
  type = partner->type; xmin = partner->xmin; xmax = partner->xmax; nbins = partner->nbins; sel = partner->sel;
  name  = _prefix+std::string("multi.dat");
  histo = new AMATOOLS::Histogram(type,xmin,xmax,nbins);
}


void Multiplicity::Evaluate(const APHYTOOLS::Parton_List & pl,double value) {
  //  cout<<" multi = "<<pl.size()<<endl;
  histo->Insert(pl.size()-2,value);
}

void  Multiplicity::Evaluate(AMATOOLS::Vec4D *,APHYTOOLS::Flavour *,double value) {
  //  histo->Insert(nout, value);  
}

void  Multiplicity::Evaluate(const APHYTOOLS::Blob_List & blobs,double value) {
  Parton_List pl;
  for (Blob_Const_Iterator blit=blobs.begin();blit!=blobs.end();++blit) {
    if ((*blit)->Type()[0]=='H') {
      for (int i=0;i<(*blit)->NInP();++i) {
	pl.push_back((*blit)->InParton(i));
      }
    }
  }

  // looging for Final State Shower blob
  for (Blob_Const_Iterator blit=blobs.begin();blit!=blobs.end();++blit) {
    for (int i=0;i<(*blit)->NOutP();++i) {
      Parton * p = (*blit)->OutParton(i);
      if (p->Info()=='F')
	pl.push_back(p);
    }
  }
  Evaluate(pl,value);

}







void ME_Rate::Evaluate(const APHYTOOLS::Blob_List & blobs,double value) {
  for (Blob_Const_Iterator bit=blobs.begin();bit!=blobs.end();++bit) {
    if ((*bit)->Type().find("Hard") !=-1 ) {
      // fill jet number
      histo->Insert((*bit)->NOutP(),value);
      sum+=value;

      // fill detailed process history:
      int jets=(*bit)->NOutP();
      Flavour flavs[3];
      flavs[0]=(*bit)->OutParton(0)->Flav();
      if (jets>=3)       flavs[1]=(*bit)->OutParton(2)->Flav();
      if (jets>=5)       flavs[2]=(*bit)->OutParton(4)->Flav();
      int nc=all_rates.size();
      for (int i=0;i<nc;++i) {
	if (all_rates[i].flavs[0]==flavs[0]
	    && all_rates[i].flavs[1]==flavs[1] && all_rates[i].flavs[2]==flavs[2] ) {
	  nc=i;
	}
      }

      if (nc==all_rates.size()) {
	//	cout<<" create rate channel :"<<(*(*bit))<<endl;
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
	me_data.histo    = new AMATOOLS::Histogram(type,xmin,xmax,nbins);
	cout<<" create rate channel :"<<me_data.name<<endl;
	all_rates.push_back(me_data);
      }
      all_rates[nc].histo->Insert((*bit)->NOutP(),value);
      all_rates[nc].sum+=value;
    }
  }
}


void  ME_Rate::EndEvaluation() {
  cout<<" Scale intern "<<name<<endl;
  Primitive_Observable_Base::EndEvaluation();
  for (int i=0; i<all_rates.size();++i) {
    all_rates[i].histo->Finalize();
    all_rates[i].histo->Scale(all_rates[i].sum/sum);
    all_rates[i].histo->Output();
  }
}

void  ME_Rate::EndEvaluation(double scale) {
  cout<<" Scale extern "<<name<<endl;
  Primitive_Observable_Base::EndEvaluation(scale);
  for (int i=0; i<all_rates.size();++i) { 
    all_rates[i].histo->Finalize();
    all_rates[i].histo->Scale(scale);
    all_rates[i].histo->Output();
  }
}


void ME_Rate::Output(std::string pname) {
  // write histo
  cout<<" Output "<<pname<<endl;
  Primitive_Observable_Base::Output(pname);

  for (int i=0;i<all_rates.size();++i) {
    all_rates[i].histo->Output((pname+std::string("/")+all_rates[i].name).c_str());
  }
}

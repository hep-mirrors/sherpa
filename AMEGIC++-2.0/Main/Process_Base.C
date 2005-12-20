#include "Process_Base.H"
#include "Run_Parameter.H"
//#include "Standard_Selector.H"
#include "Combined_Selector.H"
#include "Data_Collector.H"
#include "Message.H"

#include "Running_AlphaS.H"
#include "Running_AlphaQED.H"

#include <algorithm>
#include <stdio.h>

using namespace AMEGIC;
using namespace PHASIC;
using namespace MODEL;
using namespace BEAM;
using namespace PDF;
using namespace ATOOLS;
using namespace std;


/*------------------------------------------------------------------------------
  
  Constructor

  ------------------------------------------------------------------------------*/

Process_Base::Process_Base(): 
  Integrable_Base(0,0),
  m_gen_str(3),p_b(0),p_flin(0),p_flout(0),
  p_pl(0),p_plin(0),p_plout(0), 
  m_maxfac(1.), m_ycut(-1.),
  m_maxerror(-1.), p_psgen(0), m_print_graphs(false), p_pinfo(0)
{
  m_atoms=1;
  m_analyse=m_tables=0;

  m_n=m_kfactorscheme=m_scalescheme=0;
  m_nstrong=m_neweak=m_orderQCD=m_orderEW=0;
  m_totalxs=m_totalerr=m_totalsum=m_totalsumsqr=m_max=0.;
  m_last=m_lastdxs=0.;
  m_lastlumi=1.;
  m_scale[stp::as]=sqr(rpa.gen.Ecms());
  m_scale[stp::fac]=sqr(rpa.gen.Ecms());
  m_threshold=0.;
  m_updatescales=false;
}



Process_Base::Process_Base(Process_Info* pinfo,int _nin,int _nout,ATOOLS::Flavour * _fl,
			   PDF::ISR_Handler * _isr,BEAM::Beam_Spectra_Handler * _beam,
			   int _gen_str, int _orderQCD, int _orderEW,
			   int _scalescheme,int _kfactorscheme,double _scale,
			   Pol_Info * _pl,
			   int _nex,ATOOLS::Flavour * _ex_fl,double ycut,double error) :
  Integrable_Base(_nin,_nout,_scalescheme,_kfactorscheme,_beam,_isr),
  m_gen_str(_gen_str), m_nex(_nex),
  p_ex_fl(_ex_fl),
  m_atoms(0), m_analyse(0), m_tables(0), 
  m_maxfac(1.),
  m_orderQCD(_orderQCD), m_orderEW(_orderEW),
  m_ycut(ycut),m_maxerror(error), 
  p_psgen(0), m_print_graphs(false), p_pinfo(pinfo)
{

  m_scale[stp::as]=m_scale[stp::fac]=_scale;

  if (m_scale[stp::as]<0.) {
    m_scale[stp::as]=m_scale[stp::fac]=sqr(rpa.gen.Ecms());
  }

  m_updatescales=false;
  if (m_scalescheme<0 || m_kfactorscheme<0) m_updatescales = true;
  
  if (m_scalescheme<0)   m_scalescheme   = -m_scalescheme;
  if (m_kfactorscheme<0) m_kfactorscheme = -m_kfactorscheme;
  
  p_flavours = 0;
  p_pl = 0;
  p_pshandler = 0;
  p_selector = 0;
  p_flin  = 0;
  p_flout = 0;
  p_plin  = 0;
  p_plout = 0;

  if (m_nin==0) return;

  p_flin    = new Flavour[m_nin];
  p_flout   = new Flavour[m_nout];  
  p_plin    = new Pol_Info[m_nin];
  p_plout   = new Pol_Info[m_nout]; 
  
  for (size_t i=0;i<m_nin;i++) {
    p_flin[i]  = _fl[i];
    if (_pl!=0) p_plin[i] = _pl[i];
           else p_plin[i] = Pol_Info(p_flin[i]); 
    if (p_flin[i].Strong())        m_nstrong++;
    if (p_flin[i].IntCharge()!=0)  m_neweak++;
  }

//   pinfo->Print();
//   pinfo->FullPrint();
  pinfo->Reshuffle();
  pinfo->GetTotalFlavList(p_flout);
  pinfo->GetTotalPolList(p_plout);
  for (size_t i=0;i<m_nout;i++) { 
    if (p_flout[i].Strong())       m_nstrong++;
    if (p_flout[i].IntCharge()!=0) m_neweak++;
  }
}

Process_Base::~Process_Base() {
  if (p_flin)     { delete [] p_flin;  p_flin     = 0; }
  if (p_flout)    { delete [] p_flout; p_flout    = 0; }
  if (p_pl)       { delete [] p_pl;    p_pl       = 0; }
  if (p_plin)     { delete [] p_plin;  p_plin     = 0; }
  if (p_plout)    { delete [] p_plout; p_plout    = 0; }
  if (p_pshandler)       { delete p_pshandler;       p_pshandler       = 0; }
}


/*------------------------------------------------------------------------------
  
  Naming

  ------------------------------------------------------------------------------*/

string * Process_Base::GenerateNames(int _nin, Flavour * _flin, Pol_Info * _plin,
				     string & _name,string & _ptype, string & _lib, Process_Info* pi)
{
  Process_Info *ppi = pi;
  if (ppi==NULL) ppi=p_pinfo;
  Reshuffle(_nin, _flin, _plin);

  if (_flin[0].IsAnti() && !_flin[1].IsAnti()) {
    Flavour flhelp  = _flin[0];
    _flin[0] = _flin[1];
    _flin[1] = flhelp;
    Pol_Info plhelp  = _plin[0];
    _plin[0] = _plin[1];
    _plin[1] = plhelp;    
  }
  int _nout=ppi->TotalNout();
  
  _name=ToString(_nin)+"_"+ToString(_nout);
  if (m_gen_str>1) _ptype      = string("P")+_name;
  else _ptype      = string("N")+_name;
  _lib        = _name;
  _name      += string("_");
  
  for (size_t i=0;i<m_nin;i++) {
    _name += _flin[i].IDName();
    // polinfo for fully polarised incommings
    if (_plin[i].pol_type=='c' && _plin[i].num==1) {
      if (_plin[i].type[0]==-1) _name += string("m");
      if (_plin[i].type[0]==+1) _name += string("p");
      if (_plin[i].type[0]==0)  _name += string("z");
    }
    else if (_plin[i].pol_type=='h' && _plin[i].num==1) {
      if (_plin[i].type[0]==-1) _name += string("m");
      if (_plin[i].type[0]==+1) _name += string("p");
    }


    _name += string("_");
  }
  _name += string("__");

  _name +=ppi->GenerateName();
  return &_name;
}


class S_Data {
public:
  S_Data(const int _i, Flavour _fl, Pol_Info _pl) :
  i(_i),fl(_fl),pl(_pl) {}
  int      i;
  Flavour  fl;
  Pol_Info pl;
};

class Order_FVST {
public:
  int operator()(const S_Data & a, const S_Data & b) {
    //    if "a < b" return 1  else 0;
    if (a.fl.IsFermion() && !b.fl.IsFermion()) return 1;
    if (a.fl.IsVector() && !b.fl.IsFermion() && !b.fl.IsVector()) return 1;
    if (a.fl.IsScalar() && !b.fl.IsScalar() && 
	 !b.fl.IsFermion() && !b.fl.IsVector()) return 1;
    return 0;
  }
};

class Order_SVFT {
public:
  int operator()(const S_Data & a, const S_Data & b) {
    //    if "a < b" return 1  else 0;
    if (a.fl.IsScalar() && !b.fl.IsScalar()) return 1;
    if (a.fl.IsVector() && !b.fl.IsScalar() && !b.fl.IsVector()) return 1;
    if (a.fl.IsFermion() && !b.fl.IsFermion() && 
	 !b.fl.IsScalar() && !b.fl.IsVector()) return 1;
    return 0;
  }
};

class Order_Mass {
public:
  int operator()(const S_Data & a, const S_Data & b) {
    //    if "a > b" return 1  else 0;
    if (a.fl.Mass() <= b.fl.Mass()) return 0;
    return 1;
  }
};

class Order_InvMass {
public:
  int operator()(const S_Data & a, const S_Data & b) {
    //    if "a < b" return 1  else 0;
    if (a.fl.Mass() < b.fl.Mass()) return 1;
    return 0;
  }
};


class Order_Kfc {
public:
  int operator()(const S_Data & a, const S_Data & b) {
    //    if "a < b" return 1  else 0;
    if (a.fl.Kfcode() < b.fl.Kfcode()) return 1;
    return 0;
  }
};


class Order_Anti {
public:
  int operator()(const S_Data & a, const S_Data & b) {
    //    if "a < b" return 1  else 0;
    if ((a.fl.IsFermion() && b.fl.IsFermion())
	&& (!a.fl.IsAnti() && b.fl.IsAnti())) return 1;
    return 0;
  }
};


class Order_Coupling {
public:
  int operator()(const S_Data & a, const S_Data & b) {
      if (!a.fl.Strong() && b.fl.Strong()) return 1;
      return 0;
  }
};

//
// Note: all order operator have to return 0 if 
//       two elements are equal!
//       Otherwise the order will change even for
//       equal elements.
//

void Process_Base::Reshuffle(int n, Flavour* flav, Pol_Info* plav)
{
  std::vector<S_Data> sd;
  Flavour heaviest(kf::photon);
  for (int i=0;i<n;++i) {
    sd.push_back(S_Data(i,flav[i],plav[i]));
    if (flav[i].Mass()>heaviest.Mass()) heaviest=flav[i];
    else if (flav[i].Mass()==heaviest.Mass() &&
	     !flav[i].IsAnti()) heaviest=flav[i];
  }

  std::stable_sort(sd.begin(),sd.end(),Order_Kfc());
  std::stable_sort(sd.begin(),sd.end(),Order_Anti());
  std::stable_sort(sd.begin(),sd.end(),Order_SVFT());
  if (heaviest.IsAnti())  std::stable_sort(sd.begin(),sd.end(),Order_InvMass());
  else   std::stable_sort(sd.begin(),sd.end(),Order_Mass());
  std::stable_sort(sd.begin(),sd.end(),Order_Coupling());
  
  
  for (int i=0;i<n;++i) {
    size_t pos=m_efunc.find("p["+ToString(sd[i].i)+"]");
    if (pos!=std::string::npos) 
      m_efunc.replace(pos,3+ToString(sd[i].i).length(),"p["+ToString(i)+"]");
    flav[i]=sd[i].fl;
    plav[i]=sd[i].pl;
  }
}


bool Process_Base::IsFile(string filename)
{
  ifstream from;
  bool     hit = 0;
  from.open(filename.c_str());
  if (from) hit = 1;
  from.close();
  return hit;
}

/*------------------------------------------------------------------------------
  
  Process initialization
  
  ------------------------------------------------------------------------------*/


void Process_Base::AddChannels(Process_Base * _proc) 
{
//   if (m_nin!=2) return;
  for (size_t i=0;i<_proc->Size();i++) {
    if ((*_proc)[i]->Partner()==NULL) AddChannels((*_proc)[i]);
    else {
      if ((*_proc)[i]->Partner()==(*_proc)[i]) {
	list<string>* clist = (*_proc)[i]->PSHandler()->GetChannelLibNames();
	list<string>* tlist = PSHandler()->GetChannelLibNames();
	for (list<string>::iterator it=clist->begin();it!=clist->end();++it) {
	  bool hit = 0;
	  for (list<string>::iterator jt=tlist->begin();jt!=tlist->end();++jt) {
	    if ((*it)==(*jt)) {
	      //cout<<"Process_Base::AddChannels: "<<(*it)<<"/"<<(*jt)<<endl;
	      hit = 1;
	      break;
	    }
	  }
	  if (!hit) tlist->push_back((*it));
	}
      }
    }
  }
}
/*void Process_Base::AddChannels(Process_Base * _proc,Multi_Channel * _fsr,
			       vector<Channel_Info> & _beamparams,
			       vector<Channel_Info> & _isrparams) {

  if (m_nin!=2) return;
  bool         addit;
  Channel_Info ci;

  for (size_t i=0;i<_proc->Size();i++) {
    if ((*_proc)[i]->Partner()==NULL) AddChannels((*_proc)[i],_fsr,_beamparams,_isrparams);
    else {
      if ((*_proc)[i]->Partner()==(*_proc)[i]) {
	Single_Channel * sc;
	int next; string chname;
	for (int j=0;j<(*_proc)[i]->NumberOfFSRIntegrators();j++) { 
	  chname = ((*_proc)[i]->FSRIntegrator(j))->Name();
	  if ( (chname!=string("Rambo")) && (chname!=string("RamboKK")) 
	       && (chname!=string("Sarge")) && (chname!=string("Decay2-Channel 1")) ) { 
	    next   = chname.find(string("--"));
	    chname = chname.substr(next+2);
	    sc     = (*_proc)[i]->PSGenerator()->SetChannel(m_nin,m_nout,p_flavours,chname);
	    sc->SetName(((*_proc)[i]->FSRIntegrator(j))->Name());
	    _fsr->Add( sc );
	  }
	}

	if ((*_proc)[i]->Beam()) {
	  if ((*_proc)[i]->Beam()->On()>0) {
	    for (int j=0;j<(*_proc)[i]->NumberOfBeamIntegrators()/3;j++) {
	      (*_proc)[i]->BeamChannels(j,ci);
	      addit = 1;
	      for (size_t k=0;k<_beamparams.size();k++) {
		if (_beamparams[k]==ci) { addit = 0; break; }
	      }
	      if (addit) _beamparams.push_back(ci);
	    }
	  }
	}
	if ((*_proc)[i]->ISR()) {
	  if ((*_proc)[i]->ISR()->On()>0) {
	    for (int j=0;j<(*_proc)[i]->NumberOfISRIntegrators()/3;j++) {
	      (*_proc)[i]->ISRChannels(j,ci);
	      addit = 1;
	      for (size_t k=0;k<_isrparams.size();k++) {
		if (_isrparams[k]==ci) { addit = 0; break; }
	      }
	      if (addit) _isrparams.push_back(ci);
	    }
	  }
	}
      }
    }
  }
  }*/



/*------------------------------------------------------------------------------

  Process management
  
  ------------------------------------------------------------------------------*/

void Process_Base::SetName(string _name)                { m_name    = _name;   }
void Process_Base::SetResDir(string _resdir)            { m_resdir  = _resdir; }
void Process_Base::SetAtoms(bool _atoms)                { m_atoms   = _atoms;  }
void Process_Base::SetTables(bool _tables)              { m_tables  = _tables; }
void Process_Base::SetBeam(Beam_Spectra_Handler * _beam){ p_beamhandler    = _beam;   }
void Process_Base::SetISR(ISR_Handler * _isr)           { p_isrhandler     = _isr;    }
void Process_Base::SetSelector(Selector_Base * _sel)    { p_selector     = _sel;    }
void Process_Base::SetMomenta(ATOOLS::Vec4D * _moms)  { p_momenta    = _moms;   }
void Process_Base::SetNStrong(int _nstrong)             { m_nstrong = _nstrong;}
void Process_Base::SetNEWeak(int _neweak)               { m_neweak  = _neweak; }
void Process_Base::SetMax(const double max, int depth)  
{
  if (max!=0.) m_max     = max;     
} 
void Process_Base::SetMaxJetNumber(int max)             { m_maxjetnumber  = max;    } 
void Process_Base::SetScales(double q2_fac, double q2_ren)
{ 
  //  std::cout<<"Process_Base::SetScales("<<q2_fac<<","<<q2_ren<<") : "<<Name()<<std::endl;

  m_scale[stp::fac] = rpa.gen.FactorizationScaleFactor() * q2_fac;  
  m_scale[stp::as]  = rpa.gen.RenormalizationScaleFactor() * q2_ren;
} 

void Process_Base::AddToDataCollector(int i)
{
  std::string name;
  for (size_t j=0; j<m_nin; ++j) name+=p_flavours[j].TexName()+"\\,";
  name+="\\to\\,";
  for (size_t j=m_nin; j<m_nin+m_nout; ++j) name+=p_flavours[j].TexName()+"\\,";
  ATOOLS::Process_Info pi(name,m_totalxs*rpa.Picobarn(),m_totalerr*rpa.Picobarn());
  Data_Collector::AddData("PROCESS"+ToString(i),new Blob_Data<ATOOLS::Process_Info>(pi));
}

/*------------------------------------------------------------------------------

  Calculating total cross sections and single event generation

  ------------------------------------------------------------------------------*/

void Process_Base::RescaleXSec(double fac) {
  m_rfactor  *= fac;
  m_totalxs  *= fac;
  m_totalsum *= fac;
  m_totalerr *= fac;
  m_max      *= fac;
  m_totalsumsqr *= fac*fac; // only an estimate
}

void Process_Base::SetupEnhance() {
  if (m_enhancefac==1. && m_maxfac==1. && m_maxeps==0.) return;
  if (m_enhancefac!=1.) {
    double xs=TotalXS();
    SetTotalXS(xs*m_enhancefac);
  }
  if (m_maxeps>0.) {
    double max = GetMaxEps(m_maxeps);
    msg.Info()<<"Maximum reduction factor for "<<m_name<<": "<<Max()/max<<" (epsilon="<<m_maxeps<<")"<<endl;
    SetMax(max);
  }
  else if (m_maxfac!=1.) {
    double max=Max();
    SetMax(max*m_maxfac);
  }
}

void Process_Base::SetEnhance(double enhancefac, double maxfac, double epsilon) 
{
  m_enhancefac = enhancefac;
  m_maxfac     = maxfac;
  m_maxeps     = epsilon;
}


/*------------------------------------------------------------------------------
  
  Access methods
  
  ------------------------------------------------------------------------------*/

void Process_Base::SetScale(const double scale)         
{ 
  //  std::cout<<"Process_Base::SetScale("<<scale<<")"<<std::endl;
  m_scale[stp::as]=m_scale[stp::fac]=scale; 
}

string                  Process_Base::ResDir()                       { return m_resdir; }
string                  Process_Base::LibName()                      { return string("error"); }
int                     Process_Base::NumberOfDiagrams()             { return 0; }
Point                 * Process_Base::Diagram(int i)                 { return 0; }
bool                    Process_Base::IsFreeOfFourVertex(Point * _p) { return 1; }
Phase_Space_Generator * Process_Base::PSGenerator()                  { return p_psgen; }
double                  Process_Base::FactorisationScale()           { return m_scale[stp::fac]; }

int                     Process_Base::ISRNumber()                                        { return 0; }
int                     Process_Base::BeamNumber()                                       { return 0; }
void                    Process_Base::ISRInfo(int,int &,double &,double &)              { return; }
void                    Process_Base::BeamInfo(int,int &,double &,double &)             { return; }

void Process_Base::BeamChannels(int i,Channel_Info & ci) { p_pshandler->BeamChannels(i,ci); }
void Process_Base::ISRChannels(int i,Channel_Info & ci)  { p_pshandler->ISRChannels(i,ci); }
int              Process_Base::NumberOfBeamIntegrators() { return p_pshandler->NumberOfBeamIntegrators(); }
int              Process_Base::NumberOfISRIntegrators()  { return p_pshandler->NumberOfISRIntegrators(); }
int              Process_Base::NumberOfFSRIntegrators()  { return p_pshandler->NumberOfFSRIntegrators(); }
Multi_Channel *  Process_Base::FSRIntegrator()           { return p_pshandler->FSRIntegrator(); }
Single_Channel * Process_Base::FSRIntegrator(int i)      { return p_pshandler->FSRIntegrator(i); }

void Process_Base::SwapInOrder() {
  if (m_nin!=2) return;
  Flavour help = p_flavours[0];
  p_flavours[0] = p_flavours[1];
  p_flavours[1] = help;
  Vec4D mom = p_momenta[0];
  p_momenta[0] = p_momenta[1];
  p_momenta[1] = mom;
  m_swaped = 1;
}

void Process_Base::RestoreInOrder() {
  if (m_nin==2 && m_swaped) {
    Flavour help = p_flavours[0];
    p_flavours[0] = p_flavours[1];
    p_flavours[1] = help;
    m_swaped = 0;
  }
}

void Process_Base::SetPrintGraphs(bool print_graphs) 
{
 m_print_graphs=print_graphs; 
}

void Process_Base::SetWEventMode(int mode) {}

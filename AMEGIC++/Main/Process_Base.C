#include "Process_Base.H"
#include "Run_Parameter.H"
#include "Combined_Selector.H"
#include "Message.H"

#include "Running_AlphaS.H"
#include "Running_AlphaQED.H"

#include <algorithm>
#include <sys/stat.h>
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
  m_maxfac(1.), 
  m_maxerror(-1.), p_psgen(0), m_print_graphs(false), p_pinfo(0)
{
  m_atoms=1;
  m_analyse=m_tables=0;

  m_n=m_kfactorscheme=0;
  m_scalescheme=PHASIC::scl::unknown;
  m_nstrong=m_neweak=m_orderQCD=m_orderEW=0;
  m_totalxs=m_totalerr=m_totalsum=m_totalsumsqr=m_max=0.;
  m_last=m_lastdxs=0.;
  m_lastlumi=1.;
  m_scale[stp::ren]=sqr(rpa.gen.Ecms());
  m_scale[stp::fac]=sqr(rpa.gen.Ecms());
  m_threshold=0.;
}



Process_Base::Process_Base(Process_Info* pinfo,int _nin,int _nout,ATOOLS::Flavour * _fl,
			   PDF::ISR_Handler * _isr,BEAM::Beam_Spectra_Handler * _beam,
			   int _gen_str, int _orderQCD, int _orderEW,
			   PHASIC::scl::scheme _scalescheme,int _kfactorscheme,double _scale,
			   Pol_Info * _pl,
			   int _nex,ATOOLS::Flavour * _ex_fl,std::string cuttag,double error) :
  Integrable_Base(_nin,_nout,_scalescheme,_kfactorscheme,_beam,_isr),
  m_gen_str(_gen_str), m_nex(_nex), p_b(NULL),
  p_ex_fl(_ex_fl),
  m_atoms(0), m_analyse(0), m_tables(0), 
  m_maxfac(1.),
  m_maxerror(error), 
  p_psgen(0), m_print_graphs(false), p_pinfo(pinfo)
{
  if (pinfo!=NULL) m_corenout=pinfo->Nout();
  m_orderQCD=_orderQCD;
  m_orderEW=_orderEW;
  if (cuttag.length()>0) m_cuttag=cuttag;
  m_scale[stp::ren]=m_scale[stp::fac]=_scale;

  if (m_scale[stp::ren]<0.) {
    m_scale[stp::ren]=m_scale[stp::fac]=sqr(rpa.gen.Ecms());
  }

  p_flavours = 0;
  p_b  = 0;
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

  GenerateName(m_nin,p_flin,p_plin,m_name, pinfo);
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
  if (p_pinfo) delete p_pinfo;
  if (p_b) delete [] p_b;
}


/*------------------------------------------------------------------------------
  
  Naming

  ------------------------------------------------------------------------------*/

string * Process_Base::GenerateName(int _nin, Flavour * _flin, Pol_Info * _plin,
				     string & _name, Process_Info* pi)
{
  Process_Info *ppi = pi;
  if (ppi==NULL) ppi=p_pinfo;

  Process_Info ini(0,0);
  ini.AddSubList(_nin,_flin,_plin);
  ini.Reshuffle(ppi);
  ppi->Reshuffle(&ini);

  if (_nin==2) {
    if (ini.m_sublist[0][0]->Flav()->IsAnti() && !ini.m_sublist[0][1]->Flav()->IsAnti()) {
      _flin[0] = *(ini.m_sublist[0][1]->Flav());
      _flin[1] = *(ini.m_sublist[0][0]->Flav());
      _plin[0] = *(ini.m_sublist[0][1]->Pol());
      _plin[1] = *(ini.m_sublist[0][0]->Pol());
    }
    else {
      _flin[0] = *(ini.m_sublist[0][0]->Flav());
      _flin[1] = *(ini.m_sublist[0][1]->Flav());
      _plin[0] = *(ini.m_sublist[0][0]->Pol());
      _plin[1] = *(ini.m_sublist[0][1]->Pol());
    }
  }
  int _nout=ppi->TotalNout();
  _name = ToString(_nin)+"_"+ToString(_nout)+"__";

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


    _name += string("__");
  }

  _name +=ppi->GenerateName();
  return &_name;
}




bool Process_Base::IsFile(string filename)
{
  struct stat fst;
  if (stat(filename.c_str(),&fst)!=-1) {
    return (fst.st_mode&S_IFMT)==S_IFREG;
  }
  return false;
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
void Process_Base::SetCoreMaxJetNumber(int max)         { m_coremaxjetnumber = max; } 

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
    if (m_enhancefac>0.0) SetTotalXS(xs*m_enhancefac);
    else SetTotalXS(-m_enhancefac);
  }
  if (m_maxeps>0.) {
    double max = GetMaxEps(m_maxeps);
    msg_Info()<<"Maximum reduction factor for "<<m_name<<": "<<Max()/max<<" (epsilon="<<m_maxeps<<")"<<endl;
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

string                  Process_Base::ResDir()                       { return m_resdir; }
string                  Process_Base::LibName()                      { return string("error"); }
int                     Process_Base::NumberOfDiagrams()             { return 0; }
Point                 * Process_Base::Diagram(int i)                 { return 0; }
bool                    Process_Base::IsFreeOfFourVertex(Point * _p) { return 1; }
Phase_Space_Generator * Process_Base::PSGenerator()                  { return p_psgen; }

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

void Process_Base::TestPoint(Vec4D *tp)
{
  Flavour* flavs=new Flavour[m_nin+m_nout];
  Vec4D *hmom=new Vec4D[m_nin+m_nout];
  vector<Process_Info*> decaylist;
  for (size_t i=0;i<m_nin;i++) flavs[i]=p_flavours[i];
  size_t n=p_pinfo->GetOnshellFlavList(&flavs[m_nin],decaylist);
  p_pshandler->TestPoint(hmom,m_nin,n,flavs);
  for (size_t i=0;i<m_nin;i++) tp[i]=hmom[i];
  size_t cnt=m_nin;
  for (size_t i=0;i<n;i++) {
    if (decaylist[i]==0) tp[cnt++]=hmom[i+m_nin];
    else {
      tp[cnt]=hmom[i+m_nin];
      DecayPoint(tp,decaylist[i],cnt);
    }
  }
  delete[] flavs;
  delete[] hmom;
}

void Process_Base::DecayPoint(Vec4D *tp,Process_Info* pinfo,size_t &cnt)
{
  size_t nout=pinfo->TotalNout();
  Flavour* flavs=new Flavour[1+nout];
  Vec4D *hmom=new Vec4D[1+nout];
  flavs[0]=*(pinfo->p_fl);
  hmom[0]=Vec4D(tp[cnt].Abs(),0.,0.,0.);
  vector<Process_Info*> decaylist;
  size_t n=pinfo->GetOnshellFlavList(&flavs[1],decaylist);
  p_pshandler->TestPoint(hmom,1,n,flavs);
  Poincare bst(tp[cnt]);
  for (size_t i=0;i<n;i++) {
    bst.BoostBack(hmom[i+1]);
    if (decaylist[i]==0) tp[cnt++]=hmom[i+1];
    else {
      tp[cnt]=hmom[i+1];
      DecayPoint(tp,decaylist[i],cnt);
    }
  }
  delete[] flavs;
  delete[] hmom;  
}


void Process_Base::FillOnshellConditions()
{
  if (!p_selector) return;
  int cnt=m_nin;
  vector<pair<string,double> > osc;
  p_pinfo->GetOSConditions(osc,cnt);
  for(size_t i=0;i<osc.size();i++) p_selector->AddOnshellCondition(osc[i].first,osc[i].second);  
}

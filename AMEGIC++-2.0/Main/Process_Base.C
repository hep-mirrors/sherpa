#include "Process_Base.H"
#include "Run_Parameter.H"
#include "Standard_Selector.H"
#include "Message.H"

#include "Running_AlphaS.H"
#include "Running_AlphaQED.H"

#include <algorithm>

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
  m_gen_str(3),p_b(0),m_nvec(0),m_nin(0),m_nout(0),p_fl(0),p_flin(0),p_flout(0),
  p_pl(0),p_plin(0),p_plout(0),p_moms(0),p_ps(0),p_beam(0),p_isr(0),p_cuts(0),
  p_sel(0),p_analysis(0),p_selected(0),
  m_rfactor(1.) , m_swaped(0)
{
  m_atoms=1;
  m_analyse=m_tables=0;

  m_n=m_kfactorscheme=m_scalescheme=0;
  m_nstrong=m_neweak=m_orderQCD=m_orderEW=0;
  m_totalxs=m_totalerr=m_totalsum=m_totalsumsqr=m_max=0.;
  m_last=m_lastdxs=0.;
  m_lastlumi=1.;
  m_scale=0.;
  m_isrthreshold=0.;
}



Process_Base::Process_Base(int _nin,int _nout,ATOOLS::Flavour * _fl,
			   PDF::ISR_Handler * _isr,BEAM::Beam_Spectra_Handler * _beam,
			   int _gen_str, int _orderQCD, int _orderEW,
			   int _scalescheme,int _kfactorscheme,double _scalefactor,
			   Pol_Info * _pl,
			   int _nex,ATOOLS::Flavour * _ex_fl) :
  m_nin(_nin), m_nout(_nout), m_nvec(_nin+_nout), m_nex(_nex),m_gen_str(_gen_str),
  m_orderQCD(_orderQCD), m_orderEW(_orderEW),m_nstrong(0),m_neweak(0),
  p_isr(_isr), p_beam(_beam), p_cuts(NULL), p_analysis(NULL),
  p_selected(this), p_ex_fl(_ex_fl),p_moms(NULL), 
  m_scalescheme(_scalescheme), m_kfactorscheme(_kfactorscheme), m_scalefactor(_scalefactor),
  m_n(0), m_totalxs(0.), m_totalerr(0.), m_totalsum(0.), m_totalsumsqr(0.), m_rfactor(1.),
  m_last(0.), m_lastdxs(0.), m_max(0.), m_lastlumi(1.),
  m_atoms(0), m_analyse(0), m_tables(0)
{
  p_flin    = new Flavour[m_nin];
  p_flout   = new Flavour[m_nout];  
  p_plin    = new Pol_Info[m_nin];
  p_plout   = new Pol_Info[m_nout]; 
  
  p_fl = 0;
  p_pl = 0;
  p_ps = 0;
  p_sel = 0;
  
  for (short int i=0;i<m_nin;i++) {
    p_flin[i]  = _fl[i];
    if (_pl!=0) p_plin[i] = _pl[i];
           else p_plin[i] = Pol_Info(p_flin[i]); 
    if (p_flin[i].Strong())        m_nstrong++;
    if (p_flin[i].IntCharge()!=0)  m_neweak++;
  }
  for (short int i=0;i<m_nout;i++) { 
    p_flout[i] = _fl[i+m_nin]; 
    if (_pl!=0) p_plout[i] = _pl[i+m_nin];
           else p_plout[i] = Pol_Info(p_flout[i]); 
    if (p_flout[i].Strong())       m_nstrong++;
    if (p_flout[i].IntCharge()!=0) m_neweak++;
  }
}


Process_Base::~Process_Base() {
  if (p_fl)       { delete [] p_fl;    p_fl       = 0; }
  if (p_flin)     { delete [] p_flin;  p_flin     = 0; }
  if (p_flout)    { delete [] p_flout; p_flout    = 0; }
  if (p_pl)       { delete [] p_pl;    p_pl       = 0; }
  if (p_plin)     { delete [] p_plin;  p_plin     = 0; }
  if (p_plout)    { delete [] p_plout; p_plout    = 0; }
  if (p_moms)     { delete [] p_moms;  p_moms     = 0; }
  //if (p_sel)      { delete p_sel;      p_sel      = 0; }
  if (p_cuts)     { delete p_cuts;     p_cuts     = 0; }
  if (p_ps)       { delete p_ps;       p_ps       = 0; }
  if (p_analysis) { delete p_analysis; p_analysis = 0; }
}


/*------------------------------------------------------------------------------
  
  Naming

  ------------------------------------------------------------------------------*/

string * Process_Base::GenerateNames(int _nin, Flavour * _flin, Pol_Info * _plin,
				     int _nout,Flavour * _flout,Pol_Info * _plout,
				     string & _name,string & _ptype, string & _lib)
{
  Reshuffle(_nin, _flin, _plin);
  
  if (_flin[0].IsAnti() && !_flin[1].IsAnti()) {
    Flavour flhelp  = _flin[0];
    _flin[0] = _flin[1];
    _flin[1] = flhelp;
    Pol_Info plhelp  = _plin[0];
    _plin[0] = _plin[1];
    _plin[1] = plhelp;    
  }
  
  Reshuffle(_nout, _flout, _plout);
  
  char help[20];

  sprintf(help,"%i",_nin);
  _name       = string(help);
  _name      += string("_");
  sprintf(help,"%i",_nout);
  _name      += string(help);
  _ptype      = string("P")+_name;
  _lib        = _name;
  _name      += string("_");

  short int i;

  for (i=0;i<m_nin;i++) {
    _name += string(_flin[i].Name());
    if ((_flin[i].Kfcode()==kf::e)   ||
	(_flin[i].Kfcode()==kf::mu)  ||
	(_flin[i].Kfcode()==kf::tau) ||
	(_flin[i].Kfcode()==kf::Hmin)) {
      _name.erase(_name.length()-1,1);
      if (_flin[i].IsAnti()) _name += string("+");
                        else _name += string("-");      
    }
    else {
      if (_flin[i].IsAnti()) _name += string("b"); 
    }
    // polinfo for fully polarised incommings
    if (_plin[i].p_type=='c' && _plin[i].num==1) {
      if (_plin[i].type[0]==-1) _name += string("m");
      if (_plin[i].type[0]==+1) _name += string("p");
    }
    else if (_plin[i].p_type=='h' && _plin[i].num==1) {
      if (_plin[i].type[0]==-1) _name += string("m");
      if (_plin[i].type[0]==+1) _name += string("p");
    }


    _name += string("_");
  }
  _name += string("__");

  for (i=0;i<m_nout;i++) {
    _name += string(_flout[i].Name());
    if (_flout[i].Kfcode()==kf::e  ||
	_flout[i].Kfcode()==kf::mu ||
	_flout[i].Kfcode()==kf::tau ||
	_flout[i].Kfcode()==kf::Hmin ||
	_flout[i].Kfcode()==kf::W) {
      _name.erase(_name.length()-1,1);
      if (_flout[i].IsAnti()) _name += string("+");
                         else _name += string("-");      
    }
    else {
      if (_flout[i].IsAnti()) _name += string("b"); 
    }
    // perhaps polinfo for in


    _name += string("_");
  }
  _name.erase(_name.length()-1,1);

  string hname;
  // erase _quark
  for (;;) {
    i = _name.find("_quark");
    if (i==-1) break;
    hname = _name;
    _name = hname.substr(0,i) + hname.substr(i+6); 
  }
  // gluon -> G
  for (;;) {
    i = _name.find("gluon");
    if (i==-1) break;
    hname = _name;
    _name = hname.substr(0,i) + string("G") + hname.substr(i+5); 
  }
  // jet -> j
  for (;;) {
    i = _name.find("jet");
    if (i==-1) break;
    hname = _name;
    _name = hname.substr(0,i) + string("j") + hname.substr(i+3); 
  }
  // photon -> P
  for (;;) {
    i = _name.find("photon");
    if (i==-1) break;
    hname = _name;
    _name = hname.substr(0,i) + string("P") + hname.substr(i+6); 
  }
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
    //  sort FVST
    if (a.fl.IsFermion() && !b.fl.IsFermion()) return 1;
    if (a.fl.IsVector() && !b.fl.IsFermion() && !b.fl.IsVector()) return 1;
    if (a.fl.IsScalar() && !b.fl.IsScalar() && 
	 !b.fl.IsFermion() && !b.fl.IsVector()) return 1;
    /*
    //  sort SVFT lead to problems: with jj -> W- j !!!
    if (a.fl.IsScalar() && !b.fl.IsScalar()) return 1;
    if (a.fl.IsVector() && !b.fl.IsScalar() && !b.fl.IsVector()) return 1;
    if (a.fl.IsFermion() && !b.fl.IsFermion() && 
	 !b.fl.IsScalar() && !b.fl.IsVector()) return 1;
    // redundant    if (!a.fl.IsTensor() && b.fl.Tensor()) return 1;
    */
    return 0;
  }
};

class Order_Mass {
public:
  int operator()(const S_Data & a, const S_Data & b) {
    //    if "a < b" return 0  else 1;
    if (a.fl.Mass() < b.fl.Mass()) return 0;
    return 1;
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

//     if ((a.fl.Kfcode()== b.fl.Kfcode())
// 	&& (!a.fl.IsAnti() && b.fl.IsAnti())) return 1;
//     return 0;
  }
};


class Order_Coupling {
public:
  int operator()(const S_Data & a, const S_Data & b) {
      if (!a.fl.Strong() && b.fl.Strong()) return 1;
      return 0;
  }
};

//    std::stable_sort(prea_table.begin(),prea_table.end(),Compare_Pre_Amplitudes());


void Process_Base::Reshuffle(int n, Flavour* flav, Pol_Info* plav)
{
    
  std::vector<S_Data> sd;
  for (int i=0;i<n;++i) sd.push_back(S_Data(i,flav[i],plav[i]));

  std::stable_sort(sd.begin(),sd.end(),Order_Kfc());
  std::stable_sort(sd.begin(),sd.end(),Order_Anti());
  std::stable_sort(sd.begin(),sd.end(),Order_FVST());
  std::stable_sort(sd.begin(),sd.end(),Order_Mass());
  std::stable_sort(sd.begin(),sd.end(),Order_Coupling());
  
  
  for (int i=0;i<n;++i) {
    flav[i]=sd[i].fl;
    plav[i]=sd[i].pl;
  }
  
// old method
/*
  Flavour flhelp;
  Pol_Info plhelp;
  bool hit,shuffle;

  for (;;) {
    hit = 0;
    for (short int i=0;i<n-1;i++) {
      for (short int j=i+1;j<n;j++) {
	shuffle = 0;
        if ( (flav[i].IsVector()) && !(flav[j].IsVector()))      shuffle = 1;
	else if ( (flav[i].Kfcode()) > (flav[j].Kfcode()) ) {
	  if (!( !(flav[i].IsVector()) && (flav[j].IsVector()))) shuffle = 1;
	}
	else if ( (flav[i].IsAnti()) && !(flav[j].IsAnti()) &&
		  (flav[i].Kfcode()  == flav[j].Kfcode()) )      shuffle = 1;
	
	if (shuffle) {
	  flhelp  = flav[j];
          flav[j] = flav[i];
          flav[i] = flhelp;
	  plhelp  = plav[j];
	  plav[j] = plav[i];
	  plav[i] = plhelp;
          hit = 1;
        }
      }
    }
    if (!hit) break;
  }          
*/
}

bool Process_Base::CheckExternalFlavours(int _nin,Flavour * _in,
					 int _nout,Flavour * _out) {
  // first : sum over all invariants and compare
  for (int i=0;i<_nin;i++)  { if (_in[i].Size()>1)  return 1; }
  for (int i=0;i<_nout;i++) { if (_out[i].Size()>1) return 1; }
  
  int    cin  = 0, cout  = 0;
  int    sin  = 0, sout  = 0;
  int    qin  = 0, qout  = 0;
  int    lin  = 0, lout  = 0;
  int    qfin = 0, qfout = 0;  
  int    lfin = 0, lfout = 0;  
  double bin  = 0, bout  = 0;
  for (int i=0;i<_nin;i++) {
    cin   += _in[i].IntCharge();
    sin   += _in[i].IntSpin();
    bin   += _in[i].BaryonNumber();
    lin   += _in[i].LeptonNumber();
    qin   += _in[i].StrongCharge();
    qfin  += int(pow(-1.,_in[i].IsAnti())*pow(10.,_in[i].QuarkFamily()-1));
    lfin  += int(pow(-1.,_in[i].IsAnti())*pow(10.,_in[i].LeptonFamily()-1));
  }
  for (int i=0;i<_nout;i++) {
    cout  += _out[i].IntCharge();
    sout  += _out[i].IntSpin();
    bout  += _out[i].BaryonNumber();
    lout  += _out[i].LeptonNumber();
    qout  += _out[i].StrongCharge();
    qfout += int(pow(-1.,_out[i].IsAnti())*pow(10.,_out[i].QuarkFamily()-1));
    lfout += int(pow(-1.,_out[i].IsAnti())*pow(10.,_out[i].LeptonFamily()-1));
  }
  sin = sin%2; sout = sout%2;
  qin = qin%9; qout = qout%9;

  /*
    ATOOLS::msg.Debugging()<<"Check "<<_in[0]<<" "<<_in[1]<<" -> "
    <<_out[0]<<" "<<_out[1]<<"  "
    <<cin<<" <-> "<<cout<<"  "
    <<sin<<" <-> "<<sout<<"  "
    <<bin<<" <-> "<<bout<<"  "
    <<lin<<" <-> "<<lout<<"  "
    <<qin<<" <-> "<<qout<<"  "
    <<qfin<<" <-> "<<qfout<<"  "
    <<lfin<<" <-> "<<lfout<<endl;
  */

  if (cin  != cout) return 0;    // electric charge violation
  if (sin  != sout) return 0;    // spin/fermion number violation
  if (bin  != bout) return 0;    // baryon number violation
  //if (lin  != lout) return 0;    // lepton number violation
  //if (qin  != qout) return 0;    // strong charge violation
  //if (qfin != qfout) return 0;   // quark family violation
  //if (lfin != lfout) return 0;   // lepton family violation
  return 1;
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


void Process_Base::AddChannels(Process_Base * _proc,Multi_Channel * _fsr,
			       vector<Channel_Info> & _beamparams,
			       vector<Channel_Info> & _isrparams) {

  bool         addit;
  Channel_Info ci;

  for (int i=0;i<_proc->Size();i++) {
    if ((*_proc)[i]->Partner()==NULL) AddChannels((*_proc)[i],_fsr,_beamparams,_isrparams);
    else {
      if ((*_proc)[i]->Partner()==(*_proc)[i]) {
	Single_Channel * sc;
	int next; string chname;
	for (int j=0;j<(*_proc)[i]->NumberOfFSRIntegrators();j++) { 
	  chname = ((*_proc)[i]->FSRIntegrator(j))->Name();
	  if ( (chname!=string("Rambo")) && (chname!=string("Sarge")) ) { 
	    next   = chname.find(string("--"));
	    chname = chname.substr(0,next);
	    sc     = (*_proc)[i]->PSGenerator()->SetChannel(m_nin,m_nout,p_fl,
							    ((*_proc)[i]->FSRIntegrator(j))->ChNumber(),chname);
	    sc->SetName(((*_proc)[i]->FSRIntegrator(j))->Name());
	    _fsr->Add( sc );
	  }
	}

	if ((*_proc)[i]->Beam()) {
	  if ((*_proc)[i]->Beam()->On()>0) {
	    for (int j=0;j<(*_proc)[i]->NumberOfBeamIntegrators()/3;j++) {
	      (*_proc)[i]->BeamChannels(j,ci);
	      addit = 1;
	      for (int k=0;k<_beamparams.size();k++) {
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
	      for (int k=0;k<_isrparams.size();k++) {
		if (_isrparams[k]==ci) { addit = 0; break; }
	      }
	      if (addit) _isrparams.push_back(ci);
	    }
	  }
	}
      }
    }
  }
}
















void Process_Base::InitCuts() {
  if (p_cuts == 0) {
    p_cuts = new Cut_Data();
    p_cuts->Init(m_nin+m_nout,p_fl);
  }
}


void Process_Base::UpdateCuts(double sprime,double y)
{
  p_cuts->Update(sprime,y);
  p_sel->UpdateCuts(sprime,y,p_cuts);
}





/*------------------------------------------------------------------------------

  Process management
  
  ------------------------------------------------------------------------------*/

void Process_Base::SetName(string _name)                { m_name    = _name;   }
void Process_Base::SetResDir(string _resdir)            { m_resdir  = _resdir; }
void Process_Base::SetAtoms(bool _atoms)                { m_atoms   = _atoms;  }
void Process_Base::SetTables(bool _tables)              { m_tables  = _tables; }
void Process_Base::SetBeam(Beam_Spectra_Handler * _beam){ p_beam    = _beam;   }
void Process_Base::SetISR(ISR_Handler * _isr)           { p_isr     = _isr;    }
void Process_Base::SetCuts(Cut_Data * _cuts)            { p_cuts    = _cuts;   }
void Process_Base::SetSelector(Selector_Base * _sel)    { p_sel     = _sel;    }
void Process_Base::SetMomenta(ATOOLS::Vec4D * _moms)  { p_moms    = _moms;   }
void Process_Base::SetNStrong(int _nstrong)             { m_nstrong = _nstrong;}
void Process_Base::SetNEWeak(int _neweak)               { m_neweak  = _neweak; }
void Process_Base::SetTotal(double _total)              { m_totalxs = _total;  } 
void Process_Base::SetMax(double _max)                  { m_max     = _max;    } 
void Process_Base::SetMaxJetNumber(int max)             { m_maxjetnumber  = max;    } 
void Process_Base::SetScale(double _scale)              { m_scale   = _scale;  } 
void Process_Base::SetISRThreshold(double _isrthreshold){ m_isrthreshold  = _isrthreshold;}

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


double Process_Base::Scale(ATOOLS::Vec4D * _p) {
  if (m_nin==1) return _p[0].Abs2();
  if (m_nin!=2) {
    ATOOLS::msg.Error()<<"Error in Process_Base::Scale. "
			  <<"Do not know how to handle more than 2 incoming particles."<<endl;
    abort();
  }
  double s;
  if (m_nin==2) s = (_p[0]+_p[1]).Abs2();
  if (m_nin==1) s = (_p[0]).Abs2();
  double pt2 = 0.;

  //new
  switch (m_scalescheme) {
  case 1  :
    if (m_nin+m_nout==4) {
      double t = (_p[0]-_p[2]).Abs2()-(ATOOLS::sqr(p_fl[2].PSMass())+ATOOLS::sqr(p_fl[3].PSMass()))/2.;
      double u = (_p[0]-_p[3]).Abs2()-(ATOOLS::sqr(p_fl[2].PSMass())+ATOOLS::sqr(p_fl[3].PSMass()))/2.;
      pt2 = 4.*s*t*u/(s*s+t*t+u*u);
      return pt2;
    }
    pt2 = s;
  break;
  case 2  :
    pt2 = s;
    double pt2i;
    for (int i=m_nin;i<m_nin+m_nout;i++) {
      pt2i = ATOOLS::sqr(_p[i][1])+ATOOLS::sqr(_p[i][2]);
      if (pt2i<pt2) pt2 = pt2i;
    }
    break;
  case 3 :   // swaped with 7 (to match ME.dat)
    pt2 = m_scale;
    break;
  case 4  :
    pt2 = ATOOLS::sqr(175.);
    break;
  case 5  :
    pt2 = ATOOLS::sqr(10.);
    break;
  case 6  :
    pt2 = ATOOLS::sqr(91.188);
    break;
  case 7  :
    pt2 = 0.;
    for (int i=m_nin;i<m_nin+m_nout;i++) {
      pt2 += ATOOLS::sqr(_p[i][1])+ATOOLS::sqr(_p[i][2]);
    }
    break;
  case  8 :
    pt2 = ATOOLS::sqr(3.162);
    break;
  case 42 :
    pt2 = m_scale;
    break;
  case 63 :
    if (m_nout!=m_maxjetnumber) {
     pt2 = m_scale;
    }
    else {
      pt2 = ATOOLS::sqr(_p[m_nin][1])+ATOOLS::sqr(_p[m_nin][2]);
      for (int i=m_nin+1;i<m_nin+m_nout;++i) {
	// better would be probably min yij*s with yij given by the covariant E scheme clustering algorithm
	if (p_fl[i].Strong())
	  pt2 =  ATOOLS::Min(pt2,ATOOLS::sqr(_p[i][1])+ATOOLS::sqr(_p[i][2]));
      }
    }
    break;
  case 64 :
    pt2 = ATOOLS::sqr(_p[m_nin][1])+ATOOLS::sqr(_p[m_nin][2]);
    for (int i=m_nin+1;i<m_nin+m_nout;++i) {
      // better would be min yij*s with yij given by the covariant E scheme clustering algorithm
      pt2 =  ATOOLS::Min(pt2,ATOOLS::sqr(_p[i][1])+ATOOLS::sqr(_p[i][2]));
    }
    break;
  default :
    pt2 = s;
  }
  return pt2;
}

double Process_Base::KFactor(double _scale) {
  switch (m_kfactorscheme) {
  case 2  :
    if (m_nstrong>2) {
      double f= m_rfactor*pow(as->AlphaS(m_scale * m_scalefactor)/
			 as->AlphaS(ATOOLS::sqr(ATOOLS::rpa.gen.Ecms())),m_nstrong-2);
      return f;
    } 
    else 
      return m_rfactor;
  case 1  :
    if (m_nstrong>2) {
      return m_rfactor*pow(as->AlphaS(_scale * m_scalefactor)/
			   as->AlphaS(ATOOLS::sqr(ATOOLS::rpa.gen.Ecms())),m_nstrong-2);
    } 
    else 
      return m_rfactor;
  default :
    return 1.;
  }
}

/*------------------------------------------------------------------------------
  
  Access methods
  
  ------------------------------------------------------------------------------*/

int                     Process_Base::Nin()                          { return m_nin; }
int                     Process_Base::Nout()                         { return m_nout; }
int                     Process_Base::Nvec()                         { return m_nvec; }
Flavour               * Process_Base::Flavs()                        { return p_fl; }
Vec4D                 * Process_Base::Momenta()                      { return p_moms; }
int                     Process_Base::NStrong()                      { return m_nstrong; }
int                     Process_Base::NEWeak()                       { return m_neweak; }
string                  Process_Base::Name()                         { return m_name; }
string                  Process_Base::ResDir()                       { return m_resdir; }
string                  Process_Base::LibName()                      { return string("error"); }
bool                    Process_Base::Atoms()                        { return m_atoms; }
bool                    Process_Base::Tables()                       { return m_tables; }
int                     Process_Base::NumberOfDiagrams()             { return 0; }
Point                 * Process_Base::Diagram(int i)                 { return 0; }
bool                    Process_Base::IsFreeOfFourVertex(Point * _p) { return 1; }
Beam_Spectra_Handler  * Process_Base::Beam()                         { return p_beam;     }
ISR_Handler           * Process_Base::ISR()                          { return p_isr;      }
Cut_Data              * Process_Base::Cuts()                         { return p_cuts;     }
Selector_Base         * Process_Base::Selector()                     { return p_sel;      }
Primitive_Analysis    * Process_Base::Analysis()                     { return p_analysis; }
Phase_Space_Generator * Process_Base::PSGenerator()                  { return p_psgen; }
double                  Process_Base::Scale()                        { return m_scale;    }
double                  Process_Base::Total()                        { return m_totalxs; }
double                  Process_Base::TotalError()                   { return m_totalerr; }
double                  Process_Base::Max()                          { return m_max; }
double                  Process_Base::Last()                         { return m_last; }
double                  Process_Base::LastXS()                       { return m_lastdxs; }
double                  Process_Base::LastLumi()                     { return m_lastlumi; }

double                  Process_Base::ISRThreshold()                 { return m_isrthreshold;}
int                     Process_Base::InSwaped()                     { return m_swaped;}

int                     Process_Base::ISRNumber()                                        { return 0; }
int                     Process_Base::BeamNumber()                                       { return 0; }
void                    Process_Base::ISRInfo(int,int &,double &,double &)              { return; }
void                    Process_Base::BeamInfo(int,int &,double &,double &)             { return; }

void Process_Base::BeamChannels(int i,Channel_Info & ci) { p_ps->BeamChannels(i,ci); }
void Process_Base::ISRChannels(int i,Channel_Info & ci)  { p_ps->ISRChannels(i,ci); }
int              Process_Base::NumberOfBeamIntegrators() { return p_ps->NumberOfBeamIntegrators(); }
int              Process_Base::NumberOfISRIntegrators()  { return p_ps->NumberOfISRIntegrators(); }
int              Process_Base::NumberOfFSRIntegrators()  { return p_ps->NumberOfFSRIntegrators(); }
Multi_Channel *  Process_Base::FSRIntegrator()           { return p_ps->FSRIntegrator(); }
Single_Channel * Process_Base::FSRIntegrator(int i)      { return p_ps->FSRIntegrator(i); }

void Process_Base::SwapInOrder() {
  Flavour help = p_fl[0];
  p_fl[0] = p_fl[1];
  p_fl[1] = help;
  Vec4D mom = p_moms[0];
  p_moms[0] = p_moms[1];
  p_moms[1] = mom;
  m_swaped = 1;
}

void Process_Base::RestoreInOrder() {
  if (m_swaped) {
    Flavour help = p_fl[0];
    p_fl[0] = p_fl[1];
    p_fl[1] = help;
    m_swaped = 0;
  }
}





#include "Process_Base.H"
#include "Run_Parameter.H"
#include "Standard_Selector.H"
#include "Message.H"

#include "Running_AlphaS.H"
#include "Running_AlphaQED.H"

using namespace AMEGIC;
using namespace PHASIC;
using namespace APHYTOOLS;
using namespace BEAM;
using namespace ISR;
using namespace std;


/*------------------------------------------------------------------------------
  
  Naming
  
  ------------------------------------------------------------------------------*/

string * Process_Base::GenerateNames(int _nin, Flavour * _flin, Pol_Info * _plin,
				     int _nout,Flavour * _flout,Pol_Info * _plout,
				     string & _name,string & _ptype, string & _lib)
{
  Reshuffle(_nin, _flin,  _plin);
  Reshuffle(_nout,_flout, _plout);

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

  for (i=0;i<nin;i++) {
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

  for (i=0;i<nout;i++) {
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
}

void Process_Base::Reshuffle(int n, Flavour* flav, Pol_Info* plav)
{
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
}

bool Process_Base::CheckExternalFlavours(int _nin,Flavour * _in,
					 int _nout,Flavour * _out) {
  // first : sum over all invariants and compare
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

  AORGTOOLS::msg.Debugging()<<cin<<" <-> "<<cout<<"  "
			    <<sin<<" <-> "<<sout<<"  "
			    <<bin<<" <-> "<<bout<<"  "
			    <<lin<<" <-> "<<lout<<"  "
			    <<qin<<" <-> "<<qout<<"  "
			    <<qfin<<" <-> "<<qfout<<"  "
			    <<lfin<<" <-> "<<lfout<<endl;


  if (cin  != cout) return 0;    // electric charge violation
  if (sin  != sout) return 0;    // spin/fermion number violation
  if (bin  != bout) return 0;    // baryon number violation
  //  if (lin  != lout) return 0;    // lepton number violation
  //if (qin  != qout) return 0;    // strong charge violation
  if (qfin != qfout) return 0;   // quark family violation
  //  if (lfin != lfout) return 0;   // lepton family violation
  return 1;
}

/*------------------------------------------------------------------------------
  
  Process initialization
  
  ------------------------------------------------------------------------------*/

void Process_Base::UpdateCuts(double sprime,double y)
{
  cuts->Update(sprime,y);
  sel->UpdateCuts(sprime,y,cuts);
}





/*------------------------------------------------------------------------------

  Process management
  
  ------------------------------------------------------------------------------*/

void Process_Base::SetName(string _name)               { name    = _name;   }
void Process_Base::SetResDir(string _resdir)           { resdir  = _resdir; }
void Process_Base::SetAtoms(bool _atoms)               { atoms   = _atoms;  }
void Process_Base::SetTables(bool _tables)             { tables  = _tables; }

void Process_Base::SetBeam(Beam_Handler * _beam)       { beam    = _beam;   }
void Process_Base::SetISR(ISR_Handler * _isr)          { isr     = _isr;    }
void Process_Base::SetCuts(Cut_Data * _cuts)           { cuts    = _cuts;   }
void Process_Base::SetSelector(Selector_Base * _sel)   { sel     = _sel;    }
void Process_Base::SetMomenta(AMATOOLS::Vec4D * _moms) { moms    = _moms;   }
void Process_Base::SetNStrong(int _nstrong)            { nstrong = _nstrong; }
void Process_Base::SetNEWeak(int _neweak)              { neweak  = _neweak; }

void Process_Base::SetMax(double _max)                 { max     = _max;    } 
void Process_Base::SetScale(double _scale)             { scale   = _scale;
 cout<<" new   scale="<<scale<<endl; } 

/*------------------------------------------------------------------------------

  Calculating total cross sections and single event generation

  ------------------------------------------------------------------------------*/

void Process_Base::RescaleXSec(double fac) {
  //   totalxs, totalerr, totalsum, totalsumsqr, max, rfactor;
  cout<<" in RescaleXSec of "<<name<<endl;
  cout<<"  "<<rfactor<<" "<<totalxs<<" "<<max<<" "<<endl;

  rfactor*=fac;

  totalxs*=fac;
  totalsum*=fac;
  totalerr*=fac;
  max*=fac;

  cout<<"  "<<rfactor<<" "<<totalxs<<" "<<max<<" "<<endl;

  totalsumsqr*=fac*fac; // only an estimate
}


double Process_Base::Scale(AMATOOLS::Vec4D * _p) {
  if (nin==1) return _p[0].Abs2();
  if (nin!=2) {
    AORGTOOLS::msg.Error()<<"Error in Process_Base::Scale. "
			  <<"Do not know how to handle more than 2 incoming particles."<<endl;
    abort();
  }
  double s = (_p[0]+_p[1]).Abs2();
  double t = (_p[0]-_p[2]).Abs2();
  double u = (_p[0]-_p[3]).Abs2();

  double pt2;
//    cout<<" scalescheme = "<<scalescheme<<endl;
//   cout<<" scale="<<scale<<endl;
  switch (scalescheme) {
  case 1  :
    if (nin+nout==4) pt2 = 2.*s*t*u/(s*s+t*t+u*u);
    break;
  case 2  :
    pt2 = scale;
    break;
  default :
    pt2 = s;
  }
  // cout<<" Q2="<<pt2<<endl;
  //  AORGTOOLS::msg.Out()<<"as is : "<<(*as).AlphaS(pt2)<<std::endl;
  return pt2;
}

double Process_Base::KFactor(double _scale) {
  //double ratio;
  switch (kfactorscheme) {
  case 1  :
    if (nstrong>2) {
//       double   ratio = (*APHYTOOLS::as)(_scale)/
// 	(*APHYTOOLS::as)(AMATOOLS::sqr(AORGTOOLS::rpa.gen.Ecms()));
//        double   bratio= pow(APHYTOOLS::as->AlphaS(_scale)/
// 		 APHYTOOLS::as->AsFixed(),nstrong-2);
//       AORGTOOLS::msg.Out()<<"Scale : "<<_scale<<" : "<<bratio<<" / "
//  			  <<(*APHYTOOLS::as)(_scale)<<" : "<<nstrong<<endl;
      return rfactor*pow(APHYTOOLS::as->AlphaS(_scale)/
		 APHYTOOLS::as->AsFixed(),nstrong-2);
    } 
    else {
    return rfactor;}
  default :
    return rfactor;
  }
}

/*------------------------------------------------------------------------------
  
  Access methods
  
  ------------------------------------------------------------------------------*/

inline int           Process_Base::KFactorScheme()       { return kfactorscheme; }
inline int           Process_Base::ScaleScheme()         { return scalescheme;   }

int                  Process_Base::Nin()                 { return nin; }
int                  Process_Base::Nout()                { return nout; }
int                  Process_Base::Nvec()                { return nvec; }
Flavour            * Process_Base::Flavs()               { return fl; }
Flavour            * Process_Base::FlIn()                { return flin;  }
Flavour            * Process_Base::FlOut()               { return flout; }

AMATOOLS::Vec4D    * Process_Base::Momenta()             { return moms; }

int                  Process_Base::NStrong()             { return nstrong; }
int                  Process_Base::NEWeak()              { return neweak; }

string               Process_Base::Name()                { return name; }
string               Process_Base::ResDir()              { return resdir; }
string               Process_Base::LibName()             { return string("error"); }
bool                 Process_Base::Atoms()               { return atoms; }
bool                 Process_Base::Tables()              { return tables; }

int                  Process_Base::NumberOfDiagrams()             { return 0; }
Point              * Process_Base::Diagram(int i)                 { return 0; }
bool                 Process_Base::IsFreeOfFourVertex(Point * _p) { return 1; }


Beam_Handler                  * Process_Base::Beam()     { return beam;     }
ISR_Handler                   * Process_Base::ISR()      { return isr;      }
Cut_Data                      * Process_Base::Cuts()     { return cuts;     }
Selector_Base                 * Process_Base::Selector() { return sel;      }
APHYTOOLS::Primitive_Analysis * Process_Base::Analysis() { return analysis; }
double                          Process_Base::Scale()    { return scale;    }


void Process_Base::BeamChannels(int i,Channel_Info & ci) { ps->BeamChannels(i,ci); }
void Process_Base::ISRChannels(int i,Channel_Info & ci)  { ps->ISRChannels(i,ci); }
int              Process_Base::NumberOfBeamIntegrators() { return ps->NumberOfBeamIntegrators(); }
int              Process_Base::NumberOfISRIntegrators()  { return ps->NumberOfISRIntegrators(); }
int              Process_Base::NumberOfFSRIntegrators()  { return ps->NumberOfFSRIntegrators(); }
Multi_Channel *  Process_Base::FSRIntegrator()           { return ps->FSRIntegrator(); }
Single_Channel * Process_Base::FSRIntegrator(int i)      { return ps->FSRIntegrator(i); }

double               Process_Base::Total()               { return totalxs; }
double               Process_Base::Max()                 { return max; }
double               Process_Base::Last()                { return last; }
double               Process_Base::LastXS()              { return lastdxs; }
double               Process_Base::LastLumi()            { return lastlumi; }

int  Process_Base::ISRNumber()  { return 0; }
void Process_Base::ISRInfo(int a,int & b,double & c,double & d)  { return; }


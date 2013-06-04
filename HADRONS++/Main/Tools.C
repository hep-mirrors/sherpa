#include "HADRONS++/Main/Tools.H"
#include "ATOOLS/Org/Message.H"
#include "ATOOLS/Org/MyStrStream.H"
#include "ATOOLS/Org/Run_Parameter.H"
#include "ATOOLS/Math/Algebra_Interpreter.H"

#include <vector>

//%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
//  general tools  %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
//%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

using namespace HADRONS;
using namespace ATOOLS;

PHASIC::Decay_Table * Tools::partonic_b = 
  new PHASIC::Decay_Table(Flavour(kf_b),NULL);
PHASIC::Decay_Table * Tools::partonic_c = 
  new PHASIC::Decay_Table(Flavour(kf_c),NULL);

const double Tools::Vud = 0.97377;
const double Tools::Vus = 0.2257;
const double Tools::Vub = 4.31e-3;
const double Tools::Vcd = 0.230;
const double Tools::Vcs = 0.957;
const double Tools::Vcb = 41.6e-3;
const double Tools::Vtd = 7.4e-3;
const double Tools::Vts = Tools::Vtd/0.208;
const double Tools::Vtb = 1.0;

// 3 particle phase space function lambda
double Tools::Lambda( double a, double b, double c )
{
  double L = sqr(a-b-c)-4.*b*c;
  if (L>0.) return L;
  return 0.;
}

// standard Breit Wigner with given Mass * Width
Complex Tools::BreitWigner( double s, double Mass2, double MassWidth )
{
  return Mass2/Complex(Mass2-s,-1.*MassWidth );
}

// standard Breit Wigner with given Mass * Width
Complex Tools::BreitWignerFix( double s, double Mass2, double MassWidth )
{
  return Complex(Mass2,-1.*MassWidth)/Complex(Mass2-s,-1.*MassWidth );
}

// off shell mass * width (2 particle final state with same mass)
double Tools::OffShellMassWidth( double s, double Mass2, double Width, double ms )
{
  if (s>4.*ms && Mass2>4.*ms)
    return( sqrt(s)*Width*Mass2/s * pow( (s-4.*ms)/(Mass2-4.*ms), 1.5 ) );
  return 0.;	
}

// off shell mass * width (2 particle final state with different masses)
double Tools::OffShellMassWidth( double s, double Mass2, double Width, double ms1, double ms2 )
{
  double threshold = ms1+ms2+2.*sqrt(ms1*ms2);
  if (Mass2>threshold && s>threshold)
	  return( sqrt(s)*Width*Mass2/s * pow( Mass2/s*Lambda(s,ms1,ms2)/Lambda(Mass2,ms1,ms2), 1.5 ) );
  return 0;
}

bool Tools::ExtractFlavours(std::vector<int> & helpkfc,std::string help)
{
  if (help==std::string("SPECTATORS")) return false;
  helpkfc.clear();    
  size_t pos = help.find("{");
  bool             hit;
  if (pos!=std::string::npos) help = help.substr(pos+1);
  else {
    msg_Error()<<"WARNING in "<<METHOD<<": \n"
	       <<"   Something wrong with final state of decay "
	       <<"({ - Bracket missing) : |"<<help<<"|\n"
	       <<"   Will skip it.\n";
    return false;
  }
  pos    = help.find("}");
  if (pos!=std::string::npos) help = help.substr(0,pos);
  else {
    msg_Error()<<"WARNING in "<<METHOD<<": \n"
	       <<"   Something wrong with final state of decay "
	       <<"(} - Bracket missing): |"<<help<<"|\n"
	       <<"   Will skip it.\n";
    return false;
  }
  hit    = true;
  while (hit) {
    pos      = help.find(",");
    if (pos!=std::string::npos) {
      helpkfc.push_back(atoi((help.substr(0,pos)).c_str()));
      help  = help.substr(pos+1);
    }
    else {
      helpkfc.push_back(atoi(help.c_str()));
      hit = false;
    }
  }
  if (helpkfc.size()<1) {
    msg_Error()<<"WARNING in "<<METHOD<<": \n"
           <<"   Something wrong with final state of decay. (no particles?)\n"
           <<"   Will skip it and hope for the best.\n";
    return false;
  } 
  return true;
}

bool Tools::ExtractSpecs(std::string entry,std::vector<int> & specs,
			 std::vector<double> & specweights) {
  specs.clear();
  specweights.clear();
  size_t pos = entry.find("{");
  bool             hit;
  if (pos!=std::string::npos) entry = entry.substr(pos+1);
  else {
    msg_Error()<<"WARNING in "<<METHOD<<": \n"
	       <<"   Something wrong with spectators "
	       <<"({-Bracket missing) : |"<<entry<<"|\n"
	       <<"   Will skip it.\n";
    return false;
  }
  pos    = entry.find("}");
  if (pos!=std::string::npos) entry = entry.substr(0,pos);
  else {
    msg_Error()<<"WARNING in "<<METHOD<<": \n"
	       <<"   Something wrong with spectators "
	       <<"(}-Bracket missing) : |"<<entry<<"|\n"
	       <<"   Will skip it.\n";
    return false;
  }
  hit    = true;
  std::string help;
  while (hit) {
    pos      = entry.find(",");
    if (pos!=std::string::npos) { 
      help   = entry.substr(0,pos);
      entry  = entry.substr(pos+1);
    }
    else {
      help = entry;
      hit  = false;
    }
    pos = help.find(":");
    specs.push_back(atoi((help.substr(0,pos)).c_str()));
    specweights.push_back(atof((help.substr(pos+1)).c_str()));
  }
  return (specs.size()==specweights.size() && specs.size()>0);
}

void Tools::ExtractBRInfo( std::string entry, double & br, 
			   double & dbr, std::string & origin )
{
  size_t posa, posb;        // start and end of things b/w brackets
  size_t posmin;            // start of first bracket

  std::string sbr, sdbr;

  // extract Delta BR
  posa = entry.find("(");
  posb = entry.find(")");
  posmin = posa;
  if(posa!=std::string::npos && posb!=std::string::npos)
    sdbr = entry.substr(posa+1,posb-posa-1);
  else sdbr = "-1.0";

  // extract Origin
  posa = entry.find("[");
  posb = entry.find("]");
  if(posmin==std::string::npos || 
     (posmin!=std::string::npos && posmin>posa)) posmin=posa;
  if(posa!=std::string::npos && posb!=std::string::npos)
    origin = entry.substr(posa+1,posb-posa-1);
  else origin = std::string("");

  // extract BR
  if( posmin!=std::string::npos ) sbr = entry.substr(0,posmin);
  else                            sbr = entry.substr(0);

  Algebra_Interpreter ip;
  sdbr=ip.Interprete(sdbr);
  sbr=ip.Interprete(sbr);

  dbr=ToType<double>(sdbr);
  br=ToType<double>(sbr);

  if (dbr==-1.0) dbr = br;
}

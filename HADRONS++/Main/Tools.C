#include "HADRONS++/Main/Tools.H"
#include "ATOOLS/Org/Message.H"
#include "ATOOLS/Org/MyStrStream.H"
#include "ATOOLS/Org/Run_Parameter.H"
#include "ATOOLS/Org/Data_Reader.H"
#include "ATOOLS/Math/Algebra_Interpreter.H"

#include <vector>

//%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
//  general tools  %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
//%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

using namespace HADRONS;
using namespace ATOOLS;
using namespace std;

PHASIC::Decay_Table * Tools::partonic_b = 
  new PHASIC::Decay_Table(Flavour(kf_b),NULL);
PHASIC::Decay_Table * Tools::partonic_c = 
  new PHASIC::Decay_Table(Flavour(kf_c),NULL);

const double Tools::Vud = 0.97377;
const double Tools::Vus = 0.2257;
const double Tools::Vub = 4.31e-3;
const double Tools::Vcd = 0.225;
const double Tools::Vcs = 0.973;
const double Tools::Vcb = 41.6e-3;
const double Tools::Vtd = 7.4e-3;
const double Tools::Vts = Tools::Vtd/0.208;
const double Tools::Vtb = 1.0;
const double Tools::GF = 1.16639e-5;

std::map<kf_code, kf_code> Tools::aliases;

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
  helpkfc.clear();
  if (help.find("{")!=string::npos)
    help = help.replace(help.find("{"), sizeof("{") - 1, ""); // TODO remove
  if (help.find("}")!=string::npos)
    help = help.replace(help.find("}"), sizeof("}") - 1, ""); // TODO remove
  bool hit    = true;
  while (hit) {
    size_t pos      = help.find(",");
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

#include "Spectrum_Reader.H"
#include "MyStrStream.H"
#include "Function_Base.H"
#include "MathTools.H"

using namespace BEAM;
using namespace ATOOLS;

class Polynom0 : public Function_Base {
  double m_a0;
public:
  Polynom0(const double a0) : m_a0(a0) { 
    SetDefault(a0);
  }
  double operator()(double) {
    return m_a0;
  }
};

class Polynom1 : public Function_Base {
  double m_a0, m_a1;
public:
  Polynom1(const double a0, const double a1) : m_a0(a0), m_a1(a1) { 
    SetDefault(a0);  }
  double operator()(double x) {
    return m_a0 + m_a1*x;
  }
};

class Polynom2 : public Function_Base {
  double m_a0, m_a1, m_a2;
public:
  Polynom2(const double a0, const double a1, const double a2) : 
    m_a0(a0), m_a1(a1), m_a2(a2) { 
    SetDefault(a0);  }
  double operator()(double x) {
    return m_a0 + (m_a1 + m_a2*x)*x;
  }
};

class Polynom3 : public Function_Base {
  double m_a0, m_a1, m_a2, m_a3;
public:
  Polynom3(const double a0, const double a1, const double a2, const double a3) : 
    m_a0(a0), m_a1(a1), m_a2(a2), m_a3(a3)  { 
    SetDefault(a0);  }
  double operator()(double x) {
    return m_a0 + (m_a1 + (m_a2 + m_a3*x)*x)*x;
  }
};

class Polynom4 : public Function_Base {
  double m_a0, m_a1, m_a2, m_a3, m_a4;
public:
  Polynom4(const double a0, const double a1, const double a2, const double a3, const double a4) : 
    m_a0(a0), m_a1(a1), m_a2(a2), m_a3(a3), m_a4(a4) { 
    SetDefault(a0);  }
  double operator()(double x) {
    return m_a0 + (m_a1 + (m_a2 + (m_a3 + m_a4*x)*x)*x)*x;  }
};

class Polynom5 : public Function_Base {
  double m_a0, m_a1, m_a2, m_a3, m_a4, m_a5;
public:
  Polynom5(const double a0, const double a1, const double a2, const double a3, const double a4, const double a5) : 
    m_a0(a0), m_a1(a1), m_a2(a2), m_a3(a3), m_a4(a4), m_a5(a5) { 
    SetDefault(a0);  }
  double operator()(double x) {
    return m_a0 + (m_a1 + (m_a2 + (m_a3 + (m_a4 + m_a5*x)*x)*x)*x)*x;  }
};

class Pow : public Function_Base {
  double m_a0, m_a1, m_ex;
public:
  Pow(const double a0, const double a1, const double ex) : 
    m_a0(a0), m_a1(a1), m_ex(ex) { 
    SetDefault(a0);  }
  double operator()(double x) {
    return m_a0 + m_a1*pow(x,m_ex);
  }
};

class Lin_Exp_Gauss : public Function_Base {
  double m_a0, m_a1, m_a2, m_a3, m_b2, m_sig, m_mean1, m_mean2;
public:
    Lin_Exp_Gauss(const double a0, const double a1, const double a2, const double a3,
	      const double b2, const double sig, const double mean1, const double mean2) : 
	m_a0(a0), m_a1(a1), m_a2(a2), m_a3(a3),
    m_b2(b2), m_sig(sig), m_mean1(mean1), m_mean2(mean2) { 
    SetDefault(a0);  }
  double operator()(double x) {
    return m_a0 + m_a1*x +
      m_a2*exp(m_b2*(x-m_mean1)) + m_a3*exp(-0.5*sqr((x-m_mean2)/m_sig));
  }
};

Spectrum_Reader::Spectrum_Reader(const Flavour beam, const double energy, 
   const double polarisation, const double energy_laser, 
   const double polarisation_laser,const std::string fname, const int dir) :
  Beam_Base(std::string("Spectrum_Reader"),beam,energy,polarisation,dir),
  m_fname(fname), m_energy_laser(energy_laser), m_polarisation_laser(polarisation_laser)
{
  m_weight=1.;
  m_mode=0;
  ReadFromFile();

  m_bunch        = Flavour(kf::photon);  // to be read in!
  double disc    = 1.-sqr(m_bunch.PSMass()/m_energy);
  m_vecout       = Vec4D(m_energy,0.,0.,dir*m_energy*sqrt(disc));
}

void Spectrum_Reader::ReadFromFile()
{
  m_spectrum_histo.clear();
  m_spectrum_funcs.clear();

  std::ifstream f(m_fname.c_str());

  if (!(f.good())) {
    msg.Error()<<"ERROR: Spectrum file "<<m_fname<<" not found"<<std::endl;
    return;
  }

  msg.Out()<<" reading : "<<m_fname<<" now "<<std::endl;
  double lastx=0.;
  double fbeampol=m_polarisation;
  double flaserpol=m_polarisation_laser;
  double fac=1.;

  m_upper=0.3;
  m_peak =0.3;
  double ymax=0.;
  for (;;) {
    if (!f) break;
    std::string buffer;
    getline(f,buffer);
    if (buffer.size()==0) continue;
    if (buffer[0]=='#') {
      //      msg.Out()<<buffer<<std::endl;
      if (buffer.find("parameterization")!=std::string::npos ||
	  buffer.find("parameterisation")!=std::string::npos)
	m_mode=1;
      if (buffer.find("histo")!=std::string::npos ||
	  buffer.find("histogram")!=std::string::npos)
	m_mode=0;
      if (buffer.find("beampol")!=std::string::npos) {
	  unsigned int hit=buffer.find("beampol");
	  buffer=buffer.substr(hit+8);
	  MyStrStream str;
	  str<<buffer;
	  str>>fbeampol;
      }
      if (buffer.find("laserpol")!=std::string::npos) {
	  unsigned int hit=buffer.find("laserpol");
	  buffer=buffer.substr(hit+9);
	  MyStrStream str;
	  str<<buffer;
	  str>>flaserpol;
      }
    }
    else if (m_mode==0) {
	if (fbeampol==m_polarisation &&  flaserpol==m_polarisation_laser) fac=1.;
	else if (fbeampol==-m_polarisation &&  flaserpol==-m_polarisation_laser) fac=-1.;
	else {
	    msg.Out()<<"Warning: beam spectrum file does match beam parameter!"<<std::endl;
	    msg.Out()<<" fbeampol="<<fbeampol<<std::endl;
	    msg.Out()<<"  beampol="<<m_polarisation<<std::endl;
	    msg.Out()<<" flaserpol="<<flaserpol<<std::endl;
	    msg.Out()<<"  laserpol="<<m_polarisation_laser<<std::endl;
	    fac=1.;
	}

	//      msg.Out()<<"data:"<<buffer<<std::endl;
      MyStrStream str;
      str<<buffer;
      double x1,x2,y,pol;
      // assuming mean x of bin, y value and polarisation degree
      str>>x1>>y>>pol;


      if (lastx==0. && x1!=0.) {
	lastx=x2=2.*x1;
	x1=0.;
      }
      else {
	x2=x1+(x1-lastx);
	x1=lastx;
	lastx=x2;
      } 

      if (x1>0.3 && y>ymax) {
	ymax=y;
	m_peak=x2;
      }
      if (y>0 && x2>m_upper) {
	m_upper=x2;
      }
      pol*=fac;
      //      msg.Out()<<" Data: "<<x1<<" , "<<x2<<" , "<<y<<" , "<<pol<<std::endl;
      m_spectrum_histo.push_back(Spectrum_Point(x1,x2,y,pol));
    }
    else {
	if (fbeampol==m_polarisation &&  flaserpol==m_polarisation_laser) fac=1.;
	else if (fbeampol==-m_polarisation &&  flaserpol==-m_polarisation_laser) fac=-1.;
	else {
	    msg.Out()<<"Warning: beam spectrum file does match beam parameter!"<<std::endl;
	    fac=1.;
	}
	//      msg.Out()<<"function:"<<buffer<<std::endl;
      MyStrStream str;
      str<<buffer;
      double x1, x2;
      unsigned int id;
      str>>x1>>x2;
      id=m_spectrum_funcs.size();
      m_spectrum_funcs.push_back(GetFunc(str,1.));
      m_spectrum_funcs.push_back(GetFunc(str,fac));
      m_spectrum_histo.push_back(Spectrum_Point(x1,x2,id));

      double x=0.5*(x1 + x2);
      double y  =(*m_spectrum_funcs[id])(x1);
      double pol=(*m_spectrum_funcs[id+1])(x1);
      if (x1>0.3 && y>ymax) {
	ymax=y;
	m_peak=x1;
      }
      //      msg.Out()<<" Function: "<<x1<<" , "<<y<<" , "<<pol<<std::endl;
      y  =(*m_spectrum_funcs[id])(x);
      pol=(*m_spectrum_funcs[id+1])(x);
      if (y>0 && x2>m_upper) {
	m_upper=x2;
      }

      if (x1>0.3 && y>ymax) {
	ymax=y;
	m_peak=x;
      }
      //      msg.Out()<<" Function: "<<x<<" , "<<y<<" , "<<pol<<std::endl;
      y  =(*m_spectrum_funcs[id])(x2);
      pol=(*m_spectrum_funcs[id+1])(x2);
      if (x1>0.3 && y>ymax) {
	ymax=y;
	m_peak=x2;
      }
      //      msg.Out()<<" Function: "<<x2<<" , "<<y<<" , "<<pol<<std::endl;
    }
  }

  msg.Out()<<" m_upper="<<m_upper<<std::endl;
  msg.Out()<<" m_peak ="<<m_peak<<std::endl;

  // --- this is a test ---
  PrintSpectra("tspec");
  //  exit(0);
}


Beam_Base * Spectrum_Reader::Copy()
{
  return new Spectrum_Reader(m_beam,m_energy,m_polarisation,m_energy_laser,m_polarisation_laser,m_fname,m_dir);
}

ATOOLS::Flavour Spectrum_Reader::Remnant() {
  return m_beam;
}

bool Spectrum_Reader::CalculateWeight(const double x,const double scale)
{
  if (x<0. || 1.<x) {
    msg.Out()<<" Error: x out of range! ("<<x<<")"<<std::endl;
    return 0.;
  }

  if (m_mode==0) {
    for (Spectrum_Data_List::const_iterator sd=m_spectrum_histo.begin();
      sd!=m_spectrum_histo.end();++sd) {
      if (sd->X1()<=x && x<=sd->X2()) {
	  m_weight=sd->Y();
	  m_polar =sd->Pol();
	  if (m_polar>1.) m_polar=1.;
	  if (m_polar<-1.) m_polar=-1.;
	  if (m_weight<0.) m_weight=0.;
	  return true;
      }
    }
  }
  else {
    for (Spectrum_Data_List::const_iterator sd=m_spectrum_histo.begin();
      sd!=m_spectrum_histo.end();++sd) {
      if (sd->X1()<=x && x<=sd->X2()) {
	  unsigned int id=sd->Id();
	  m_weight=(*m_spectrum_funcs[id])(x);
	  m_polar =(*m_spectrum_funcs[id+1])(x);
	  if (m_polar>1.) m_polar=1.;
	  if (m_polar<-1.) m_polar=-1.;
	  if (m_weight<0.) m_weight=0.;
	  return true;
      }
    }
  }

  msg.Out()<<" Error: x="<<x<<" not found in data points"<<std::endl;
  m_weight=1.;
  return false;
}

double Spectrum_Reader::Weight(Flavour flav)
{
  return m_weight;
}

void Spectrum_Reader::PrintSpectra(std::string name,int) {
  MyStrStream str;
  str<<name<<m_dir<<".dat";
  str>>name;
  
  std::ofstream f(name.c_str());
  for (double x=0.001; x<=1.; x+=0.003) {
    CalculateWeight(x,1.);
    double y=Weight(Flavour(kf::photon));
    double p=Polarisation();
    f<<x<<"  "<<y<<"  "<<p<<std::endl;
  }
  f.close();
}


Function_Base * Spectrum_Reader::GetFunc(MyStrStream & str, const double fac) 
{
  std::string name;
  str>>name;
  if (name=="Pol0") {
    double a0;
    str>>a0;
    return new Polynom0(a0*fac);
  }
  if (name=="Pol1") {
    double a0,a1;
    str>>a0>>a1;
    return new Polynom1(a0*fac,a1*fac);
  }
  if (name=="Pol2") {
    double a0,a1,a2;
    str>>a0>>a1>>a2;
    return new Polynom2(a0*fac,a1*fac,a2*fac);
  }
  if (name=="Pol3") {
    double a0,a1,a2,a3;
    str>>a0>>a1>>a2>>a3;
    return new Polynom3(a0*fac,a1*fac,a2*fac,a3*fac);
  }
  if (name=="Pol4") {
    double a0,a1,a2,a3,a4;
    str>>a0>>a1>>a2>>a3>>a4;
    return new Polynom4(a0*fac,a1*fac,a2*fac,a3*fac,a4*fac);
  }
  if (name=="Pol5") {
    double a0,a1,a2,a3,a4,a5;
    str>>a0>>a1>>a2>>a3>>a4>>a5;
    return new Polynom5(a0*fac,a1*fac,a2*fac,a3*fac,a4*fac,a5*fac);
  }
  if (name=="Pow") {
      double a0,a1,ex;
      str>>a0>>a1>>ex;
      return new Pow(a0*fac,a1*fac,ex);
  }
  if (name=="LEG") {
      double a0, a1, a2, a3, b2, sig, mean1, mean2;
      str>>a0>>a1>>a2>>a3>>b2>>sig>>mean1>>mean2;
      return new Lin_Exp_Gauss(a0*fac,a1*fac,a2*fac,a3*fac,b2,sig,mean1,mean2);
  }

  msg.Out()<<" ERROR: unknown function type in Spectrum_Reader : "<<name<<std::endl;

  return new Polynom0(0.);
}

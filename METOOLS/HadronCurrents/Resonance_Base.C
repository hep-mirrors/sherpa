#include "METOOLS/HadronCurrents/Resonance_Base.H"
#include "ATOOLS/Org/Message.H"
#include "ATOOLS/Org/MyStrStream.H"

using namespace METOOLS;
using namespace ATOOLS;
using namespace std;

/////////////////////////////////////////////////////////////
//
// Partial width base class
//
/////////////////////////////////////////////////////////////

Partial_Width_Base::
Partial_Width_Base(const Flavour & inflav,
		   const vector<Flavour> & outflavs,
		   const double & BR) :
  m_inflav(inflav), m_outflavs(outflavs), m_BR(BR), 
  m_mass(inflav.HadMass()), m_mass2(m_mass*m_mass),
  m_mmin(0.), m_smin(0.), 
  m_width(inflav.Width()),  m_width2(m_width*m_width),
  m_partialwidth(m_BR*m_width),
  m_prefactor(1.)
{
  for (size_t i=0;i<m_outflavs.size();i++) {
    double mass = m_outflavs[i].HadMass();
    m_decmasses.push_back(mass);
    m_decmasses2.push_back(mass*mass);
    m_mmin += mass;
  }
  m_smin = sqr(m_mmin);
  msg_Out()<<METHOD<<" for "<<m_inflav<<" -->";
  for (size_t i=0;i<m_outflavs.size();i++) msg_Out()<<" "<<m_outflavs[i];
  msg_Out()<<".\n";
}

const double Partial_Width_Base::operator()(const double & s) {
  if (s<m_smin) return 0.;
  return Calculate(s);
}

void Partial_Width_Base::FixPrefactor() {
  double raw_width = (*this)(m_mass2);
  m_prefactor = m_partialwidth / raw_width;
}

void Partial_Width_Base::Output() {
  double gamma = (*this)(m_mass2);
  msg_Out()<<"  * [";
  for (size_t i=0;i<m_outflavs.size();i++) msg_Out()<<" "<<m_outflavs[i];
  msg_Out()<<", BR = "<<m_BR<<"]: "
	   <<(100.*gamma/m_width)<<"%.\n";
}


/////////////////////////////////////////////////////////////
//
// Total width base class
//
/////////////////////////////////////////////////////////////

Total_Width_Base::Total_Width_Base(const ATOOLS::Flavour & inflav) :
  m_inflav(inflav),
  m_mass(inflav.HadMass()), m_mass2(m_mass*m_mass),
  m_width(inflav.Width()),  m_width2(m_width*m_width) {}

Total_Width_Base::~Total_Width_Base() {
  while (!m_channels.empty()) {
    delete *m_channels.begin();
    m_channels.erase(m_channels.begin());
  }
}

const double Total_Width_Base::operator()(const double & s) {
  double width = (m_channels.size()==0 ? m_width : 0.);
  for (set<Partial_Width_Base *>::iterator pwit=m_channels.begin();
       pwit!=m_channels.end();pwit++) {
    width += (**pwit)(s);
  }
  return width;
}

void Total_Width_Base::PrintBRs() const {
  msg_Out()<<"======== "<<METHOD<<" ======= \n"
	   <<"Checking "<<m_inflav<<" decays:\n";
  for (set<Partial_Width_Base *>::iterator pwit=m_channels.begin();
       pwit!=m_channels.end();pwit++) {
    (*pwit)->Output();
  }
}

void Total_Width_Base::OutputLineshape(const double & mmin,const double &mmax,
				       const size_t & steps) {
  ofstream output;
  string   filename = string(m_inflav.IDName()+"_Lineshape.dat");
  output.open(filename.c_str());
  double stepsize = (mmax-mmin)/double(steps);
  for (size_t i=0;i<=steps;i++) {
    double Q = mmin+i*stepsize;
    double Gamma = (*this)(Q*Q);
    output<<std::setprecision(6)<<std::setw(12)<<Q<<"   "
	  <<std::setprecision(6)<<std::setw(12)<<Gamma<<"\n";
    //msg_Out()<<std::setprecision(6)<<std::setw(12)<<Q<<"   "
    //	     <<std::setprecision(6)<<std::setw(12)<<Gamma<<"\n";
  }
  output.close();
}


Partial_Width_Base * Total_Width_Base::SelectChannel(const double & s) {
  return NULL;
}





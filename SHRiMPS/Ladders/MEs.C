#include "SHRiMPS/Ladders/MEs.H"
#include "ATOOLS/Org/Run_Parameter.H"

using namespace SHRIMPS;
using namespace ATOOLS;

MEs::MEs(const double & smin, const double & tmin) :
  p_sigma(NULL), m_shatmin(smin), m_thatmin(tmin)
{ }


void MEs::SetPartonic(Sigma_Partonic * sigma) {
  p_sigma = sigma; 
  if (m_shatmin<0.)  m_shatmin = p_sigma->Smin(); 
  if (m_thatmin==0.) m_thatmin = p_sigma->Tmin(); 
}
  
double MEs::PDFratio(const Vec4D & qprev,const ATOOLS::Flavour & fprev,
		     const Vec4D & qact,const ATOOLS::Flavour & fact,
		     const size_t & dir) {
  double xprev, xact;
  if (dir==0) {
    xprev = qprev.PPlus()/rpa->gen.PBeam(0).PPlus();
    xact  = qact.PPlus()/rpa->gen.PBeam(0).PPlus();
  }
  else {
    xprev = qprev.PMinus()/rpa->gen.PBeam(1).PMinus();
    xact  = qact.PMinus()/rpa->gen.PBeam(1).PMinus();
  }
  if (xprev<0. || xact<0.) return 0.;
  return (p_sigma->PDF(dir, xact,  dabs(qact.Abs2()),  fact)/
	  p_sigma->PDF(dir, xprev, dabs(qprev.Abs2()), fprev));
}

double MEs::operator()(Ladder * ladder) {
  TPropList::iterator winner;
  if (!ladder->ExtractHardest(winner,0.) ||
      dabs(winner->Q2())<dabs(m_thatmin)) {
    return 1.;
  }
  return (winner->Col()==colour_type::octet ?
	  m_thatmin/dabs(winner->Q2()) :
	  sqr(m_thatmin/dabs(winner->Q2())) );
}

double MEs::operator()(Ladder * ladder, const double & qt2min) {
  if (ladder->InPart(0)->Momentum().PPlus()>rpa->gen.Ecms() ||
      ladder->InPart(1)->Momentum().PMinus()>rpa->gen.Ecms()) return 0.;
  if (ladder->GetProps()->size()==1)                          return 1.;
  
  TPropList::iterator winner;
  if (!ladder->ExtractHardest(winner,qt2min) ||
      dabs(winner->Q2())<dabs(m_thatmin))                     return 1.;
  double weight = (winner->Col()==colour_type::octet ?
		   m_thatmin/dabs(winner->Q2()) :
		   sqr(m_thatmin/dabs(winner->Q2())) );
  
  Vec4D q[2];
  ladder->HardestIncomingMomenta(winner, q[0], q[1]);
  for (size_t i=0;i<2;i++) {
    double x =  (i==0 ?
		 q[i].PPlus()/rpa->gen.PBeam(i).PPlus():
		 q[i].PMinus()/rpa->gen.PBeam(i).PMinus() );
    weight *= p_sigma->PDF(0,x,winner->QT2());
  }
  return weight;
}


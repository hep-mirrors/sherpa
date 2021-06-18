#include "SHRiMPS/Ladders/MEs.H"
#include "ATOOLS/Org/Run_Parameter.H"

using namespace SHRIMPS;
using namespace ATOOLS;

MEs::MEs(Sigma_Partonic * sigma,const double & smin, const double & tmin) :
  p_sigma(sigma), m_shatmin(smin), m_thatmin(tmin)
{
  if (m_shatmin<0.)  m_shatmin = p_sigma->Smin(); 
  if (m_thatmin==0.) m_thatmin = p_sigma->Tmin(); 
}

double MEs::operator()(Ladder * ladder, const double & qt2min) {
  if (ladder->InPart(0)->Momentum().PPlus()>rpa->gen.Ecms() ||
      ladder->InPart(1)->Momentum().PMinus()>rpa->gen.Ecms()) return 0.;
  if (ladder->GetProps()->size()==1)                          return 1.;
  
  TPropList::iterator winner;
  if (!ladder->ExtractHardest(winner,qt2min) ||
      dabs(winner->Q2())<dabs(m_thatmin))                     return 1.;
  
  Vec4D q[2];
  ladder->HardestIncomingMomenta(winner, q[0], q[1]);
  if (q[0].PPlus()>rpa->gen.Ecms() ||
      q[1].PMinus()>rpa->gen.Ecms() ||
      (q[0]+q[1]).Abs2()<m_shatmin)                           return 0.;
  //msg_Out()<<METHOD<<" for qt^2 = "<<qt2max<<" vs. "<<win->Q().Abs2()<<"\n";
  double weight = (winner->Col()==colour_type::octet?
		   qt2min/dabs(winner->Q2()) :
		   sqr(qt2min/dabs(winner->Q2())));

  for (size_t i=0;i<2;i++) {
    double x0 = (i==0?
		 ladder->InPart(i)->Momentum().PPlus()/rpa->gen.PBeam(i).PPlus():
		 ladder->InPart(i)->Momentum().PMinus()/rpa->gen.PBeam(i).PMinus());
    double x1 = (i==0?
		 q[i].PPlus()/rpa->gen.PBeam(i).PPlus():
		 q[i].PMinus()/rpa->gen.PBeam(i).PMinus());
    double pdf0 = p_sigma->PDF(i,x0,0.);
    double pdf1 = p_sigma->PDF(i,x1,winner->QT2());
    if (pdf1/pdf0<0.) {
      msg_Out()<<METHOD<<" pushes unphysical result: "<<(pdf1/pdf0)<<" "
	       <<"for x1 = "<<x1<<" & x0 = "<<x0<<"\n";
    }
    weight *= (pdf1>0.?pdf1/pdf0:0.);
  }
  return weight;
}


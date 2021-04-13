#include "SHRiMPS/Ladders/Ladder_Generator_LDC.H"
#include "SHRiMPS/Tools/MinBias_Parameters.H"
#include "SHRiMPS/Tools/Special_Functions.H"
#include "MODEL/Main/Model_Base.H"
#include "ATOOLS/Math/Random.H"
#include "ATOOLS/Org/Run_Parameter.H"
#include <list>

using namespace SHRIMPS;
using namespace MODEL;
using namespace ATOOLS;
using namespace std;

Ladder_Generator_LDC::Ladder_Generator_LDC() :
  Ladder_Generator_Base(),
  m_partonic(Sigma_Partonic(xs_mode::Regge)),
  m_kt2min(m_Q02)
{
  for (size_t i=0;i<2;i++) {
    m_Ebeam[i] = rpa->gen.PBeam(i)[0];
    m_Pbeam[i] = rpa->gen.PBeam(i);
  }
}

Ladder * Ladder_Generator_LDC::operator()(const Vec4D & pos) {
  m_weight = 0.;
  p_ladder = new Ladder(pos);
  if (FixInitialPartons() &&
      FillZs()) {
    exit(1);
  }
  else { delete p_ladder; p_ladder = NULL; }
  exit(1);
  return p_ladder;
}

bool Ladder_Generator_LDC::FixInitialPartons() {
  do {
    for (size_t beam=0;beam<2;beam++) {
      m_y[beam]   = (beam==0 ? 1.: -1.) * (m_Ymax + (1.-2.*ran->Get()*m_deltaY));
      m_kt2[beam] = p_eikonal->FF(beam)->SelectQT2(sqr(m_Ebeam[beam]/cosh(m_y[beam])));
      double phi  = 2.*M_PI*ran->Get();
      m_k[beam]   = sqrt(m_kt2[beam]) * Vec4D(cosh(m_y[beam]),cos(phi),sin(phi),sinh(m_y[beam]));
      m_zp[beam]  = m_k[beam].PPlus()/m_Pbeam[0].PPlus();
      m_zm[beam]  = m_k[beam].PMinus()/m_Pbeam[1].PMinus();;
    }
  } while (m_k[0][0]>m_Ebeam[0] || m_k[1][0]>m_Ebeam[1]);
  for (size_t beam=0;beam<2;beam++) {
    m_zps.insert(m_zp[beam]);
    m_zms.insert(m_zm[beam]);    
  }
  Output();
  return true;
}

bool Ladder_Generator_LDC::FillZs() {
  double a = sqrt(3.*AlphaS(m_kt2[0])/M_PI) * log(1./m_zp[1]); 
  double b = sqrt(3.*AlphaS(m_kt2[0])/M_PI) * log(1./m_zm[0]);
  //double a = sqrt(3.*AlphaS(m_kt2[0])/M_PI) * log(m_zp[0]/m_zp[1]); 
  //double b = sqrt(3.*AlphaS(m_kt2[0])/M_PI) * log(m_zm[1]/m_zm[0]);
  double zp, zm;
  size_t N = SelectN(a,b);
  for (size_t i=0;i<N;i++) {
    do {
      zp = a*ran->Get();
      zm = b*ran->Get();
      msg_Out()<<"z+/z- = "<<(zp/zm)<<" in ["<<(m_zp[1]/m_zm[0])<<", "<<(m_zp[0]/m_zm[1])<<"]\n";
    } while (zp/zm < m_zp[1]/m_zm[0] || zp>m_zp[0]/m_zm[1]);
    m_zps.insert(zp);
    m_zms.insert(zm);
  }
  Output();
  return true;
}

size_t Ladder_Generator_LDC::SelectN(const double & a,const double & b) {
  double G_over = sqrt(a/b) * SF.In(1,2.*sqrt(a*b)), f = sqrt(a*b)/G_over, summed = f;
  double random = ran->Get();
  size_t n      = 1;
  msg_Out()<<METHOD<<"(a = "<<a<<" and b = "<<b<<"): f = "<<f<<", G = "<<G_over<<", "
	   <<"random = "<<random<<".\n";
  while (random>summed && n<20) {
    f      *= (a*b)/(n+1)/n++;
    summed += f;
    msg_Out()<<"   f = "<<f<<", summed = "<<summed<<"\n";
    //n++;
  }
  msg_Out()<<METHOD<<" yields "<<(n-1)<<" emissions.\n";
  return n-1;
}


void Ladder_Generator_LDC::CalculateWeight() {
  m_weight = 1.;
}

void Ladder_Generator_LDC::Output() {
  msg_Out()<<METHOD<<":\n";
  for (size_t beam=0;beam<2;beam++) {
    msg_Out()<<"   * "<<m_k[beam]<<" : kt^2 = "<<m_k[beam].PPerp2()<<", y = "<<m_y[beam]<<", "
	     <<"z+ = "<<m_zp[beam]<<", z- = "<<m_zm[beam]<<".\n";
  }
  msg_Out()<<"z+: ";
  for (set<double>::iterator sit=m_zps.begin();sit!=m_zps.end();sit++) msg_Out()<<(*sit)<<" ";
  msg_Out()<<"\n";
  msg_Out()<<"z-: ";
  for (set<double>::iterator sit=m_zms.begin();sit!=m_zms.end();sit++) msg_Out()<<(*sit)<<" ";
  msg_Out()<<"\n";
  msg_Out()<<"y:  ";
  set<double>::reverse_iterator spit=m_zps.rbegin();
  set<double>::iterator smit=m_zms.begin();
  while (spit!=m_zps.rend() && smit!=m_zms.end()) {
    msg_Out()<<(1./2.*log((*spit)/(*smit)))<<" ";
    spit++; smit++;
  }
  msg_Out()<<"\n";
}

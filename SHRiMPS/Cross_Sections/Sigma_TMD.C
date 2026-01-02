#include "SHRiMPS/Cross_Sections/Sigma_TMD.H"
#include "PHASIC++/Channels/Channel_Elements.H"
#include "ATOOLS/Math/Random.H"
#include "ATOOLS/Org/Run_Parameter.H"
#include "ATOOLS/Org/Message.H"


using namespace SHRIMPS;
using namespace ATOOLS;

Sigma_TMD::Sigma_TMD(const double & pt02) :
  m_xsec(0.), m_dxsmax(0.), 
  m_Ymax(MBpars.GetEikonalParameters().originalY),
  m_S(sqr(rpa->gen.Ecms())),
  m_pt02(pt02),m_pt2max_norm(m_S/(4.*m_pt02)+1.),
  m_xmin(0.2/rpa->gen.Ecms())
{
  m_P[0] = sqrt(m_S/4.) * Vec4D(1.,0.,0.,1.);
  m_P[1] = sqrt(m_S/4.) * Vec4D(1.,0.,0.,-1.);
}

Sigma_TMD::~Sigma_TMD() {}
    
void Sigma_TMD::Initialise(REMNANTS::Pseudo_TMD * tmd[2]) {
  for (size_t i=0;i<2;i++) p_tmd[i] = tmd[i];
  Calculate();
}

const double Sigma_TMD::MakeEvent() {
  double dxs;
  do {
    dxs = MakePoint()*MakeMEWeight();
  } while (dxs<m_dxsmax*ran->Get());
  //msg_Out()<<METHOD<<" y = "<<m_y<<" in ["<<m_ymin<<", "<<m_ymax<<", "
  //	   <<"pt = "<<m_pt2[0]<<", "<<m_pt2[1]<<"] "
  //	   <<"from tau = "<<m_tau<<", xmin = "<<m_xmin<<" -->"<<(dxs/m_dxsmax)<<"\n";
  return m_xsec;
}

const double Sigma_TMD::xTMD(const size_t beam,
			    const double & x,const double & Q2,
			    const double & pt2,const double & Y,
			    const ATOOLS::Flavour & flav) {
  p_tmd[beam]->Calculate(x,Q2,pt2,Y);
  m_xtmd[beam] = p_tmd[beam]->XTMD(flav); 
  return m_xtmd[beam];
}

const double Sigma_TMD::MakePoint() {
  double dummy, expo = 2., weight=0.;
  for (size_t i=0;i<2;i++) {
    m_pt2[i]  = PHASIC::PeakedDist(m_pt02,expo,0.,m_S/4.,1,ran->Get());
    m_phi[i]  = 2.*M_PI*ran->Get();
    m_qvec[i] = sqrt(m_pt2[i])*Vec4D(0.,cos(m_phi[i]),sin(m_phi[i]),0.);
  }
  m_tau  = (m_qvec[0]+m_qvec[1]).PPerp2()/m_S;
  if (m_tau>sqr(m_xmin)) {
    m_ymax = log(1./m_tau)/2.;
    m_ymin = -m_ymax;
    m_y    = m_ymin + ran->Get()*(m_ymax-m_ymin);
    m_x[0] = sqrt(m_tau*exp(m_y));
    m_x[1] = sqrt(m_tau*exp(-m_y));
    if (m_x[0]>1.-m_xmin || m_x[0]<m_xmin ||
	m_x[1]>1.-m_xmin || m_x[1]<m_xmin) return 0.;
    for (size_t i=0;i<2;i++) m_qvec[i] += m_x[i]*m_P[i];
    m_kvec = m_qvec[0] + m_qvec[1];
    // Jacobeans of two azimuthal angles, one rapidity, and two transverse momentum
    // integrations, times factors 1/4 from the transformation d^2pt -> dpt^2,
    // 2 from dx1 dx2 -> dy dtau, and 1/4E^2=1/S) from the delta function in k^2
    weight = ( sqr(2.*M_PI) * (m_ymax-m_ymin) *
	       PHASIC::PeakedWeight(m_pt02,expo,0.,m_S/4.,m_pt2[0],1,dummy)*
	       PHASIC::PeakedWeight(m_pt02,expo,0.,m_S/4.,m_pt2[1],1,dummy)/
	       (4.*m_S) );
    //msg_Out()<<METHOD<<"[y = "<<m_y<<", qt = "<<m_pt2[0]<<", "<<m_pt2[1]<<"] "
    //	     <<"yields wt = "<<weight<<"\n";
  }
  return weight;
}

const double Sigma_TMD::MakeMEWeight() {
  //double disc  = sqr(2.*m_tau*m_S)-m_pt2[0]*m_pt2[1];
  double disc  = sqr(m_qvec[0]*m_qvec[1])-m_qvec[0].Abs2()*m_qvec[1].Abs2();
  if (disc<0.) return 0.;
  double flux  = sqrt(disc), colfac = 3./8., tmds = 1.;
  double scale = dabs((m_qvec[0]-m_x[0]*m_P[0]+m_qvec[1]-m_x[1]*m_P[1]).Abs2());
  double alpha = 4.*M_PI*(*p_alphaS)(m_kvec.PPerp2());
  for (size_t beam=0;beam<2;beam++) {
    double tmd = xTMD(beam,m_x[beam],scale,m_pt2[beam],m_y);
    tmds      *= tmd;
  }
  double me2    = gg2g();
  double weight = colfac * alpha * tmds/m_tau * me2 / flux;
  //if (weight>1.)
  //msg_Out()<<METHOD<<"[y = "<<m_y<<"] yields wt = "
  //<<(tmds/m_tau)<<" * "<<me2<<" / "<<flux<<" = "<<weight<<"\n";
  return weight;
}

const double Sigma_TMD::gg2g() {
  double me2 = (sqr(m_tau*m_S/2.)/((m_pt2[0]+m_pt02)*(m_pt2[1]+m_pt02)) *
		(2.*m_tau*m_S + 4.*m_pt2[0]*m_pt2[1]/(m_tau*m_S) +
		 (m_qvec[0]+m_qvec[1]).Abs2()) );
  return me2;
}


const bool Sigma_TMD::Calculate() {
  double PS, ME, sum = 0., sigma = 0., sum2 = 0., error = 1.;
  long int N = 0., maxiter = 5000;
  while (error>1.e-2 && N<100000) {
    for (size_t n=0;n<maxiter;n++) {
      PS    = MakePoint();
      ME    = (PS>1.e-12 ? MakeMEWeight() : 0.);
      sum  += sigma = PS*ME;
      if (sigma>m_dxsmax) m_dxsmax = sigma;
      sum2 += sqr(sigma);
    }
    N    += maxiter;
    error = sqrt(sum2-sqr(sum/double(N)))/sum;
  }
  m_xsec = sum/double(N)*rpa->Picobarn();
  msg_Out()<<"=======================================================================\n"
	   <<"=======================================================================\n"
	   <<METHOD<<"(n = "<<N<<"): "
	   <<"sigma = "<<m_xsec<<" pb +/- "<<(100.*error)<<"%, max = "<<m_dxsmax<<".\n"
	   <<"=======================================================================\n"
	   <<"=======================================================================\n";
  return true;
}

void Sigma_TMD::SelectFlavours(const bool & fixflavour) {
}

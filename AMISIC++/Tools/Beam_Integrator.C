#include "AMISIC++/Tools/NonPerturbative_XSecs.H"

using namespace AMISIC;
using namespace BEAM;
using namespace ATOOLS;
using namespace std;

void Beam_Integrator::Init(Hadronic_XSec_Calculator * xsecs,
			   Beam_Base * beam1,Beam_Base * beam2,
			   evt_type::code & type) {
  p_xsecs = xsecs; p_beams[0] = beam1; p_beams[1] = beam2;
  for (size_t i=0;i<2;i++) {
    m_variable   += (p_beams[i]->Type()==beamspectrum::EPA)*(i+1);
    m_beammoms[i] = p_beams[i]->InMomentum();
  }
  m_smin    = Max(m_smin,1.0001*xsecs->Smin()); 
  m_smax    = (m_beammoms[0]+m_beammoms[1]).Abs2();
  m_swt     = log(m_smax/m_smin);
  m_sign    = m_beammoms[0][3]>0. ? 1: -1;
  m_evttype = type;
  CalculateXSecs();
}

bool Beam_Integrator::operator()() {
  long int trials = 0;
  do {
    if (trials++>1000000) {
      msg_Error()<<"Error in "<<METHOD<<": no good event kinematics found.\n";
      return false;
    }
  } while (TrialEvent() < ran->Get()*m_max);
  for (size_t beam=0;beam<2;beam++) p_beams[beam]->SetOutMomentum(m_inmoms[beam]);
  return true;
}

double Beam_Integrator::TrialEvent() {
  double weight = MakePoint(), meweight;
  double x[2];
  for (size_t beam=0;beam<2;beam++) {
    x[beam] = m_inmoms[beam][0]/m_beammoms[beam][0]; 
    p_beams[beam]->CalculateWeight(x[beam],0.);
  }
  weight *= ( m_weight = (x[0]*x[1]*p_beams[0]->Weight()*p_beams[1]->Weight() ));
  (*p_xsecs)(m_sprime);
  weight   *= meweight = MEweight();
  m_weight *= meweight;
  return weight;
}

void Beam_Integrator::CalculateXSecs() {
  size_t cycles = 0;
  double sum    = 0.,      sum2 = 0.;
  msg_Debugging()<<METHOD<<" will integrate over s and y in E range "
		 <<"["<<sqrt(m_smin)<<", "<<sqrt(m_smax)<<"]\n";
  do {
    for (size_t i=0;i<m_Npoints;i++) {
      double wt = TrialEvent();
      sum  += wt; sum2 += wt*wt;
      if (wt>m_max) m_max = wt;
    }
    cycles++;
    m_total = sum/double(cycles*m_Npoints);
    m_error = sqrt(sum2/double(cycles*m_Npoints)-sqr(m_total))/double(cycles*m_Npoints)/m_total;
    if (cycles>100) THROW(fatal_error,"No convergence achieved.")
  } while (m_error>1.e-3);
  msg_Info()<<METHOD<<" converges after "<<double(cycles*m_Npoints)<<" PS points with\n"
	    <<"    * xs = "<<m_total<<" mb +/- "<<(m_error*100.)<<" %, max = "<<m_max<<".\n";
}

double Beam_Integrator::MakePoint() {
  m_sprime   = m_smin * pow(m_smax/m_smin,ran->Get());
  double ywt = 1., logtau;
  switch (m_variable) {
  case 3:
    logtau      = ywt = log(m_smax/m_sprime);
    m_yprime    = logtau * (ran->Get()-0.5);
    m_inmoms[0] = sqrt(m_sprime)/2. * exp( m_yprime) * Vec4D(1.,0,0., m_sign);
    m_inmoms[1] = sqrt(m_sprime)/2. * exp(-m_yprime) * Vec4D(1.,0,0.,-m_sign);
    break;
  case 2:
    m_inmoms[0] = m_beammoms[0];
    m_inmoms[1] = m_sprime/(4.*m_beammoms[0][0]) * Vec4D(1.,0.,0.,-m_sign);
    m_yprime    = 0.5 * log((m_inmoms[0][0]+m_inmoms[0][3] + m_inmoms[1][0]+m_inmoms[1][3])/
			    (m_inmoms[0][0]-m_inmoms[0][3] + m_inmoms[1][0]-m_inmoms[1][3]));
    break;
  case 1:
    m_inmoms[0] = m_sprime/(4.*m_beammoms[1][0]) * Vec4D(1.,0.,0.,m_sign);
    m_inmoms[1] = m_beammoms[1];
    m_yprime    = 0.5 * log((m_inmoms[0][0]+m_inmoms[0][3] + m_inmoms[1][0]+m_inmoms[1][3])/
			    (m_inmoms[0][0]-m_inmoms[0][3] + m_inmoms[1][0]-m_inmoms[1][3]));
    break;
  case 0:
  default:
    msg_Error()<<METHOD<<" should not arrive here!\n";
    exit(1);
  }
  return m_swt*ywt;
}

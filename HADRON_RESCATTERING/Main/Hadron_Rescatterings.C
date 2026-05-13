#include "HADRON_RESCATTERING/Main/Hadron_Rescatterings.H"
#include "ATOOLS/Math/Random.H"
#include "ATOOLS/Math/Poincare.H"
using namespace HADRON_RESCATTERING;
using namespace ATOOLS;

Hadron_Rescatterings::Hadron_Rescatterings(const bool & on) :
  m_pref(3./4.), m_on(on), m_nfails(0), 
  p_mm(nullptr), p_bm(nullptr), p_bb(nullptr)
{
  hrpars = new HR_Parameters();
}

void Hadron_Rescatterings::Initialize() {
  p_mm = new MesonMeson();
  p_bm = new BaryonMeson();
  p_bb = new BaryonBaryon();
}

Hadron_Rescatterings::~Hadron_Rescatterings() {
  delete p_mm;
  delete p_bm;
  delete p_bb;
  if (hrpars) { delete hrpars; hrpars = NULL; }
}

bool Hadron_Rescatterings::SelectInteraction(Particle * inpart1,Particle * inpart2) {
  m_flav1 = inpart1->Flav();
  m_flav2 = inpart2->Flav();
  p_interaction = nullptr;
  if      (m_flav1.IsMeson() && m_flav2.IsMeson())    p_interaction = p_mm;
  else if (m_flav1.IsBaryon() && m_flav2.IsBaryon())  p_interaction = p_bb;
  else if ((m_flav1.IsMeson() && m_flav2.IsBaryon()) ||
	   (m_flav1.IsBaryon() && m_flav2.IsMeson())) p_interaction = p_bm;
  return (p_interaction!=nullptr);
}

bool Hadron_Rescatterings::WillRescatter(const double & s,const double & b2) {
  // "Geometric" cross section: sigma_geometric = pi b^2 to be compared with
  // the physical, flavour- and energy-dependent cross section.
  // The extra factor of 10^25 covers the translation from mm^2 to mbarn.
  double xsgeom     = m_pref * M_PI * b2 * 1.e25;
  double xstot      = p_interaction->XStot(m_flav1,m_flav2,s);
  //msg_Out()<<METHOD<<"("<<m_flav1<<", "<<m_flav2<<"): "
  //	   <<"E = "<<sqrt(s)<<", b = "<<sqrt(b2)<<", "
  //	   <<"xs = "<<xsgeom<<" vs. "<<xstot<<" mb "
  //	   <<"--> P = "<<exp(-xsgeom/xstot)<<".\n";
  if (xstot<=1.e-24 || exp(-xsgeom/xstot)<ran->Get()) return false;
  return true;
}

Blob * Hadron_Rescatterings::ElasticRescatter(Hadron_Collision * candidate) {
  //why is b_slope = 2.4 ? https://arxiv.org/pdf/hep-ph/0312187 --> page 5, eq 4a.
  double b_slope  = 2.4;
  double Ecms1    = candidate->MomCMS(0)[0], E2 = sqr(Ecms1);
  double pcms1    = candidate->MomCMS(0)[3], p2 = sqr(pcms1);
  double max_t    = -2.*(E2+p2), min_t = -2.*(E2-p2);
  double expBmax  = exp(b_slope*max_t), expBmin = exp(b_slope*min_t);
  double t        = 1./b_slope * log(ran->Get()*(expBmin-expBmax)+expBmax);
  double costheta = (t+2.*E2)/(2.*p2), sintheta = sqrt(1.-sqr(costheta));
  double phi      = 2.*M_PI*ran->Get();
  Blob * blob     = new Blob();
  blob->SetPosition(candidate->X_Lab());
  blob->SetType(btp::Soft_Collision);
  blob->SetTypeSpec("Rescattering: Sherpa");
  blob->SetStatus(blob_status::needs_hadronRescatter |
		  blob_status::needs_hadrondecays);
  blob->SetId();
  for (size_t i=0;i<2;i++) {
    blob->AddToInParticles(candidate->InPart(i));
    Vec4D mom = ( candidate->InvBoost() *
		  Vec4D(candidate->MomCMS(i)[0],
			candidate->MomCMS(i)[3]*Vec3D(sintheta*cos(phi),sintheta*sin(phi),costheta)) );
    Particle * out = new Particle(-1,candidate->InPart(i)->Flav(),mom,'F');
    out->SetNumber();
    blob->AddToOutParticles(out);
  }
  msg_Out()<<METHOD<<" for "<<min_t<<" < t = "<<t<<" < "<<max_t<<", s = "<<candidate->Ecms2()<<", "
	   <<"p = "<<candidate->MomCMS(0)<<"\n    --> "
	   <<"cos(theta) = "<<(t+2.*E2)<<"/"<<(2.*p2)<<" = "<<costheta<<"\n";
  return blob;
}

Blob * Hadron_Rescatterings::operator()(Hadron_Collision * candidate) {
  if (!SelectInteraction(candidate->InPart(0),candidate->InPart(1)) ||
      !WillRescatter((m_S = candidate->Ecms2()),candidate->B2())) return nullptr;
  return ElasticRescatter(candidate);
}

void Hadron_Rescatterings::Reset() {}


#include "BEAM/Main/Collider_Kinematics.H"
#include "ATOOLS/Org/Run_Parameter.H"
#include "ATOOLS/Org/Exception.H"

using namespace BEAM;
using namespace ATOOLS;
using namespace std;

Collider_Kinematics::Collider_Kinematics(Beam_Base ** beams) :
  Kinematics_Base(beams), m_mode(collidermode::unknown)
{
  if (p_beams[0]->Type()==beamspectrum::monochromatic &&
      p_beams[1]->Type()==beamspectrum::monochromatic)
    m_mode = collidermode::monochromatic;
  else if (p_beams[0]->Type()!=beamspectrum::monochromatic &&
	   p_beams[1]->Type()==beamspectrum::monochromatic)
    m_mode = collidermode::spectral_1;
  else if (p_beams[0]->Type()==beamspectrum::monochromatic &&
	   p_beams[1]->Type()!=beamspectrum::monochromatic)
    m_mode = collidermode::spectral_2;
  else if (p_beams[0]->Type()!=beamspectrum::monochromatic &&
	   p_beams[1]->Type()!=beamspectrum::monochromatic)
    m_mode = collidermode::both_spectral;
  if (m_mode==collidermode::unknown)
    THROW(fatal_error, "Bad settings for collider mode.");
  InitSystem();
  InitIntegration();
}

Collider_Kinematics::~Collider_Kinematics() {}

void Collider_Kinematics::InitSystem() {
  // cms system from beam momenta - this is for potentially asymmetric collisions.
  m_E1   = p_beams[0]->Energy();
  m_E2   = p_beams[1]->Energy();
  m_Ecms = sqrt(m_S);

  rpa->gen.SetEcms(m_Ecms);
  rpa->gen.SetPBeam(0,p_beams[0]->InMomentum());
  rpa->gen.SetPBeam(1,p_beams[1]->InMomentum());
  
  double x  = (m_S+m_m2[0]-m_m2[1])/(2.*m_S);
  double E1 = x*m_Ecms, pz = sqrt(sqr(E1)-m_m2[0]);
  m_fixp_cms[0] = Vec4D(E1,       0.,0., pz);
  m_fixp_cms[1] = Vec4D(m_Ecms-E1,0.,0.,-pz);
  m_asymmetric  =
    ((dabs((m_fixp_cms[0]-p_beams[0]->InMomentum()).Abs2())>0.0000001) ||
     (dabs((m_fixp_cms[1]-p_beams[1]->InMomentum()).Abs2())>0.0000001) );
  m_on = (m_mode != collidermode::monochromatic);
}
  
void Collider_Kinematics::InitIntegration() {
  Beam_Parameters parameters;
  //check for if they have been initialised to other values
  double sminratio = parameters("BEAM_SMIN");
  double smaxratio = parameters("BEAM_SMAX");
  m_xmin = p_beams[0]->Xmin()*p_beams[1]->Xmin();
  m_xmax = p_beams[0]->Xmax()*p_beams[1]->Xmax();
  m_smin = m_S*Max(m_xmin, sminratio);
  m_smax = m_S*Min(m_xmax, smaxratio);
  // the rapidity interval can be done in a better way. ==> to do
  m_ymin = -10.;
  m_ymax = 10.;
  m_exponent[0] = .5;
  m_exponent[1] = .98 * ( p_beams[0]->Exponent() + p_beams[1]->Exponent());
}
  
bool Collider_Kinematics::operator()(ATOOLS::Vec4D * moms) {
  switch (m_mode) {
  case collidermode::monochromatic:
    return MakeMonochromaticBeams(moms);
  case collidermode::spectral_1:
  case collidermode::spectral_2:
  case collidermode::both_spectral:
    return MakeCollinearBeams(moms);
  }
  return false;
}

bool Collider_Kinematics::MakeMonochromaticBeams(ATOOLS::Vec4D * moms) {
  for (size_t beam=0;beam<2;beam++) m_xkey[3*beam+2] = 1.;
  moms[0] = m_fixp_cms[0];
  moms[1] = m_fixp_cms[1];
  return true;
}

bool Collider_Kinematics::MakeCollinearBeams(ATOOLS::Vec4D * moms) {  
  m_sprime = m_sprimekey[3];
  double y = m_ykey[2];
  if ( (m_sprime<m_smin) ||
       //(m_sprime>1.00000001*m_smax) ||
       m_sprimekey[0]==m_sprimekey[1] ) return false;
  double E  = sqrt(m_sprimekey[2]), Eprime = sqrt(m_sprime);
  double x  = (m_sprime+m_m2[0]-m_m2[1])/(2.*m_sprime);
  double E1 = x*Eprime, E2 = Eprime-E1;
  // c.m. momenta
  moms[0] = Vec4D(E1,0.,0.,sqrt(sqr(E1)-m_m2[0]));
  moms[1] = Vec4D(E2,(-1.)*Vec3D(moms[0]));
  // defining the boost
  double coshy = exp(y)+exp(-y), sinhy = exp(y)-exp(-y);
  m_CMSBoost   = Poincare(Vec4D(coshy,0.,0.,sinhy));
  for (size_t i=0;i<2;i++) CalculateAndSetX(i,moms[i]);
  return true;
}

void Collider_Kinematics::CalculateAndSetX(size_t beam,const ATOOLS::Vec4D & p) {
  Vec4D q = p;
  m_CMSBoost.BoostBack(q);
  m_xkey[3*beam+2] = 2.*q[0]/p_beams[beam]->Energy();
  p_beams[beam]->SetX(m_xkey[3*beam+2]);
}

void Collider_Kinematics::AssignKeys(Integration_Info *const info) {
  m_sprimekey.Assign(m_keyid+string("s'"),5,0,info);
  m_ykey.Assign(m_keyid+string("y"),3,0,info);
  m_xkey.Assign(m_keyid+string("x"),6,0,info);
  m_sprimekey[0] = Max(m_smin, m_sminPS);
  m_sprimekey[1] = m_sprimekey[2] = m_smax;
  m_sprimekey[2] = m_S;
  m_sprimekey[3] = m_S;
  m_ykey[0]      = m_ymin;
  m_ykey[1]      = m_ymax;
  m_ykey[2]      = 0.;
}

void Collider_Kinematics::SetLimits() {
  m_sprimekey[0] = Max(m_smin, m_sminPS);
  m_sprimekey[1] = m_sprimekey[2] = m_smax;
  m_sprimekey[3] = m_S;
  m_ykey[0]      = m_ymin;
  m_ykey[1]      = m_ymax;
  m_ykey[2]      = 0.;
  for (size_t i=0;i<2;i++) {
    double p = i==0?
      p_beams[0]->OutMomentum().PPlus():
      p_beams[1]->OutMomentum().PMinus();
    double e  = p_beams[i]->OutMomentum()[0];
    m_xkey[3*i]   = -0.5*(IsZero(m_m[i],1.e-13)?
			  numeric_limits<double>::max():log(p/m_m[i]));
    m_xkey[3*i+1] = log (Min(p_beams[i]->Xmax(),
			     (e/p*(1.0+sqrt(1.0-sqr(m_m[i]/e))))));
  }
  // sprime's with masses - still need to check for masses
  double sprimemin = Max(m_sprimekey[0],m_S*exp(m_xkey[0]+m_xkey[3]));
  if (sprimemin>sqr(m_m[0]+m_m[1])) m_sprimekey[0] = sprimemin;
  double sprimemax = Min(m_smax,m_S*exp(m_xkey[1]+m_xkey[4]));
  if (sprimemax>sqr(m_m[0]+m_m[1])) m_sprimekey[1] = sprimemax;
}


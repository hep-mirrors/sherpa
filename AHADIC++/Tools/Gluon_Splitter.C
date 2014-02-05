#include "AHADIC++/Tools/Gluon_Splitter.H"
#include "AHADIC++/Tools/Hadronisation_Parameters.H"
#include "ATOOLS/Math/Poincare.H"
#include "ATOOLS/Math/Random.H"

using namespace AHADIC;
using namespace ATOOLS;
using namespace std;

Gluon_Splitter::Gluon_Splitter() : 
  Splitter_Base(),
  m_etay(hadpars->Get(string("G2QQ_Exponent"))),
  m_etay_lead(hadpars->Get(string("G2QQ_LeadExponent"))),
  p_out1(NULL), p_out2(NULL)
{
  m_anapath = std::string("glue");
}

bool Gluon_Splitter::
operator()(Dipole * dip,const bool & first,const bool & vetodiquark) {
  Reset();
  if (!SelectSplitter(dip->Triplet(),dip->AntiTriplet())) exit(0);
  DefineTags();
  dip->SetSwitched(!m_swap);
  ConstructTrafos();
  if (ConstructLightC() && ConstructSystem()) {
    if (m_ana) Analysis();
    return true;
  }
  UndoTrafos();
  return false;
}

bool Gluon_Splitter::
SelectSplitter(Proto_Particle * part1,Proto_Particle * part2) {
  Flavour tflav(part1->m_flav), aflav(part2->m_flav);
  bool q1(tflav.IsQuark() || tflav.IsDiQuark());
  bool q2(aflav.IsQuark() || aflav.IsDiQuark());
  bool hit1(part1->m_info=='L' || part1->m_info=='B');
  bool hit2(part2->m_info=='L' || part2->m_info=='B');
  if (q1 && q2) return false;
  // m_swap = true if anti-triplet splits
  m_swap = (q1 && !q2 ||
	    (!q1 && !q2 && 
	     ((((hit1 && hit2) || (!hit1 && !hit2)) && ran->Get()<0.5))));

  p_split = m_swap?part2:part1;
  p_spect = m_swap?part1:part2;
  return true;
}

bool Gluon_Splitter::ConstructSystem() {  
  if (m_LC.m_E-m_LC.m_mspect-sqrt(4.*m_mmin2)<0.) return false;
  double etay(FixExponent()),sqq;
  long int calls(0);
  bool vetodi(false); //m_leadspect
  double pt2max((m_leadspect?1.:1.)*m_pt2max);
  if (m_leadspect) pt2max *= m_pt2max/Max(m_pt2max,m_LC.m_mspect2);
  m_popped.push_back(new PoppedPair);
  do { 
    ConstructKinematics(etay);
  } while (calls++<=100 && 
	   (!SelectFlavour(m_popped.back()->m_sqq,vetodi) || 
	    !AcceptSystem(pt2max)));
  if (calls>100) {
    Reset();
    p_out1 = p_out2 = NULL;
    return false;
  }
  MakeKinematics();
  MakeParticles();
  Reset();
  return true;
}

bool Gluon_Splitter::PoppedMassPossible(const double & m2) {
  PoppedPair * pop(m_popped.back());
  double mhat2(m2/m_LC.m_smandel);
  double delta((mhat2*(1.-pop->m_y))/
	       (pop->m_y*((1.-pop->m_y)*m_LC.m_smandel-m_LC.m_mspect2)));
  if (delta<=0.) return false;
  double zmint(0.5*(1.-sqrt(1.-4.*delta)));
  double zmaxt(0.5*(1.+sqrt(1.-4.*delta)));
  if (pop->m_z<=zmint || pop->m_z>=zmaxt || 
      sqr(m_LC.m_smandel+4.*m2-m_LC.m_mspect2)-16.*m2*m_LC.m_mspect2<=0.) 
    return false;
  double central((m_LC.m_smandel+4.*m2-m_LC.m_mspect2)/(2.*m_LC.m_smandel));
  double lambda(Lambda(m_LC.m_smandel,4.*m2,m_LC.m_mspect2));
  if ((pop->m_y<=central-lambda)||(pop->m_y>=central+lambda)) return false; 
  return true;
}


void Gluon_Splitter::ConstructKinematics(const double & etay) {  
  double mminhat2(m_mmin2/m_LC.m_smandel); 
  double mspecthat2(m_LC.m_mspect2/m_LC.m_smandel);
  double central((m_LC.m_smandel+4.*m_mmin2-m_LC.m_mspect2)/
		 (2.*m_LC.m_smandel));
  double lambda(Lambda(m_LC.m_smandel,4.*m_mmin2,m_LC.m_mspect2));
  double ymin(central-lambda),ymax(central+lambda);
  double offset(m_pt02/m_LC.m_smandel);
  long int calls(0);
  double y,z,sqq,delta,weight;
  do {
    y      = SelectY(ymin,ymax,etay,offset);
    delta  = (mminhat2*(1.-y))/(y*(1.-y-mspecthat2));
    z      = SelectZ(delta,m_leadspect);
    sqq    = y*m_LC.m_smandel*(1.-mspecthat2/(1.-y));
    weight = 1.;//exp(-sqq/(4.*m_pt02));
    calls++;
  } while (weight<ran->Get() && calls<=100);
  if (calls<=100) {
    PoppedPair * pop(m_popped.back());
    pop->m_y = y; pop->m_z = z; pop->m_sqq = sqq;
  }
}

bool Gluon_Splitter::AcceptSystem(const double & pt2max) {
  PoppedPair * pop(m_popped.back());
  pop->m_kt2 = 
    pop->m_y*pop->m_z*(1.-pop->m_z)*
    (m_LC.m_smandel-m_LC.m_mspect2/(1.-pop->m_y))-pop->m_mpop2;
  if (pop->m_kt2 < 0.)     return false;
  if (pop->m_kt2 > pt2max) return false;
  return exp(-(pop->m_sqq-4.*pop->m_mpop2)/(4.*m_pt02)) > ran->Get();
  //  return (*p_as)(pop->m_sqq,false)/p_as->MaxValue() > ran->Get();
}

double Gluon_Splitter::FixExponent() {
  if (m_isbeam) return m_etay_lead+1.;
  if (m_leadspect && p_spect->m_flav.IsQuark()) return m_etay_lead;
  return m_etay;
}

void Gluon_Splitter::MakeKinematics() {
  PoppedPair * pop(m_popped.back());
  double kt(sqrt(pop->m_kt2));
  double phi(2.*M_PI*ran->Get());
  Vec4D kperp(0.,kt*cos(phi),kt*sin(phi),0.);
  double beta(m_LC.m_mspect2/(m_LC.m_smandel*(1.-pop->m_y)));
  pop->m_outmom[0]  = 
    pop->m_y*pop->m_z*m_LC.m_pB + 
    (1.-beta)*(1.-pop->m_z)*m_LC.m_pA + kperp;
  pop->m_outmom[1]  = 
    pop->m_y*(1.-pop->m_z)*m_LC.m_pB + 
    (1.-beta)*pop->m_z*m_LC.m_pA - kperp;
  m_spectmom = (1.-pop->m_y)*m_LC.m_pB + beta*m_LC.m_pA;
  for (size_t i(0);i<2;i++) {
    m_rotat.RotateBack(pop->m_outmom[i]);
    m_boost.BoostBack(pop->m_outmom[i]);
  }
  m_rotat.RotateBack(m_spectmom);
  m_boost.BoostBack(m_spectmom);
}

void Gluon_Splitter::MakeParticles() {
  PoppedPair * pop(m_popped.back());
  double m1((m_spectmom+pop->m_outmom[0]).Abs2());
  double m2((m_spectmom+pop->m_outmom[1]).Abs2());
  bool swap(/*ran->Get()>0.5);//*/m_leadspect?m1>m2:m2<m1);
  char info(p_split->m_info=='B'||p_spect->m_info=='B'?'B':'l');
  p_out1 = new Proto_Particle(pop->m_flav.Bar(),pop->m_outmom[int(swap)],info);
  p_out2 = new Proto_Particle(pop->m_flav,pop->m_outmom[1-int(swap)],info);
  p_out1->p_partner = p_out2;
  p_out2->p_partner = p_out1;
  p_out1->m_kt2max  = p_out2->m_kt2max = pop->m_kt2;
  p_spect->m_mom    = m_spectmom;
  delete p_split;
}



#include "AHADIC++/Formation/Colour_Reconnections.H"
#include "ATOOLS/Org/Message.H"
#include "ATOOLS/Math/Random.H"

using namespace AHADIC;
using namespace ATOOLS;
using namespace std;


Colour_Reconnections::Colour_Reconnections(int kinmode,int colmode,double t0) : 
  m_kinmode(kinmode), m_colmode(colmode), m_t0(t0)
{ }

Colour_Reconnections::~Colour_Reconnections() 
{ }

void Colour_Reconnections::Singlet_CR(Cluster_List * clin)
{
  if (clin->size()<2) return;
  Cluster_Iterator cit1,cit2;
  int gen=1;
  bool direction = (ran.Get()>0.5);
  if (direction) clin->reverse();

  int lead1,lead2;
  cit1 = cit2 = clin->begin(); cit2++;
  do {
    if (TestClusters((*cit1),(*cit2),gen)) {
      Proto_Particle * help = (*cit1)->GetAnti();

      (*cit1)->SetAnti((*cit2)->GetAnti());
      (*cit2)->SetAnti(help);

      (*cit1)->Update();
      (*cit2)->Update();
      if (m_w14/(m_w23+m_w14)>ran.Get()) { cit2++; gen++; }
                                    else { cit1 = cit2; cit2++; gen=1; }
    } 
    else {
      if (m_w12/(m_w12+m_w34)>ran.Get()) { cit2++; gen++; }
                                    else { cit1 = cit2; cit2++; gen=1; }
    }
  } while (cit2!=clin->end());
}

void Colour_Reconnections::Two_Singlet_CR(Cluster_List * cl1,Cluster_List * cl2)
{
  // To be filled.
}

bool Colour_Reconnections::TestClusters(Cluster * cl1,Cluster * cl2,int gen)
{
  double kinweight = KinematicWeight(cl1->GetTrip()->m_mom,cl1->GetAnti()->m_mom,
				     cl2->GetTrip()->m_mom,cl2->GetAnti()->m_mom);
  double colweight = ColourWeight(gen);
  if (kinweight*colweight>ran.Get()) return true;
  return false;
}

double Colour_Reconnections::KinematicWeight(const Vec4D & mom1,const Vec4D & mom2,
					     const Vec4D & mom3,const Vec4D & mom4)
{
  switch (m_kinmode) {
  case 2:
  case 1:
    m_w12 = sqrt((mom1+mom2).Abs2());
    m_w34 = sqrt((mom3+mom4).Abs2());
    m_w14 = sqrt((mom1+mom4).Abs2());
    m_w23 = sqrt((mom2+mom3).Abs2());
    break;
  default:
    m_w12 = m_w34 = 0;
    m_w14 = m_w23 = 1.e64;
    break;
  }
  double w1234 = m_t0/(m_t0+4.*sqr(m_w12+m_w34));
  double w1423 = m_t0/(m_t0+4.*sqr(m_w14+m_w23));

  return w1423/(w1234+w1423);
}


double Colour_Reconnections::ColourWeight(int gen) { 
  switch (m_colmode) {
  case 1:
    if (gen==0) return 0.;
    if (gen==1) return 1./9.;
    return pow(1./9.,gen); 
  default: return 0.;
  }
}

#include "Cluster_Formation_Handler.H"

using namespace AHADIC;
using namespace ATOOLS;

Cluster_Formation_Handler::Cluster_Formation_Handler(Clusters_2_Hadrons * _transformer) :
  p_newmasses(NULL), p_oldens2(NULL), p_ens(NULL), m_xmt(0.), m_number(0),
  p_transformer(_transformer)
{
  p_partlist     = new Part_List;
  p_gludecayer   = new Gluon_Decayer(); 
  p_cformer      = new Cluster_Former();
}


Cluster_Formation_Handler::~Cluster_Formation_Handler()
{
  p_partlist->clear(); delete p_partlist;

}

Blob *  Cluster_Formation_Handler::FormClusters(Blob_List * bl) 
{
  if (bl==NULL) return false;
  p_partlist->clear();
  Blob * blob = new Blob();
  blob->SetType(btp::Cluster_Formation);
  blob->SetId();
  bl->push_back(blob);
  Shift(bl,blob);
  p_gludecayer->DecayList(p_partlist);
  p_cformer->FormClusters(p_partlist);
  return blob;
}

void Cluster_Formation_Handler::Shift(Blob_List * bl,Blob * blob)
{
  Particle * part1, * part2;
  for (Blob_Iterator blit=bl->begin();blit!=bl->end();++blit) {
    if ((*blit)->Type()==btp::FS_Shower || (*blit)->Type()==btp::IS_Shower) {
      for (int i=0;i<(*blit)->NOutP();i++) {
	part2 = (*blit)->OutParticle(i); 
	if (part2->Status()==1) {
	  blob->AddToInParticles(part2);
	  part1 = new Particle(-1,part2->Flav(),part2->Momentum(),'L');
	  part1->SetNumber(0);
	  part1->SetFlow(1,part2->GetFlow(1));
	  part1->SetFlow(2,part2->GetFlow(2));
	  p_partlist->push_back(part1);
	}
      }
    }
  }
  if (p_partlist->size()<2) {
    return; 
  }
  if (p_partlist->size()==2) {
    part1 = p_partlist->front(); part2 = p_partlist->back();
    Flavour fl1     = part1->Flav(), fl2 = part2->Flav();
    double m12      = sqr(hadpars.Mass(fl1));
    double m22      = sqr(hadpars.Mass(fl2));
    double energy   = part1->Momentum()[0] + part2->Momentum()[0];
    double energy0  = (sqr(energy)+m12-m22)/(2.*energy);
    double energy1  = (sqr(energy)-m12+m22)/(2.*energy);
    Vec3D direction = Vec3D(part1->Momentum())/(Vec3D(part1->Momentum()).Abs());
    Vec3D p0        = direction*sqrt(sqr(energy0)-m12);
    Vec3D p1        = (-1.)*direction*sqrt(sqr(energy1)-m12);
    part1->SetMomentum(Vec4D(energy0,p0));
    part2->SetMomentum(Vec4D(energy1,p1));
    return; 
  }

  m_xmt          = 0.;
  m_number       = p_partlist->size();
  p_newmasses    = new double[m_number];
  p_oldens2      = new double[m_number];
  p_ens          = new double[m_number];
  Vec4D cms      = Vec4D(0.,0.,0.,0.);
  short int k    = 0;
  Flavour flav;
  for (Part_Iterator pit=p_partlist->begin();pit!=p_partlist->end();++pit,++k) {
    flav         = (*pit)->Flav();
    m_xmt       += p_newmasses[k] = hadpars.Mass(flav);
    cms         += (*pit)->Momentum();
    p_oldens2[k] = sqr((*pit)->Momentum()[0]);
  }
  double ET  = sqrt(cms.Abs2()); 
  double x   = sqrt(1.-sqr(m_xmt/ET));
  double acc = ET*1.e-14;

  double f0,g0,x2;
  for (k=0;k<10;k++) {
    f0 = -ET;g0 = 0.;x2 = x*x;
    for (short int i=0;i<m_number;i++) {
      p_ens[i] = sqrt(sqr(p_newmasses[i])+x2*p_oldens2[i]);
      f0      += p_ens[i];
      g0      += p_oldens2[i]/p_ens[i];
    }
    if (dabs(f0)<acc) break; 
    x -= f0/(x*g0);  
  }
  // Construct Momenta
  k = 0;
  for (Part_Iterator pit=p_partlist->begin();pit!=p_partlist->end();++pit,++k) {
    (*pit)->SetMomentum(Vec4D(p_ens[k],x*Vec3D((*pit)->Momentum())));
  }

  delete [] p_oldens2;
  delete [] p_ens;
  delete [] p_newmasses;
}

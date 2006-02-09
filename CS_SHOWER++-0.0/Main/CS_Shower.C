#include "CS_Shower.H"
#include "ISR_Handler.H"
#include "Run_Parameter.H"

using namespace CS_SHOWER;
using namespace ATOOLS;
using namespace std;

CS_Shower::CS_Shower(PDF::ISR_Handler * _isr, bool _isron, bool _fsron) : 
  m_isron(_isron), m_fsron(_fsron), p_allsinglets(new All_Singlets), m_shower(_isr), 
  p_momenta(new Vec4D[4]), p_flavours(new Flavour[4])
{
}

CS_Shower::~CS_Shower() 
{
  if (p_allsinglets) { delete p_allsinglets; p_allsinglets = NULL; }
  
  delete [] p_momenta;
  delete [] p_flavours;
}



bool CS_Shower::PerformShowers(All_Singlets * allsinglets) {
  if (allsinglets!=NULL) {
    if (p_allsinglets->size()>0) {
      for (ASiter asit=p_allsinglets->begin();asit!=p_allsinglets->end();asit++) {
	Singlet * sing = (*asit);
	if (sing) { delete sing; sing = NULL; }
      }
      p_allsinglets->clear();
    }
    p_allsinglets = allsinglets;
  }
  return m_shower.EvolveShower(p_allsinglets);
}

bool CS_Shower::ExtractPartons(ATOOLS::Blob_List * bloblist) {
 
  Blob * fs=NULL, * is=NULL, * hard=NULL, * ise=NULL;

  for (deque<Blob*>::iterator blit=bloblist->begin();blit!=bloblist->end();blit++) {
    if ((*blit)->Type()==btp::Signal_Process)       hard = *blit;
    if ((*blit)->Type()==btp::Shower)               fs = *blit;
    if ((*blit)->Type()==btp::IS_Shower)
      if ((*blit)->OutParticle(0)->Flav().Strong()) is = *blit;
      else ise = *blit;
  }
  
  if (fs==NULL || hard==NULL) return false;
  if ((is==NULL || ise==NULL) && m_isron)    return false;
  
  ExtractPartons(fs,is);

  if (is && fs) {
    if (InitBreitFrame(hard)) { 
      //boost momenta from Breit Frame back to the LAB 
      Vec4D tmp;
      tmp = is->InParticle(0)->Momentum();
      BoostBack(tmp);
      is->InParticle(0)->SetMomentum(tmp);
      for (int i=0;i<fs->NOutP();i++) {
	tmp = fs->OutParticle(i)->Momentum();
	BoostBack(tmp);
	fs->OutParticle(0)->SetMomentum(tmp);
      }
      for(int i=0; i<hard->NInP(); ++i) {
	if(hard->InParticle(i)->Flav().Strong()==false) {
	  //save
	  Particle* par=new Particle(*hard->InParticle(i));
	  ise->AddToInParticles(par);
	}
      }
      /*
      for(int i=0; i<hard->NOutP(); ++i) {
	if(hard->OutParticle(i)->Flav().Strong()==false) {
	Particle* par=new Particle(*hard->OutParticle(i));
	fs->AddToOutParticles(par);
	}
      }
      */
    }
    else {
      std::cout<<" ERROR : Could not initialize the Breit Frame !!! "<<std::endl; 
    }
  }
  return true;
}

void CS_Shower::ExtractPartons(ATOOLS::Blob * fs,ATOOLS::Blob * is) {
  for (ASiter asit=p_allsinglets->begin();asit!=p_allsinglets->end();asit++) {
    (*asit)->ExtractPartons(fs,is);
  }
}

void CS_Shower::PrepareAllSinglets()
{
  for (ASiter asit=p_allsinglets->begin();asit!=p_allsinglets->end();asit++) {
    Singlet * sing = (*asit);
    if (sing) { delete sing; sing = NULL; }
  }
  p_allsinglets->clear();
}


bool CS_Shower::InitBreitFrame(Blob *blob) 
{
  if (blob==NULL) return false;
  if ((blob->NInP()!=2) || (blob->NOutP()!=2)) {
     if ((blob->NInP()!=2) || (blob->NOutP()!=2)) {
       THROW(fatal_error,"Cannot handle blobs with more than 4 legs.");
       return false;
     }
  }
  
  for (size_t i=0;i<(size_t)blob->NInP();++i) {
    p_flavours[i]=blob->InParticle(i)->Flav();
    p_momenta[i]=blob->InParticle(i)->Momentum();
  }
  
  for (size_t i=0;i<(size_t)blob->NOutP();++i) {
    p_flavours[i+blob->NInP()]=blob->OutParticle(i)->Flav();
    p_momenta[i+blob->NInP()]=blob->OutParticle(i)->Momentum();
  }
  
  int lepton=0;
  if (p_flavours[0].Strong()) lepton=1;
  
  Vec4D q = p_momenta[lepton];
  
  for (size_t i=0;i<(size_t)blob->NOutP();++i) {
    if (p_flavours[i+blob->NInP()]==p_flavours[lepton]) {
      q-=p_momenta[i+blob->NInP()];
    }
  }
    
  int hadron = rpa.gen.Beam1().Strong()==1?0:1;
  Vec4D beam = rpa.gen.PBeam(hadron);
  Vec4D store = q;
  /*
  double Eprime = ATOOLS::rpa.gen.Ecms();
  double m02 = rpa.gen.PBeam(0)*rpa.gen.PBeam(0);
  double m12 = rpa.gen.PBeam(1)*rpa.gen.PBeam(1);
  double x1 = 1./2.+(m02-m12)/(2.*sqr(Eprime));
  double E1 = x1*Eprime;
  double E2 = (1.-x1)*Eprime;
  
  Vec4D * beam = new Vec4D[2];
  
  beam[0] = Vec4D(E1,0.,0.,sqrt(sqr(E1)-m02));
  beam[1] = Vec4D(E2,(-1.)*Vec3D(beam[0]));
  
  std::cout<<" beams "<<beam[0]<<std::endl;
  std::cout<<" beams "<<beam[1]<<std::endl;
  */

  double x = -q.Abs2()/(2.*beam*q); 
  Vec4D pp = 2.*x*beam+q;
  
  double gamma = pp[0]/pp.Abs();
  Vec3D eta    = Vec3D(pp)/pp.Abs();

  m_boost = Poincare(Vec4D(gamma,eta));
  m_boost.Boost(q);
  m_zrot  = Poincare(-1.*q,Vec4D::ZVEC);
  m_zrot.Rotate(q);
  
  BoostBack(q);

  //checks
  if (dabs(q*pp)>1.e-9) {
    msg.Error()<<" ERROR: CS_Shower::InitBreitFrame could not initialize Breit frame correctly (1) : "
	       <<dabs(q*pp)<<std::endl;
    return false;
  }
  for (int i=0;i<3;i++) {
    if (dabs((q[i]-store[i]))>1.e10) {
      msg.Error()<<" ERROR: CS_Shower::InitBreitFrame could not initialize Breit frame correctly (2) : "
		 <<q-store<<std::endl;
      return false;
    }
  }
  return true;
}

void CS_Shower::BoostBack(Vec4D & p)
{
  m_zrot.RotateBack(p);
  m_boost.BoostBack(p);
}

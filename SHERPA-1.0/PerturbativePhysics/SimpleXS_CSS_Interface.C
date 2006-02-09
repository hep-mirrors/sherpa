#include "SimpleXS_CSS_Interface.H"

#include "CS_Shower.H"
#include "Run_Parameter.H"
#include "XS_Base.H"
#include "Exception.H"

using namespace SHERPA;
using namespace CS_SHOWER;
using namespace EXTRAXS;
using namespace ATOOLS;

SimpleXS_CSS_Interface::SimpleXS_CSS_Interface(Matrix_Element_Handler *_p_mehandler,
					       Shower_Handler *_p_shower) :
  Perturbative_Interface(_p_mehandler,_p_shower), p_momenta(new Vec4D[4]), 
  p_flavours(new Flavour[4]), p_blob(NULL), p_hard(NULL)
{ 
  m_hadron = rpa.gen.Beam1().Strong()==1?0:1;
  m_lepton = m_hadron==1?0:1;
}

SimpleXS_CSS_Interface::~SimpleXS_CSS_Interface() 
{ 
  delete [] p_momenta;
  delete [] p_flavours;
}

int SimpleXS_CSS_Interface::DefineInitialConditions(Blob * blob) 
{
  if (blob==NULL) return false;
  if ((blob->NInP()!=2) || (blob->NOutP()!=2)) {
    THROW(fatal_error,"Cannot handle blobs with more than 4 legs.");
  }
  p_hard = blob;
  p_shower->CleanUp();
  return InitColours(blob);
}

bool SimpleXS_CSS_Interface::InitBreitFrame(Blob *blob) 
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
  
  Vec4D store = q;
  
  double Eprime = ATOOLS::rpa.gen.Ecms();
  double m02 = rpa.gen.PBeam(0)*rpa.gen.PBeam(0);
  double m12 = rpa.gen.PBeam(1)*rpa.gen.PBeam(1);
  double x1 = 1./2.+(m02-m12)/(2.*sqr(Eprime));
  double E1 = x1*Eprime;
  double E2 = (1.-x1)*Eprime;
  
  Vec4D * beam = new Vec4D[2];
  
  beam[0] = Vec4D(E1,0.,0.,sqrt(sqr(E1)-m02));
  beam[1] = Vec4D(E2,(-1.)*Vec3D(beam[0]));
  
  m_x = -q.Abs2()/(2.*beam[m_hadron]*q); 
  Vec4D pp = 2.*m_x*beam[m_hadron]+q;

  double gamma = pp[0]/pp.Abs();
  Vec3D eta    = Vec3D(pp)/pp.Abs();

  m_boost = Poincare(Vec4D(gamma,eta));
  m_boost.Boost(q);
  m_zrot  = Poincare(-1.*q,Vec4D::ZVEC);
  m_zrot.Rotate(q);
  
  m_zrot.RotateBack(q);
  m_boost.BoostBack(q);
  
  //checks
  if (dabs(q*pp)>1.e-9) {
    msg.Error()<<" ERROR: SimpleXS_CSS_Interface::InitBreitFrame could not initialize Breit frame correctly (1) : "
	       <<dabs(q*pp)<<std::endl;
    return false;
  }
  for (int i=0;i<3;i++) {
    if (dabs((q[i]-store[i]))>1.e10) {
      msg.Error()<<" ERROR: SimpleXS_CSS_Interface::InitBreitFrame could not initialize Breit frame correctly (2) : "
		 <<q-store<<std::endl;
      return false;
    }
  }
  delete [] beam;
  return true;
}

void SimpleXS_CSS_Interface::BoostInBreitFrame(Vec4D & p)
{
  m_boost.Boost(p);
  m_zrot.Rotate(p);
}

bool SimpleXS_CSS_Interface::InitColours(Blob * blob) 
{
  XS_Base * xs = p_mehandler->GetXS();
  if (!(xs->SetColours(p_mehandler->Momenta()))) return false;
  for (int j=0;j<2;j++) {
    for (int i=0;i<blob->NInP();++i) {
      blob->InParticle(i)->SetFlow(j+1,xs->Colours()[i][j]);
    }
    for (int i=0;i<blob->NOutP();++i) blob->OutParticle(i)->SetFlow(j+1,xs->Colours()[i+blob->NInP()][j]);
  }

  std::cout<<" InitColours -> "<<*blob<<std::endl;

  if(Valid(blob)) {
    p_blob=ConstructSinglets(blob,xs->Scale(PHASIC::stp::as));
    return true;
  }
  return false;
}

bool SimpleXS_CSS_Interface::Valid(Blob* blob) {
  int strong=0,aflag=0,part=0;
  unsigned int flow=99;
  
  for(int i=0; i<blob->NInP(); ++i) {
    for(int j=1; j<3; j++) {
      if(blob->InParticle(i)->GetFlow(j)!=0) { 
	flow=blob->InParticle(i)->GetFlow(j);
	part=i;aflag=j;
      }
    }
    if(blob->InParticle(i)->Flav().Strong()) strong++;
  }
  
  std::cout<<"VALID : "<<strong<<std::endl;
  if (strong==0)
    if(blob->OutParticle(0)->Flav().IsAnti()) blob->SwapOutParticles(0,1);
  
  if (strong==0) {
    if(!(blob->OutParticle(0)->GetFlow(1)!=0 &&
	 blob->OutParticle(1)->GetFlow(2)!=0)) return false;
  }
  if (strong==1) 
    if(!(flow==blob->OutParticle(part)->GetFlow(aflag))) return false;
 
  if (strong==2) return true;

  return true;
}


Blob* SimpleXS_CSS_Interface::ConstructSinglets(Blob * meblob, const double scale)
{
  Blob * psblob = new Blob();
  psblob->SetStatus(1);
  psblob->SetType(btp::Shower);
  psblob->SetTypeSpec("CSS++0.0");

  std::map<Particle*,pst::code> particles;
  std::map<Particle*,pst::code>::iterator piter,piter1;
  for (int i=0;i<2;i++) {
    if (meblob->InParticle(i)->Flav().Strong())  particles[meblob->InParticle(i)]=pst::IS;
    if (meblob->OutParticle(i)->Flav().Strong()) particles[meblob->OutParticle(i)]=pst::FS;
  }
  int psize(particles.size());
  if (psize==0) {
    msg.Error()<<"ERROR in SimpleXS_CSS_Interface::ConstructSinglets: "<<std::endl
	       <<"   No coloured particles in blob : "<<std::endl<<(*meblob)<<"   abort."<<std::endl;
    abort();
  }
  bool notcomplete(true), initsinglet(false);
  int flowrun(0), flowend(0);
  Singlet * singlet(NULL);
  while (psize>0) {
    initsinglet = false;
    for (piter=particles.begin();piter!=particles.end();piter++) {
      // triplet FS start
      if (piter->first->GetFlow(1)>0 && piter->first->GetFlow(2)==0 &&
	  piter->second==pst::FS) {
	flowrun = piter->first->GetFlow(1);
	flowend = 0;
	initsinglet = true;
	break;
      }
      // antitriplet IS start
      else if (piter->first->GetFlow(1)==0 && piter->first->GetFlow(2)>0 &&
	  piter->second==pst::IS) {
	flowrun = piter->first->GetFlow(2);
	flowend = 0;
	initsinglet = true;
	break;
      }
    }
    // ring made of glue
    if (!initsinglet) {
      for (piter=particles.begin();piter!=particles.end();piter++) {
	if (piter->first->GetFlow(1)>0 && piter->first->GetFlow(2)>0) {
	  flowrun = piter->first->GetFlow(1);
	  flowend = piter->first->GetFlow(2);
	  // This is a tricky point, need to close the coloured ring ....
	  initsinglet = true;
	  break;
	}
      }
    }
    if (initsinglet) {
      singlet = new Singlet(flowrun,true);
      p_shower->GetAllSinglets()->push_back(singlet);
      psblob->AddToInParticles(piter->first);
      Parton  * parton  = new Parton(piter->first,piter->second);
      std::cout<<" Calc ? "<<piter->first<<" "<<int(piter->second)<<std::endl;
      if (piter->second==pst::IS) {
	std::cout<<"Calc xbj"<<std::endl;
	if (Vec3D(piter->first->Momentum())*Vec3D(rpa.gen.PBeam(0))>0.)
	  parton->SetXbj(piter->first->Momentum()[0]/rpa.gen.PBeam(0)[0]);
	else 
	  parton->SetXbj(piter->first->Momentum()[0]/rpa.gen.PBeam(1)[0]);
      }
      std::cout<<"Check xbj = "<<parton->Xbj()<<std::endl;
      parton->SetStart(scale/4.);
      parton->SetVeto(scale/4.);
      singlet->push_back(parton);
      std::cout<<"Init : "<<(*singlet)<<std::endl;
      piter1  = piter;
      particles.erase(piter);
      notcomplete = true;
      std::cout<<"Before do, psize = "<<particles.size()<<std::endl;
      do {
	for (piter1=particles.begin();piter1!=particles.end();piter1++) {
	  if (piter1->first->GetFlow(1)==flowrun) {
	    flowrun = piter1->first->GetFlow(2);
	    if (flowrun==flowend) notcomplete=false;
	    psblob->AddToInParticles(piter1->first);
	    Parton  * parton  = new Parton(piter1->first,piter1->second);
	    std::cout<<" Calc ? "<<piter1->first<<" "<<int(piter1->second)<<std::endl;
	    if (piter1->second==pst::IS) {
	      std::cout<<"Calc xbj"<<std::endl;
	      if (Vec3D(piter1->first->Momentum())*Vec3D(rpa.gen.PBeam(0))>0.)
		parton->SetXbj(piter1->first->Momentum()[0]/rpa.gen.PBeam(0)[0]);
	      else 
		parton->SetXbj(piter1->first->Momentum()[0]/rpa.gen.PBeam(1)[0]);
	    }
	    std::cout<<"Check xbj = "<<parton->Xbj()<<std::endl;
	    parton->SetStart(scale/4.);
	    parton->SetVeto(scale/4.);
	    singlet->push_back(parton);
	    particles.erase(piter1);
	    std::cout<<"Iter : "<<(*singlet)<<std::endl;
	    break;
	  }
	  else if (piter1->first->GetFlow(2)==flowrun) {
	    flowrun = piter1->first->GetFlow(1);
	    if (flowrun==flowend) notcomplete=false;
	    psblob->AddToInParticles(piter1->first);
	    Parton  * parton  = new Parton(piter1->first,piter1->second);
	    std::cout<<" Calc ? "<<piter1->first<<" "<<int(piter1->second)<<std::endl;
	    if (piter1->second==pst::IS) {
	      std::cout<<"Calc xbj"<<std::endl;
	      if (Vec3D(piter1->first->Momentum())*Vec3D(rpa.gen.PBeam(0))>0.)
		parton->SetXbj(piter1->first->Momentum()[0]/rpa.gen.PBeam(0)[0]);
	      else 
		parton->SetXbj(piter1->first->Momentum()[0]/rpa.gen.PBeam(1)[0]);
	    }
	    std::cout<<"Check xbj = "<<parton->Xbj()<<std::endl;
	    parton->SetStart(scale/4.);
	    parton->SetVeto(scale/4.);
	    singlet->push_back(parton);
	    particles.erase(piter1);
	    std::cout<<"Iter : "<<(*singlet)<<std::endl;
	    break;
	  }
	}
      } while (notcomplete);
      psize = particles.size();
      if (psize==0) break;
    }
  }
  return psblob;
}


bool SimpleXS_CSS_Interface::FillBlobs(Blob_List * blobs)
{
  p_blob->SetId();
  blobs->push_back(p_blob);

  /*
    if (p_shower->ISROn()==1 && p_shower->FSROn()==1) {
    Blob * p_is = new Blob();
    p_is->SetStatus(1);
    p_is->SetType(btp::IS_Shower);
    p_is->SetTypeSpec("CSS++0.0");
    p_is->SetId();
    p_is->SetBeam(1);
    p_is->AddToOutParticles(p_hard->InParticle(1));
    blobs->push_front(p_is);
    
    p_is = new Blob();
    p_is->SetStatus(1);
    p_is->SetType(btp::IS_Shower);
    p_is->SetTypeSpec("CSS++0.0");
    p_is->SetId();
    p_is->SetBeam(0);
    p_is->AddToOutParticles(p_hard->InParticle(0));
    blobs->push_front(p_is);
    }
  */
  return true;
}

int SimpleXS_CSS_Interface::PerformShowers()
{
  std::cout<<"Perform showers."<<std::endl;
  return p_shower->PerformShowers(false,1,1.,1.,rpa.gen.Ycut());
}










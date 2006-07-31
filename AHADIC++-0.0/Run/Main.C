#include "Ahadic.H"
#include "Apacic.H"
#include "Tree.H"
#include "Knot.H"
#include "Model_Base.H"
#include "Standard_Model.H"
#include "ISR_Base.H"
#include "ISR_Handler.H"
#include "Intact.H"
#include "Blob_List.H"
#include "Particle_List.H"
#include "Run_Parameter.H"
#include "Data_Read.H"
#include "Random.H"


using namespace AHADIC;
using namespace ATOOLS;
using namespace std;


int main(int argc,char* argv[]) 
{  
  rpa.Init(string("./"),string("Run.dat"),argc,argv);
  ParticleInit(string("./"));
  Ahadic * ahadic = new Ahadic(string("./"),string("Cluster.dat"),true);

  MODEL::Model_Base * model     = new MODEL::Standard_Model(string("./"),
							    string("Model.dat"));
  PDF::ISR_Base    ** isrbases  = new PDF::ISR_Base*[2];
  isrbases[0]                   = new PDF::Intact(Flavour(kf::e));     
  isrbases[1]                   = new PDF::Intact(Flavour(kf::e).Bar());     
  PDF::ISR_Handler  * isr       = new PDF::ISR_Handler(isrbases);

  Data_Read         * dataread  = new Data_Read(string("./")+string("Shower.dat"));
  Jet_Finder jf(rpa.gen.Ycut(),1,false);
  APACIC::Apacic    * apacic    = new APACIC::Apacic(isr,model,&jf,dataread);
  delete dataread;

  APACIC::Tree      * tree      = (*apacic->FinTrees()->begin());
  APACIC::Knot      * mo;

  Blob              * blob;
  Blob_List         * blobs     = new Blob_List;
  Particle_List     * particles = new Particle_List;

  Flavour mo_flavs[2];
  rpa.gen.SetEcms(91.2);
  double p, E = rpa.gen.Ecms()/2., E2 = 4.*E*E, costh,phi;
  Vec3D pvec;

  msg.Out()<<"---------------------- Generate "<<rpa.gen.NumberOfEvents()
	   <<" events ("<<2.*E<<" GeV) -----------------------."<<endl;
  for (int i=1;i<=rpa.gen.NumberOfEvents();i++) {
    if (i%1000==0) msg.Out()<<" "<<i<<" th event "<<std::endl;
    blobs->Clear();
    particles->clear();

    Blob::Reset();
    Particle::Reset();
    Flow::ResetCounter();

    apacic->PrepareTrees();
    tree->SetVetoScale(E2);

    mo          = tree->NewKnot();
    mo->t       = E2;
    mo->E2      = E2;
    mo->maxpt2  = 0.;
    mo->z       = 0.5;
    mo->costh   = -1.; 
    mo->thcrit  = 1.e6;
    mo->phi     = 2.*M_PI*ran.Get();
    mo->stat    = 1;  
    *(mo->part) = Particle(1,Flavour(kf::photon),Vec4D(2.*E,0.,0.,0.));
    mo->part->SetStatus(part_status::decayed);
    mo->part->SetInfo('M');
    mo->didkin  = true;

    //mo_flavs[0] = Flavour(kf::code(1+int(ran.Get()*3.)));   
    mo_flavs[0] = Flavour(kf::b);   
    mo_flavs[1] = mo_flavs[0].Bar();
    p           = sqrt(E*E-sqr(mo_flavs[0].PSMass()));
    costh       = 1.-2.*ran.Get();
    phi         = 2.*M_PI*ran.Get();
    pvec        = p*Vec3D(sqrt(1.-sqr(costh))*sin(phi),sqrt(1.-sqr(costh))*cos(phi),costh);

    mo->left             = tree->NewKnot();
    mo->left->prev       = mo;
    mo->left->stat       = 3;    
    mo->left->t          = mo->t;
    mo->left->tout       = sqr(mo_flavs[0].PSMass());
    mo->left->maxpt2     = 0.;
    mo->left->E2         = E*E;
    mo->left->thcrit     = mo->thcrit;
    *(mo->left->part)    = Particle(2,mo_flavs[0],Vec4D(E,pvec));
    mo->left->part->SetStatus(part_status::active);
    mo->left->part->SetInfo('H');
    mo->left->part->SetFlow(1,-1);
    mo->left->didkin     = true;
 
    mo->right            = tree->NewKnot();
    mo->right->prev      = mo;
    mo->right->stat      = 3;     
    mo->right->t         = mo->t;
    mo->right->tout      = sqr(mo_flavs[1].PSMass());
    mo->right->maxpt2    = 0.;
    mo->right->E2        = E*E;
    mo->right->thcrit    = mo->thcrit;
    *(mo->right->part)   = Particle(3,mo_flavs[1],Vec4D(E,-1.*pvec)); 
    mo->right->part->SetStatus(part_status::active);
    mo->right->part->SetInfo('H');
    mo->right->part->SetFlow(2,mo->left->part->GetFlow(1));
    mo->right->didkin    = true;

    blob        = new Blob();
    blob->SetType(btp::ME_PS_Interface_FS);
    blob->SetId();
    blobs->push_back(blob);
    blob->AddToInParticles(new Particle((*mo->part)));

    if (apacic->PerformShowers(false,true,1.,1.)) {
      apacic->ExtractPartons(false,true,blobs,particles);
      ahadic->Hadronize(blobs);
    }
  }

  delete ahadic;
  delete apacic;
  delete isr;
  delete model;
}


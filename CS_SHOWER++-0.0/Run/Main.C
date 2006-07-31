#include "CS_Shower.H"
#include "Singlet.H"
#include "Parton.H"
#include "Random.H"
#include "Particle.H"
#include "Blob.H"
#include "Run_Parameter.H"
#include "Standard_Model.H"
#include "Model_Base.H"
#include "Histogram.H"

using namespace CS_SHOWER;
using namespace ATOOLS;
using namespace std;


int main(int argc,char* argv[]) 
{  
  ParticleInit(string("./"));
  rpa.Init(string("./"),string("Run.dat"),argc,argv);
  MODEL::Model_Base * model = new MODEL::Standard_Model(string("./"),
							string("Model.dat"));
  rpa.gen.SetModel(model);
  Particle * part;
  Blob     * blob(NULL);
  Parton   * parton;
  Singlet  * singlet(NULL);
  All_Singlets * allsinglets = new All_Singlets;
  CS_Shower shower;
  

  ATOOLS::Histogram multi(0,2.,12.,10);
  for (int i=0;i<1000000;i++) {
    if (i%1000==0) cout<<i<<endl; 
    if (!singlet) singlet = new Singlet(true);
    if (!blob)    blob    = new Blob();
    Particle::ResetCounter();
    Particle::Reset();
    part   = new Particle(0,Flavour(kf::u),Vec4D(45.6,0.,0.,45.6),'F');
    part->SetNumber(0);
    part->SetStatus(part_status::decayed);
    blob->AddToInParticles(part);
    parton = new Parton(part,pst::FS);
    parton->SetStart(8317.44);
    parton->SetVeto(8317.44);
    singlet->push_back(parton);
    part   = new Particle(0,Flavour(kf::u).Bar(),Vec4D(45.6,0.,0.,-45.6),'F');
    part->SetNumber(0);
    part->SetStatus(part_status::decayed);
    parton = new Parton(part,pst::FS);
    blob->AddToInParticles(part);
    parton->SetStart(8317.44);
    parton->SetVeto(8317.44);
    singlet->push_back(parton);
    allsinglets->push_back(singlet);
    if (shower.PerformShowers(allsinglets)) shower.ExtractPartons(blob);
    

    cout<<(*blob);
    
    multi.Insert(singlet->size());
    for (ASiter asit=allsinglets->begin();asit!=allsinglets->end();asit++) {
      Singlet * sing = (*asit);
      if (sing) { delete sing; sing = NULL; }
    }
    allsinglets->clear();
    singlet = NULL;
    delete blob; blob = NULL;
  }
  multi.Output(string("Multi.dat"));
  delete allsinglets;
}

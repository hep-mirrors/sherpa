#include <iostream>
#include <algorithm>
#include "ATOOLS/Org/Run_Parameter.H"
#include "MODEL/Main/Standard_Model.H"
#include "PDF/Main/ISR_Handler.H"
#include "PDF/Main/Intact.H"

#include "APACIC++/Main/Apacic.H"
#include "APACIC++/Main/Tree.H"
#include "ATOOLS/Math/Random.H"
#include "APACIC++/Showers/Final_State_Shower.H"
#include "ATOOLS/Math/Histogram.H"

using namespace ATOOLS;
using namespace MODEL;
using namespace PDF;
using namespace APACIC;

class Thrust {
  std::vector<ATOOLS::Vec3D> moms;
  ATOOLS::Vec3D thrustaxis;
  double thrust;


  void   FindAxis();
  double SumP(const std::vector<ATOOLS::Vec3D> &);
  double SumNP(const std::vector<ATOOLS::Vec3D> &,const ATOOLS::Vec3D &);
  double CalculateThrust();
  Vec3D  NewAxis(const std::vector<ATOOLS::Vec3D> & p,const ATOOLS::Vec3D & ref);
  unsigned int ipow(int base, int exponent);
 public:
  int operator()(const ATOOLS::Vec3D & lhs,const ATOOLS::Vec3D & rhs) {
    return (lhs.Sqr()>rhs.Sqr()); 
  }

  double CalcThrust(const Particle_List & pl);
};


double Thrust::CalcThrust(const Particle_List & pl)
{
  moms.clear(); 
  for (Particle_List::const_iterator pit=pl.begin();pit!=pl.end();++pit) {
    moms.push_back(Vec3D((*pit)->Momentum()));
  }
  FindAxis();
  return thrust;
}

void Thrust::FindAxis()
{
  const unsigned int startaxes=4;
  const unsigned int maxidentaxes=2;
  const double accuracy = 1.e-4;

  thrust = 0.;

  Vec3D maxthrustaxis, lastaxis, curraxis;
  double maxthrust=0., lastthrust , currthrust;
  unsigned int min_generators = startaxes < moms.size() ? startaxes : moms.size();

  std::vector<Vec3D> initialaxes;
  int addsign;

  initialaxes.clear();
  std::sort(moms.begin(),moms.end(),Thrust());
    
  for(unsigned int i=1;i<=ipow(2,min_generators-1);++i) {
    Vec3D axis;
    for(unsigned int j=1;j<=min_generators;++j) {
      addsign = -1;
      if (ipow(2,j)*((i+ipow(2,j-1)-1)/ipow(2,j)) >= i) addsign = 1;
      axis = axis+addsign*moms[j-1];
    }
    initialaxes.push_back(axis);
  }
  // sort the initial axes with respect to their size ( which corresponds 
  // to their thrust because of the common denominator) 
  std::sort(initialaxes.begin(),initialaxes.end(), Thrust());
  for(unsigned int j=0;j<initialaxes.size();j++) 
    initialaxes[j] = initialaxes[j]/initialaxes[j].Abs();

  unsigned int ident = 0;
  double sump        = SumP(moms);
  maxthrust          = 0.;
  for(unsigned int j=0; (j<initialaxes.size()) && (ident<maxidentaxes); j++) {
    curraxis         = initialaxes[j];
    currthrust       = SumNP(moms,curraxis)/sump;
    lastthrust       = 0.;
    while (currthrust > lastthrust+accuracy) {
      lastthrust     = currthrust;
      lastaxis       = curraxis;
      curraxis       = NewAxis(moms,curraxis);
      currthrust     = SumNP(moms,curraxis)/sump;
    }
    // if it gets worse then skip this axis alltogether
    if (lastthrust < maxthrust-accuracy) break;
    // if it is a better solution then keep this one
    if (lastthrust > maxthrust+accuracy) {
      ident          = 0;
      maxthrustaxis  = lastaxis;
      maxthrust      = lastthrust;
    }
    ident++;
  }
  thrustaxis = maxthrustaxis; 
  thrust     = maxthrust; 
}


double Thrust::SumP(const std::vector<Vec3D> & p) { 
  double sum_p = 0.;
  for (unsigned int i=0; i<p.size(); i++) sum_p += p[i].Abs();
  return sum_p;
}

double Thrust::SumNP(const std::vector<Vec3D> & p,const Vec3D & n) { 
  double sum_np = 0.;
  for (unsigned int i=0; i<p.size(); i++) sum_np += dabs(p[i]*n);
  return sum_np;
}

Vec3D Thrust::NewAxis(const std::vector<Vec3D> & p,const Vec3D & ref) {
  Vec3D nextref = Vec3D(0.,0.,0.);
  int addsign;
  for (unsigned int i=0;i<p.size();++i) {
    addsign = 1;
    if (ref*p[i]<0.) addsign = -1;
    nextref = nextref+addsign*p[i];
  }
  return nextref/nextref.Abs();  
}


unsigned int Thrust::ipow(int base,int exponent) { 
  int result=1;
  if (exponent>0) 
    do {
      result*=base;
    } while(--exponent);
  return result;
}


// ======================================================================
//         APACIC 2.0 demonstration program  (LEP1)
// ======================================================================

int main(int argc,char **argv)
{
  std::cout<<" APACIC 2.0 demonstration program "<<std::endl;

  // initialize the framework
  //  rpa      (general runtime paramters)
  //  particle (physics framework)
  //  dataread (shower runtime parameters)
  //  model    (running alphas and alphaqed)
  //  isr      (pdf or electron structure function)
  //  apacic   (shower interface)
  ATOOLS::ParticleInit("./"); 
  rpa.Init("./","Run.dat",argc,argv);
  rpa.gen.SetEcms(91.2);

  Data_Read     * dataread     = new Data_Read("Shower.dat");
  Model_Base    * model        = new Standard_Model("./","Model.dat");

  ISR_Base ** isrbases = new ISR_Base*[2];
  isrbases[0] = new Intact(Flavour(kf_e));     
  isrbases[1] = new Intact(Flavour(kf_e).Bar());
  ISR_Handler * isr = new ISR_Handler(isrbases);

  // init analysis
  Histogram histo_thrust(0,0.,.5,50);   // (1 - thrust)
  Histogram histo_multi(0,-0.5,20.5,21);
  Thrust analysis;
  Jet_Finder jf(rpa.gen.Ycut(),1,false);
  Apacic * apacic = new APACIC::Apacic(isr,model,&jf,dataread);
  // end initialization


  // get pointer to final state tree
  Tree  * tree = apacic->FinTree();

  double E = 0.5*rpa.gen.Ecms();
  double rscale = rpa.gen.Ecms();
  double scale = rscale*rscale;

  std::cout<<" Energy : "<<E<<std::endl;

  // start event loop
  for (int i=1;i<=rpa.gen.NumberOfEvents();i++) {
    if (i%2500==0) {
      msg_Out()<<" "<<i<<" th event "<<std::endl;
    }

    // define hard momenta
    double costheta=2.*ran.Get()-1.;
    double sintheta=sqrt(1.-costheta*costheta);
    double phi=2.*M_PI*ran.Get();

    Flavour quark=Flavour((kf_code)(1+int(ran.Get()*4.)));

    Vec3D dir(sin(phi)*sintheta,cos(phi)*sintheta,costheta);

    // reset shower history
    tree->Reset();

    // define initial condition for final state shower evolution
    Knot * mo = tree->NewKnot();
    Knot * k2 = tree->NewKnot();
    Knot * k3 = tree->NewKnot();

    mo->z       = 0.5;
    *(mo->part) = Particle(1,Flavour(kf_photon),Vec4D(rscale,0,0,0));
    mo->E2    = scale;
    mo->t    = scale;  
    
    mo->stat    = 1;  
  
    mo->left  = k2;
    mo->right = k3;

    *(k2->part) = Particle(2,quark,Vec4D(E,E*dir));
    k2->part->SetFlow(1,-1);
    *(k3->part) = Particle(3,quark.Bar(),Vec4D(E,(-E)*dir));
    k3->part->SetFlow(2,k2->part->GetFlow(1));

    k2->stat = 3;
    k2->E2   = E*E;
    k2->t    = scale;  
    k2->prev = mo;

    k3->stat = 3;
    k3->E2   = E*E;
    k3->t    = scale;  
    k3->prev = mo;

    // call shower routines
    apacic->PerformShowers();

    // extract parton list
    Particle_List * pl = new Particle_List();
    apacic->FinShower()->ExtractPartons(tree->GetRoot(),NULL,NULL,pl);
    double thrust = analysis.CalcThrust(*pl);
    int npart = pl->size();
    delete pl;

    // fill histograms
    histo_thrust.Insert(1.-thrust);
    histo_multi.Insert(npart);
  }

  // output multiplicity and thrust distribution
  histo_thrust.Finalize();
  histo_thrust.Output("one-thrust.dat");
  histo_multi.Finalize();
  histo_multi.Output("multi.dat");

  // finalization
  delete apacic;
  delete isr;
  delete model;
  delete dataread;

  // 

  return 0;
}

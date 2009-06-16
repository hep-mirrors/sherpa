#include <iostream>
#include "ATOOLS/Org/Run_Parameter.H"
#include "MODEL/Main/Standard_Model.H"
#include "PDF/Main/ISR_Handler.H"
#include "PDF/Main/Structure_Function.H"
#include "PDF/CTEQ/CTEQ6_Fortran_Interface.H"

#include "APACIC++/Main/Apacic.H"
#include "APACIC++/Main/Tree.H"
#include "ATOOLS/Math/Random.H"
#include "XS_Drell_Yan.H"

using namespace ATOOLS;
using namespace MODEL;
using namespace PDF;
using namespace APACIC;
using namespace EXTRAXS;

// simple example process
class Hard_Process {
  Flavour proton;
  Flavour antiproton;
  PDF_Base  * pdf1;
  PDF_Base  * pdf2;
  XS_Base  ** xsecs;
  int     nxsecs;

  double  *sigmas;
  double  sigma_max;
  double  smin, smax, sbeam;
  double  ymin, ymax;

  double  sprime, y, costheta;
  double  x1, x2;

  Flavour parton1, parton2;

  long    ntotal;
  double  sumxs;
  double  sumxs2;
 public:
  Hard_Process();
  double CrossSection();
  void Dice();

  PDF_Base * GetPDF1() { return pdf1; }
  PDF_Base * GetPDF2() { return pdf2; }
  const Flavour & GetBeam1() { return proton; }
  const Flavour & GetBeam2() { return antiproton; }

  double GetX1() { return x1; }
  double GetX2() { return x2; }
  double GetCosTheta() { return costheta; }
  const Flavour & Parton1() { return parton1; }
  const Flavour & Parton2() { return parton2; }
  
  double GetTotal() { return sumxs/double(ntotal); }
};


Hard_Process::Hard_Process() :
  proton(kf_p_plus),antiproton(proton.Bar()),
  nxsecs(5),
  ntotal(0),sumxs(0.),sumxs2(0.)
{
  pdf1 = new CTEQ6_Fortran_Interface(proton,"cteq6l1",0,"grid");
  pdf2 = new CTEQ6_Fortran_Interface(antiproton,"cteq6l1",0,"grid");

  // initialize Drell-Yan processes
  xsecs = new XS_Base*[nxsecs];
  sigmas = new double[2*nxsecs];
  for (int i=0; i<nxsecs; ++i) {
    Flavour electron=Flavour(kf_e);
    Flavour quark=Flavour((kf_code)(i+1));
    Flavour flavours[4];
    flavours[0]=quark;
    flavours[1]=quark.Bar();
    flavours[2]=electron;
    flavours[3]=electron.Bar();
    xsecs[i]=new XS_ee_ffbar(2,2,flavours);  
  }

  smin   = sqr(66.);
  smax   = sqr(116.);
  sbeam  = sqr(rpa.gen.Ecms());

  x1=x2= 91.2/rpa.gen.Ecms();
  sprime = x1*x2*sbeam;
  ymin   = .5*log(sprime/sbeam);
  ymax   = -ymin;
  costheta  = 1.;
  sigma_max = CrossSection();
}


void Hard_Process::Dice()
{
  double sigma=0.;
  for (;;) {
    // dice point (sprime, y, costheta)
    sprime = (smax-smin)*ran.Get() + smin;

    double tau = sprime/sbeam;
    double rtau = sqrt(tau);
    ymin = .5*log(tau);
    ymax = -ymin;
    double y = (ymax-ymin)*ran.Get() + ymin;

    x1 = rtau*exp(y);
    x2 = rtau*exp(-y);

    costheta = 2.*ran.Get() - 1.;

    // calc cross section
    sigma = CrossSection();

    ntotal++;
    sumxs+=sigma;
    sumxs2+=sqr(sigma);
    
    // reject or accept point
    if (sigma>sigma_max) {
      std::cout<<" WARNING max to small \n";
      sigma_max=sigma;
    }
    if (sigma/sigma_max > ran.Get()) break;
  }

  // select channel
  double disc = sigma*ran.Get();
  int i=0;
  while (sigmas[i]<disc) {
    ++i;
    if (i>=nxsecs*2) {
      i=0; std::cout<<" WARNING: Hard_Process::Dice() cross section selection failed \n";
      break;
    }
  }
  if (i<nxsecs) {
    parton1=xsecs[i]->Flavours()[0];
    parton2=xsecs[i]->Flavours()[1];
  }
  else if (i<2*nxsecs) {
    parton1=xsecs[i-nxsecs]->Flavours()[1];
    parton2=xsecs[i-nxsecs]->Flavours()[0];
  }
}

double Hard_Process::CrossSection()
{
  pdf1->Calculate(x1,sprime);
  pdf2->Calculate(x2,sprime);
  double s = sprime;
  double t = .5*sprime*(costheta-1.);
  double u = -.5*sprime*(costheta+1.);
  double sigma=0.;
  for (int i=0;i<nxsecs;++i) {
    sigmas[i]=sigma+=pdf1->GetXPDF(xsecs[i]->Flavours()[0]) *
      pdf2->GetXPDF(xsecs[i]->Flavours()[1]) * 
      (*xsecs[i])(s,t,u);
  }
  for (int i=0;i<nxsecs;++i) {
    sigmas[i+nxsecs]=sigma+=pdf1->GetXPDF(xsecs[i]->Flavours()[1]) *
      pdf2->GetXPDF(xsecs[i]->Flavours()[0]) * 
      (*xsecs[i])(s,u,t);    
  }
  sigma/= sqr(sprime)*32.*M_PI;       
  sigma*= (ymax-ymin)*(smax-smin)*2;  
  return sigma;
}

// ======================================================================
//              ApacicTevatron -- simple testprogram 
// ======================================================================
int main(int argc,char **argv)
{
  std::cout<<" APACIC 2.0 Test program "<<std::endl;

  // initialize the framework
  //  particle (physics framework)
  //  rpa      (general runtime parameters)
  //  dataread (shower runtime parameters)
  //  model    (running alphas and alphaqed (Standard Model))
  //  isr      (parton distribution functions)
  //  histo    (simple analysis)
  //  apacic   (shower interface)

  ATOOLS::ParticleInit("./"); 
  rpa.Init("./","Run.dat",argc,argv);
  rpa.gen.SetEcms(1800.);  // Tevatron Run I
  
  Data_Read     * dataread     = new Data_Read("Shower.dat");
  Model_Base    * model        = new Standard_Model("./","Model.dat");
  MODEL::s_model=p_model;

  Hard_Process * hardprocess = new Hard_Process();

  ISR_Base ** isrbases = new ISR_Base*[2];
  isrbases[0] = new Structure_Function(hardprocess->GetPDF1(),hardprocess->GetBeam1());
  isrbases[1] = new Structure_Function(hardprocess->GetPDF2(),hardprocess->GetBeam2());
  ISR_Handler * isr = new ISR_Handler(isrbases);

  Histogram histo_pt(0,0.,200.,50);
  Histogram histo_mass(0,60.,120.,60);

  Jet_Finder jf(rpa.gen.Ycut(),4,false);
  Apacic * apacic   = new APACIC::Apacic(isr,model,&jf,dataread);
  Tree   * fintree  = apacic->FinTree();
  Tree  ** initrees = apacic->IniTrees();

  std::cout<<" CMS Energy : "<<rpa.gen.Ecms()<<std::endl;
  // end initialization

  // event loop
  for (int i=1;i<=rpa.gen.NumberOfEvents();i++) {
    if (i%2500==0) {
      msg_Out()<<" "<<i<<" th event "<<std::endl;
    }

    // determine kinematics of hard process
    hardprocess->Dice();
    double x1  = hardprocess->GetX1();
    double x2  = hardprocess->GetX2();
    double tau = x1*x2;
    double E   = .5*sqrt(tau)*rpa.gen.Ecms();
    double costheta = hardprocess->GetCosTheta();
    double sintheta = sqrt(1.-costheta*costheta);
    double phi = 2.*M_PI*ran.Get();

    double rscale = 2.*E;
    double scale  = sqr(rscale);
    double iscale = sqr(rpa.gen.Ecms());
 
    // define cms momenta and particle flavour
    Vec4D mom0(E,0,0,E);
    Vec4D mom1(E,0,0,-E);

    Flavour electron = Flavour(kf_e);
    Flavour quark    = hardprocess->Parton1();

    Vec3D dir(sin(phi)*sintheta,cos(phi)*sintheta,costheta);

    // reset shower history
    apacic->PrepareTrees();

    // define initial condition of shower evolution
    Knot * mo = fintree->NewKnot();
    Knot * k2 = fintree->NewKnot();
    Knot * k3 = fintree->NewKnot();

    // final state part
    mo->z       = 0.5;
    *(mo->part) = Particle(1,Flavour(kf_photon),Vec4D(rscale,0,0,0));
    mo->E2    = scale;
    mo->t    = scale;  
    
    mo->stat    = 1;  
  
    mo->left  = k2;
    mo->right = k3;

    *(k2->part) = Particle(2,electron,Vec4D(E,E*dir));
    *(k3->part) = Particle(3,electron.Bar(),Vec4D(E,(-E)*dir));

    k2->stat = 3;
    k2->E2   = E*E;
    k2->t    = scale;  
    k2->prev = mo;

    k3->stat = 3;
    k3->E2   = E*E;
    k3->t    = scale;  
    k3->prev = mo;

    // initial state part

    Knot * k0 = initrees[0]->NewKnot();
    Knot * k1 = initrees[1]->NewKnot();
    
    *(k0->part) = Particle(4,quark,Vec4D(E,E*Vec3D::ZVEC));
    *(k1->part) = Particle(5,quark.Bar(),Vec4D(E,(-E)*Vec3D::ZVEC));
    if (quark.IsAnti()) {
      k1->part->SetFlow(1,-1);
      k0->part->SetFlow(2,k1->part->GetFlow(1));
    }
    else {
      k0->part->SetFlow(1,-1);
      k1->part->SetFlow(2,k0->part->GetFlow(1));
    }

    k0->stat = 3;
    k0->E2   = E*E;
    k0->t    = -iscale; 
    k0->x    = x1;

    k1->stat = 3;
    k1->E2   = E*E;
    k1->t    = -iscale;  
    k1->x    = x2;

    // call shower routines
    apacic->PerformShowers(true,true,-1);

    Vec4D mom_boson = fintree->GetRoot()->part->Momentum();
    // fill histograms
    histo_pt.Insert(mom_boson.PPerp());
    histo_mass.Insert(mom_boson.Mass());
  }

  // output pt and mass distribution of Z-boson
  histo_pt.Finalize();
  histo_pt.Output("pt.dat");
  histo_mass.Finalize();
  histo_mass.Output("mass.dat");

  // output total cross section
  std::cout<<" total xs "<<hardprocess->GetTotal()*rpa.Picobarn()<<"pb \n";

  // finalization
  delete apacic;
  delete isr;
  delete model;
  delete dataread;

  return 0;
}

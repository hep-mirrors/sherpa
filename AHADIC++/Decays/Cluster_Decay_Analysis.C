#include "AHADIC++/Decays/Cluster_Decay_Analysis.H"
#include "ATOOLS/Org/Run_Parameter.H"
#include "ATOOLS/Math/Vector.H"

using namespace AHADIC;
using namespace ATOOLS;
using namespace std;

Cluster_Decay_Analysis::Cluster_Decay_Analysis() 
{
  return;
  m_histograms[string("pi+_Number")]                 = new Histogram(0,0.,20.,20);
  m_histograms[string("pi-_Number")]                 = new Histogram(0,0.,20.,20);
  m_histograms[string("pi0_Number")]                 = new Histogram(0,0.,20.,20);
  m_histograms[string("K+_Number")]                  = new Histogram(0,0.,20.,20);
  m_histograms[string("K-_Number")]                  = new Histogram(0,0.,20.,20);
  m_histograms[string("K0_Number")]                  = new Histogram(0,0.,20.,20);
  m_histograms[string("K0_Bar_Number")]              = new Histogram(0,0.,20.,20);
  m_histograms[string("eta_Number")]                 = new Histogram(0,0.,20.,20);
  m_histograms[string("etaPrime_Number")]            = new Histogram(0,0.,20.,20);
  
  m_histograms[string("rho+_Number")]                = new Histogram(0,0.,20.,20);
  m_histograms[string("rho-_Number")]                = new Histogram(0,0.,20.,20);
  m_histograms[string("rho0_Number")]                = new Histogram(0,0.,20.,20);
  m_histograms[string("KStar+_Number")]              = new Histogram(0,0.,20.,20);
  m_histograms[string("KStar-_Number")]              = new Histogram(0,0.,20.,20);
  m_histograms[string("KStar0_Number")]              = new Histogram(0,0.,20.,20);
  m_histograms[string("KStar0_Bar_Number")]          = new Histogram(0,0.,20.,20);
  m_histograms[string("omega_Number")]               = new Histogram(0,0.,20.,20);
  m_histograms[string("phi_Number")]                 = new Histogram(0,0.,20.,20);
  
  m_histograms[string("x_p_Pseudoscalars")]          = new Histogram(0,0.,1.,100);
  m_histograms[string("x_p_Vectors")]                = new Histogram(0,0.,1.,100);
  m_histograms[string("x_p_C-Hadrons")]              = new Histogram(0,0.,1.,100);
  m_histograms[string("x_p_B-Mesons")]               = new Histogram(0,0.,1.,100);
  m_histograms[string("x_p_Bstar-Mesons")]           = new Histogram(0,0.,1.,100);
  m_histograms[string("x_E_Bstar-Mesons")]           = new Histogram(0,0.,1.,100);
  m_histograms[string("x_E_B-Mesons")]               = new Histogram(0,0.,1.,100);
  m_histograms[string("x_E_B-Quarks")]               = new Histogram(0,0.,1.,100);
  m_histograms[string("x_E_B-Mesons_L")]             = new Histogram(0,0.,1.,100);
  m_histograms[string("x_E_B-Quarks_L")]             = new Histogram(0,0.,1.,100);

  m_histograms[string("Q(Lambda_Lambda)")]           = new Histogram(0,0.,5.5,11);
}



Cluster_Decay_Analysis::~Cluster_Decay_Analysis()
{ 
  return;
  Histogram * histo;
  string name;
  for (map<string,Histogram *>::iterator hit=m_histograms.begin();
       hit!=m_histograms.end();hit++) {
    histo = hit->second;
    name  = string("Fragmentation_Analysis/")+hit->first+string(".dat");
    histo->Output(name);
    delete histo;
  }
  m_histograms.clear();
}

void Cluster_Decay_Analysis::AnalyseThis(Blob * blob)
{
  return;
  int Npiplus(0),Npiminus(0),Npi0(0),NKplus(0),NKminus(0),NK0(0),NK0b(0),Neta(0),Netaprime(0),
    Nrhoplus(0),Nrhominus(0),Nrho0(0),NKstarplus(0),NKstarminus(0),NKstar0(0),NKstar0b(0),
    Nomega(0),Nphi(0);
  Particle * part;
  int kfc, LambdaCount(0), LambdaP(0), LambdaM(0);
  Vec4D QLambda(0.,0.,0.,0.);
  double x, max_x = 0., max_xB = 0.;
  for (int i=0;i<blob->NInP();i++) {
    part = blob->InParticle(i);
    kfc  = int(part->Flav());
    x    = 2.*part->Momentum()[0]/rpa->gen.Ecms();
    if (x>max_x) max_x = x;
    switch (kfc) {
    case 5:
    case -5:
      if (x>max_xB) max_xB = x;
      m_histograms[string("x_E_B-Quarks")]->Insert(x);
      break;
    }
  }
  if (max_x==max_xB) m_histograms[string("x_E_B-Quarks_L")]->Insert(max_x);

  max_x  = max_xB = 0.;

  for (int i=0;i<blob->NOutP();i++) {
    part = blob->OutParticle(i);
    x    = 2.*part->Momentum()[0]/rpa->gen.Ecms();
    if (x>max_x) max_x = x;
    kfc  = int(part->Flav());
    switch (kfc) {
    case 111:
      Npi0++;
      m_histograms[string("x_p_Pseudoscalars")]->Insert(2.*Vec3D(part->Momentum()).Abs()/rpa->gen.Ecms());
      break;
    case 211:
      Npiplus++;
      m_histograms[string("x_p_Pseudoscalars")]->Insert(2.*Vec3D(part->Momentum()).Abs()/rpa->gen.Ecms());
      break;
    case -211:
      Npiminus++;
      m_histograms[string("x_p_Pseudoscalars")]->Insert(2.*Vec3D(part->Momentum()).Abs()/rpa->gen.Ecms());
      break;
    case 221:
      Neta++;
      m_histograms[string("x_p_Pseudoscalars")]->Insert(2.*Vec3D(part->Momentum()).Abs()/rpa->gen.Ecms());
      break;
    case 311:
      NK0++;
      m_histograms[string("x_p_Pseudoscalars")]->Insert(2.*Vec3D(part->Momentum()).Abs()/rpa->gen.Ecms());
      break;
    case -311:
      NK0b++;
      m_histograms[string("x_p_Pseudoscalars")]->Insert(2.*Vec3D(part->Momentum()).Abs()/rpa->gen.Ecms());
      break;
    case 321:
      NKplus++;
      m_histograms[string("x_p_Pseudoscalars")]->Insert(2.*Vec3D(part->Momentum()).Abs()/rpa->gen.Ecms());
      break;
    case -321:
      NKminus++;
      m_histograms[string("x_p_Pseudoscalars")]->Insert(2.*Vec3D(part->Momentum()).Abs()/rpa->gen.Ecms());
      break;
    case 331:
      Netaprime++;
      m_histograms[string("x_p_Pseudoscalars")]->Insert(2.*Vec3D(part->Momentum()).Abs()/rpa->gen.Ecms());
      break;
    case 113:
      Nrho0++;
      m_histograms[string("x_p_Vectors")]->Insert(2.*Vec3D(part->Momentum()).Abs()/rpa->gen.Ecms());
      break;
    case 213:
      Nrhoplus++;
      m_histograms[string("x_p_Vectors")]->Insert(2.*Vec3D(part->Momentum()).Abs()/rpa->gen.Ecms());
      break;
    case -213:
      Nrhominus++;
      m_histograms[string("x_p_Vectors")]->Insert(2.*Vec3D(part->Momentum()).Abs()/rpa->gen.Ecms());
      break;
    case 223:
      Nomega++;
      m_histograms[string("x_p_Vectors")]->Insert(2.*Vec3D(part->Momentum()).Abs()/rpa->gen.Ecms());
      break;
    case 313:
      NKstar0++;
      m_histograms[string("x_p_Vectors")]->Insert(2.*Vec3D(part->Momentum()).Abs()/rpa->gen.Ecms());
      break;
    case -313:
      NKstar0b++;
      m_histograms[string("x_p_Vectors")]->Insert(2.*Vec3D(part->Momentum()).Abs()/rpa->gen.Ecms());
      break;
    case 323:
      NKstarplus++;
      m_histograms[string("x_p_Vectors")]->Insert(2.*Vec3D(part->Momentum()).Abs()/rpa->gen.Ecms());
      break;
    case -323:
      NKstarminus++;
      m_histograms[string("x_p_Vectors")]->Insert(2.*Vec3D(part->Momentum()).Abs()/rpa->gen.Ecms());
      break;
    case 333:
      Nphi++;
      m_histograms[string("x_p_Vectors")]->Insert(2.*Vec3D(part->Momentum()).Abs()/rpa->gen.Ecms());
      break;
    case 411:
    case 413:
    case 421:
    case 423:
    case 431:
    case 433:
      m_histograms[string("x_p_C-Hadrons")]->Insert(2.*Vec3D(part->Momentum()).Abs()/rpa->gen.Ecms());
      break;
    case 511:
    case 521:
      //case 531:
      //      if (2.*Vec3D(part->Momentum()).Abs()/rpa->gen.Ecms()<0.89 &&
      //	  2.*Vec3D(part->Momentum()).Abs()/rpa->gen.Ecms()>0.88) 
      //	cout<<(*blob)<<endl<<endl
      //	    <<"#################################################################################"<<endl<<endl;
      m_histograms[string("x_p_B-Mesons")]->Insert(2.*Vec3D(part->Momentum()).Abs()/rpa->gen.Ecms());
      m_histograms[string("x_E_B-Mesons")]->Insert(x);
      if (x>max_xB) max_xB = x;
      break;
    case 513:
    case 523:
    case 533:
      m_histograms[string("x_p_Bstar-Mesons")]->Insert(2.*Vec3D(part->Momentum()).Abs()/rpa->gen.Ecms());
      m_histograms[string("x_E_Bstar-Mesons")]->Insert(x);
      break;
    case 3122:
      QLambda += part->Momentum();
      LambdaCount++; LambdaP++;
      break;
    case -3122:
      QLambda -= part->Momentum();
      LambdaCount++; LambdaM++;
      break;
    }
  }    
  if (LambdaCount==2 && LambdaP==1 && LambdaM==1) 
    m_histograms[string("Q(Lambda_Lambda)")]->Insert(sqrt(-QLambda.Abs2()));
  if (max_x==max_xB)   m_histograms[string("x_E_B-Mesons_L")]->Insert(max_x);

  m_histograms[string("pi+_Number")]->Insert(Npiplus);
  m_histograms[string("pi-_Number")]->Insert(Npiminus);
  m_histograms[string("pi0_Number")]->Insert(Npi0);
  m_histograms[string("K+_Number")]->Insert(NKplus);
  m_histograms[string("K-_Number")]->Insert(NKminus);
  m_histograms[string("K0_Number")]->Insert(NK0);
  m_histograms[string("K0_Bar_Number")]->Insert(NK0b);
  m_histograms[string("eta_Number")]->Insert(Neta);
  m_histograms[string("etaPrime_Number")]->Insert(Netaprime);
    
  m_histograms[string("rho+_Number")]->Insert(Nrhoplus);
  m_histograms[string("rho-_Number")]->Insert(Nrhominus);
  m_histograms[string("rho0_Number")]->Insert(Nrho0);
  m_histograms[string("KStar+_Number")]->Insert(NKstarplus);
  m_histograms[string("KStar-_Number")]->Insert(NKstarminus);
  m_histograms[string("KStar0_Number")]->Insert(NKstar0);
  m_histograms[string("KStar0_Bar_Number")]->Insert(NKstar0b);
  m_histograms[string("omega_Number")]->Insert(Nomega);
  m_histograms[string("phi_Number")]->Insert(Nphi);
}

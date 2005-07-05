#include "Cluster_Decay_Analysis.H"
#include "Run_Parameter.H"

using namespace AHADIC;
using namespace ATOOLS;
using namespace std;

Cluster_Decay_Analysis::Cluster_Decay_Analysis() 
{
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
  m_histograms[string("x_p_B-Hadrons")]              = new Histogram(0,0.,1.,100);
  m_histograms[string("x_p_C-Hadrons")]              = new Histogram(0,0.,1.,100);
  m_histograms[string("x_E_B-Hadrons")]              = new Histogram(0,0.,1.,100);
}



Cluster_Decay_Analysis::~Cluster_Decay_Analysis()
{ 
  Histogram * histo;
  string name;
  for (map<string,Histogram *>::iterator hit=m_histograms.begin();
       hit!=m_histograms.end();hit++) {
    histo = hit->second;
    name  = string("Analysis/")+hit->first+string(".dat");
    histo->Output(name);
    delete histo;
  }
  m_histograms.clear();
}

void Cluster_Decay_Analysis::AnalyseThis(Blob * blob)
{
  int Npiplus=0,Npiminus=0,Npi0=0,NKplus=0,NKminus=0,NK0=0,NK0b=0,Neta=0,Netaprime=0,NPS=0;
  int Nrhoplus=0,Nrhominus=0,Nrho0=0,NKstarplus=0,NKstarminus=0,NKstar0=0,NKstar0b=0,
    Nomega=0,Nphi=0,NV=0;
  Particle * part;
  int kfc;
  for (int i=0;i<blob->NOutP();i++) {
    part = blob->OutParticle(i);
    kfc  = int(part->Flav());
    switch (kfc) {
    case 111:
      Npi0++;
      m_histograms[string("x_p_Pseudoscalars")]->Insert(2.*Vec3D(part->Momentum()).Abs()/rpa.gen.Ecms());
      break;
    case 211:
      Npiplus++;
      m_histograms[string("x_p_Pseudoscalars")]->Insert(2.*Vec3D(part->Momentum()).Abs()/rpa.gen.Ecms());
      break;
    case -211:
      Npiminus++;
      m_histograms[string("x_p_Pseudoscalars")]->Insert(2.*Vec3D(part->Momentum()).Abs()/rpa.gen.Ecms());
      break;
    case 221:
      Neta++;
      m_histograms[string("x_p_Pseudoscalars")]->Insert(2.*Vec3D(part->Momentum()).Abs()/rpa.gen.Ecms());
      break;
    case 311:
      NK0++;
      m_histograms[string("x_p_Pseudoscalars")]->Insert(2.*Vec3D(part->Momentum()).Abs()/rpa.gen.Ecms());
      break;
    case -311:
      NK0b++;
      m_histograms[string("x_p_Pseudoscalars")]->Insert(2.*Vec3D(part->Momentum()).Abs()/rpa.gen.Ecms());
      break;
    case 321:
      NKplus++;
      m_histograms[string("x_p_Pseudoscalars")]->Insert(2.*Vec3D(part->Momentum()).Abs()/rpa.gen.Ecms());
      break;
    case -321:
      NKminus++;
      m_histograms[string("x_p_Pseudoscalars")]->Insert(2.*Vec3D(part->Momentum()).Abs()/rpa.gen.Ecms());
      break;
    case 331:
      Netaprime++;
      m_histograms[string("x_p_Pseudoscalars")]->Insert(2.*Vec3D(part->Momentum()).Abs()/rpa.gen.Ecms());
      break;
    case 113:
      Nrho0++;
      m_histograms[string("x_p_Vectors")]->Insert(2.*Vec3D(part->Momentum()).Abs()/rpa.gen.Ecms());
      break;
    case 213:
      Nrhoplus++;
      m_histograms[string("x_p_Vectors")]->Insert(2.*Vec3D(part->Momentum()).Abs()/rpa.gen.Ecms());
      break;
    case -213:
      Nrhominus++;
      m_histograms[string("x_p_Vectors")]->Insert(2.*Vec3D(part->Momentum()).Abs()/rpa.gen.Ecms());
      break;
    case 223:
      Nomega++;
      m_histograms[string("x_p_Vectors")]->Insert(2.*Vec3D(part->Momentum()).Abs()/rpa.gen.Ecms());
      break;
    case 313:
      NKstar0++;
      m_histograms[string("x_p_Vectors")]->Insert(2.*Vec3D(part->Momentum()).Abs()/rpa.gen.Ecms());
      break;
    case -313:
      NKstar0b++;
      m_histograms[string("x_p_Vectors")]->Insert(2.*Vec3D(part->Momentum()).Abs()/rpa.gen.Ecms());
      break;
    case 323:
      NKstarplus++;
      m_histograms[string("x_p_Vectors")]->Insert(2.*Vec3D(part->Momentum()).Abs()/rpa.gen.Ecms());
      break;
    case -323:
      NKstarminus++;
      m_histograms[string("x_p_Vectors")]->Insert(2.*Vec3D(part->Momentum()).Abs()/rpa.gen.Ecms());
      break;
    case 333:
      Nphi++;
      m_histograms[string("x_p_Vectors")]->Insert(2.*Vec3D(part->Momentum()).Abs()/rpa.gen.Ecms());
      break;
    case 411:
    case 413:
    case 421:
    case 423:
    case 431:
    case 433:
      m_histograms[string("x_p_C-Hadrons")]->Insert(2.*Vec3D(part->Momentum()).Abs()/rpa.gen.Ecms());
      break;
    case 511:
    case 513:
    case 521:
    case 523:
    case 531:
    case 533:
      m_histograms[string("x_p_B-Hadrons")]->Insert(2.*Vec3D(part->Momentum()).Abs()/rpa.gen.Ecms());
      m_histograms[string("x_E_B-Hadrons")]->Insert(part->Momentum()[0]/rpa.gen.Ecms());
      break;
    }
  }    

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

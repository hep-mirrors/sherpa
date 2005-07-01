#include "Cluster_Decay_Handler.H"
#include "Cluster_Part.H"
#include "Hadron_Part.H"


using namespace AHADIC;
using namespace ATOOLS;
using namespace std;

Cluster_Decay_Handler::Cluster_Decay_Handler(Cluster_Transformer * transformer,bool ana) :
  m_cdm(cdm::QoverM_Isotropic), 
  p_decayer(NULL), 
  p_transformer(transformer),
  p_partlist(new Part_List), 
  m_analyse(ana)
{ 
  Cluster_Part * cp;
  Hadron_Part * hp;
  switch (int(m_cdm/10)) {
  case 2:
    cp = new FourFermion();
  case 1:
  default:
    cp = new QoverM();
  }
  switch (int(m_cdm%10)) {
  case 2:
    hp = new Retain();
  case 1:
  default:
    hp = new Isotropic();
  }
  p_decayer = new Cluster_Decayer_Base(cp,hp);
  InitializeAnalysis();
}



Cluster_Decay_Handler::~Cluster_Decay_Handler()
{ 
  if (m_analyse) {
    Histogram * histo;
    string name;
    for (map<string,Histogram *>::iterator hit=m_histograms.begin();
	 hit!=m_histograms.end();hit++) {
      histo = hit->second;
      name  = hit->first+string(".dat");
      histo->Output(name);
      delete histo;
    }
    m_histograms.clear();
  }
  if (p_decayer)  { delete p_decayer;  p_decayer=NULL;  }
  if (p_partlist) { delete p_partlist; p_partlist=NULL; }
}

void Cluster_Decay_Handler::DecayClusters(Cluster_List * clusters,Blob * blob)
{
  p_partlist->clear();
  msg.Tracking()<<"Decay the clusters ------------------------------------------------"<<endl;
  //<<(*clusters);
  Cluster_Iterator cit;
  Vec4D clumom = Vec4D(0.,0.,0.,0.), partmom = Vec4D(0.,0.,0.,0.);
  for (cit=clusters->begin();cit!=clusters->end();) {
    clumom += (*cit)->Momentum();
    if (DecayIt((*cit))) cit++;
    else cit=clusters->erase(cit);
  }
  msg.Tracking()<<"Add "<<p_partlist->size()
		<<" particles to the blob --------------------------------------"<<endl;
  for (Part_Iterator pit=p_partlist->begin();pit!=p_partlist->end();) {
    blob->AddToOutParticles((*pit));
    partmom += (*pit)->Momentum();
    pit = p_partlist->erase(pit);
  }

  if (dabs(blob->CheckMomentumConservation().Abs2())>1.e-9) {
    msg.Tracking()<<"Check this : "
		  <<blob->CheckMomentumConservation()<<", "<<blob->CheckMomentumConservation().Abs2()<<endl
		  <<"   Compare with "<<clumom<<" -> "<<partmom<<" = "<<clumom-partmom<<endl
		  <<(*blob)
		  <<"----------------------------------------------------------"<<endl;
  }
  AnalyseThis(blob);
}




bool Cluster_Decay_Handler::DecayIt(Cluster * cluster)
{
  if (p_decayer->Treat(cluster,p_partlist)) {
    if (cluster->GetLeft())  DecayIt(cluster->GetLeft());
    if (cluster->GetRight()) DecayIt(cluster->GetRight());
    return true;
  }
  p_transformer->TreatSingleCluster(cluster,p_partlist);
  return false;
}


void Cluster_Decay_Handler::InitializeAnalysis()
{
  if (!m_analyse) return;
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
}

void Cluster_Decay_Handler::AnalyseThis(Blob * blob)
{
  if (!m_analyse) return;
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
      m_histograms[string("x_p_Pseudoscalars")]->Insert(2.*Vec3D(part->Momentum()).Abs()/91.2);
      break;
    case 211:
      Npiplus++;
      m_histograms[string("x_p_Pseudoscalars")]->Insert(2.*Vec3D(part->Momentum()).Abs()/91.2);
      break;
    case -211:
      Npiminus++;
      m_histograms[string("x_p_Pseudoscalars")]->Insert(2.*Vec3D(part->Momentum()).Abs()/91.2);
      break;
    case 221:
      Neta++;
      m_histograms[string("x_p_Pseudoscalars")]->Insert(2.*Vec3D(part->Momentum()).Abs()/91.2);
      break;
    case 311:
      NK0++;
      m_histograms[string("x_p_Pseudoscalars")]->Insert(2.*Vec3D(part->Momentum()).Abs()/91.2);
      break;
    case -311:
      NK0b++;
      m_histograms[string("x_p_Pseudoscalars")]->Insert(2.*Vec3D(part->Momentum()).Abs()/91.2);
      break;
    case 321:
      NKplus++;
      m_histograms[string("x_p_Pseudoscalars")]->Insert(2.*Vec3D(part->Momentum()).Abs()/91.2);
      break;
    case -321:
      NKminus++;
      m_histograms[string("x_p_Pseudoscalars")]->Insert(2.*Vec3D(part->Momentum()).Abs()/91.2);
      break;
    case 331:
      Netaprime++;
      m_histograms[string("x_p_Pseudoscalars")]->Insert(2.*Vec3D(part->Momentum()).Abs()/91.2);
      break;
    case 113:
      Nrho0++;
      m_histograms[string("x_p_Vectors")]->Insert(2.*Vec3D(part->Momentum()).Abs()/91.2);
      break;
    case 213:
      Nrhoplus++;
      m_histograms[string("x_p_Vectors")]->Insert(2.*Vec3D(part->Momentum()).Abs()/91.2);
      break;
    case -213:
      Nrhominus++;
      m_histograms[string("x_p_Vectors")]->Insert(2.*Vec3D(part->Momentum()).Abs()/91.2);
      break;
    case 223:
      Nomega++;
      m_histograms[string("x_p_Vectors")]->Insert(2.*Vec3D(part->Momentum()).Abs()/91.2);
      break;
    case 313:
      NKstar0++;
      m_histograms[string("x_p_Vectors")]->Insert(2.*Vec3D(part->Momentum()).Abs()/91.2);
      break;
    case -313:
      NKstar0b++;
      m_histograms[string("x_p_Vectors")]->Insert(2.*Vec3D(part->Momentum()).Abs()/91.2);
      break;
    case 323:
      NKstarplus++;
      m_histograms[string("x_p_Vectors")]->Insert(2.*Vec3D(part->Momentum()).Abs()/91.2);
      break;
    case -323:
      NKstarminus++;
      m_histograms[string("x_p_Vectors")]->Insert(2.*Vec3D(part->Momentum()).Abs()/91.2);
      break;
    case 333:
      Nphi++;
      m_histograms[string("x_p_Vectors")]->Insert(2.*Vec3D(part->Momentum()).Abs()/91.2);
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

#include "Primitive_Detector.H"
#include "Primitive_Analysis.H"
#include "Particle_List.H"
#include "MathTools.H"

#ifdef PROFILE__Analysis_Phase
#include "prof.hh"
#else 
#define PROFILE_HERE {}
#define PROFILE_LOCAL(LOCALNAME) {}
#endif

using namespace ANALYSIS;
using namespace ATOOLS;

Primitive_Detector::Primitive_Detector() : 
  p_et_cal(NULL), p_costheta(NULL), p_sintheta(NULL), p_cosphi(NULL), p_sinphi(NULL)
{
  m_splitt_flag=false;

  m_mode=0;
  m_name="SimpleConeJet";
  m_dR=0.5;
  m_etjet_cut=8.;

  m_mineta=-2.5;
  m_maxeta= 2.5;

  m_neta=50;
  m_nphi=40;

  m_delta_eta=(m_maxeta-m_mineta)/double(m_neta);
  m_delta_phi=2.*M_PI/double(m_nphi);

  Init();
}

Primitive_Detector::Primitive_Detector(const double dR,const double etjet_cut,const int mode,const std::string & name): 
  p_et_cal(NULL), p_costheta(NULL), p_sintheta(NULL), p_cosphi(NULL), p_sinphi(NULL)
{
  m_splitt_flag=false;

  m_mode=mode;
  m_name=name;
  m_dR=dR;
  m_etjet_cut=etjet_cut;

  m_mineta=-2.5;
  m_maxeta= 2.5;

  m_neta=40;
  m_nphi=50;

  m_delta_eta=(m_maxeta-m_mineta)/double(m_neta);
  m_delta_phi=2.*M_PI/double(m_nphi);

  Init();
}

Primitive_Detector::~Primitive_Detector() 
{
  ClearAll();
}

void Primitive_Detector::ClearAll()
{
  if (p_et_cal) {
    for (int i=0; i<m_neta;++i) delete [] p_et_cal[i];
    p_et_cal=NULL;
  }
  if (p_jetno) {
    for (int i=0; i<m_neta;++i) delete [] p_jetno[i];
    p_jetno=NULL;
  }
  if (p_costheta) delete [] p_costheta; p_costheta=0;
  if (p_sintheta) delete [] p_sintheta; p_sintheta=0;
  if (p_cosphi) delete [] p_cosphi; p_cosphi=0;
  if (p_sinphi) delete [] p_sinphi; p_sinphi=0;
}



void Primitive_Detector::Init()
{
  PROFILE_HERE;
  if (p_et_cal==NULL) {
    p_et_cal = new double*[m_neta];
    for (int i=0; i<m_neta;++i) p_et_cal[i]=new double[m_nphi];
    p_jetno = new int*[m_neta];
    for (int i=0; i<m_neta;++i) p_jetno[i]=new int[m_nphi];

    p_costheta = new double[m_neta];
    p_sintheta = new double[m_neta];

    p_cosphi = new double[m_nphi]; 
    p_sinphi = new double[m_nphi]; 

    for (int j=0; j<m_nphi; ++j) {
      double phi=m_delta_phi*(j+0.5);
      p_cosphi[j]=cos(phi);
      p_sinphi[j]=sin(phi);
    }

    for (int i=0; i<m_neta; ++i) {
      double eta=m_mineta +m_delta_eta*(i+0.5);
      double theta=2. *atan(exp(-eta));
      p_costheta[i]=cos(theta);
      p_sintheta[i]=sin(theta);
    }
  }

  for (int i=0; i<m_neta; ++i) {
    for (int j=0; j<m_nphi; ++j) {
      p_et_cal[i][j]=0.;
      p_jetno[i][j]=0;
    }
  }

  m_jets.clear();
}

double Primitive_Detector::PseudoRapidityNAzimuthalAngle(const Vec4D & p, double & phi)
{
  PROFILE_HERE;
  double pt2=sqr(p[1])+sqr(p[2]);
  double pp =sqrt(pt2+sqr(p[3]));
  double pz =dabs(p[3]);
  double sn =p[3]/pz;
  if (pt2<1.e-10*pp*pp) {
    phi=0.;
    return sn*20.;
  }
  if (dabs(p[1])<1.e-10*dabs(p[2])) {
    if (p[2]>0) phi=0.5*M_PI;
    else phi =1.5*M_PI;
  }
  else {
    phi=atan2(p[2],p[1]);
    if (phi<0) phi+=2.*M_PI;
  }
  return sn*0.5*log(sqr(pp+pz)/pt2);
}

void Primitive_Detector::SmearEnergy(const Flavour & fl, double & E)
{
  // perform gauss smearing if necessary
}

void Primitive_Detector::FillCalorimeter(const Particle_List & pl)
{
  PROFILE_HERE;
  for (Particle_List::const_iterator it=pl.begin(); it!=pl.end();++it) {
    //    if (!((*it)->Flav().IsLepton() && (*it)->Flav().IntCharge()==0) ) { //ignore neutrino
    if (!((*it)->Flav().IsLepton()) ) { //ignore all leptons
      double phi=0;
      double y=PseudoRapidityNAzimuthalAngle((*it)->Momentum(),phi);
      
      double E=(*it)->Momentum()[0];
      SmearEnergy((*it)->Flav(),E);

      int i=int((y-m_mineta)/m_delta_eta);
      int j=int(phi/m_delta_phi);
      if (i>=0 && i<m_neta) 
	p_et_cal[i][j]+=p_sintheta[i]*E;
    }
  }  
}

void Primitive_Detector::CalcJets() 
{
  PROFILE_HERE;
  // a simplied version of the UA1  jet algorithm with jetradius and 
  //   minimum scalar transverse energy
  //        (RJET=1., EJCUT=5. FOR UA1)
  //  getjet.f obtained by MLM Feb 2004
  //    "Calorimeter simulation obtained from Frank Paige 23 March 1988"
  //    (RESEM=.15, RESHAD=.35 FOR URANIUM CALORIMETER)
  //    (RESEM=.15, RESHAD=.70 FOR IRON CALORIMETER)
  //    (RESEM=.11, RESHAD=.70 FOR CDF)


  //
  double dR = 0.7;
  double etjet_cut=8.;
   
  double dR2=sqr(dR);
  double etstop=1.5;
  //  double eccut=0.1;

  int dneta=int(m_neta*dR/(m_maxeta-m_mineta));
  int dnphi=int(m_nphi*dR/(2.*M_PI));

  for (;;) {
  
    // find highest tower
    double maxet=0;
    int ii=-1,jj=-1;
    for (int i=0; i<m_neta; ++i) {
      for (int j=0; j<m_nphi; ++j) {
	if (p_jetno[i][j]>0) continue;
	if (p_et_cal[i][j]<maxet) continue;
	maxet=p_et_cal[i][j];
	ii=i;
	jj=j;
      }
    }

    if (ii==-1) break;
    if (maxet<etstop) break;


    // add jet:
    Vec4D   jetmom;
    double  jetet=0.;

    for (int i=ii-dneta; i<=ii+dneta; ++i) {
      if (i<0) i=0;
      if (i>=m_neta) break; 
      for (int jp=jj-dnphi; jp<=jj+dnphi; ++jp) {
	int j=jp;
	if (j<0) j+=m_nphi;
	else if (j>=m_nphi) j-=m_nphi;
	double dr2=sqr(m_delta_phi*(jp-jj))+sqr(m_delta_eta*(i-ii));
	if (dr2>dR2) continue;
	if (p_jetno[i][j]>0) continue;
	
	p_jetno[i][j]=m_jets.size()+1;
	// add to jet
	double pt=p_et_cal[i][j];
	double px=pt/p_sintheta[i];
	jetmom[0]+=px;
	jetmom[1]+=pt*p_cosphi[j];
	jetmom[2]+=pt*p_sinphi[j];
	jetmom[3]+=px*p_costheta[i];
	jetet+=pt;
      }
    }
    
    if (jetet<etjet_cut) break;
    m_jets.push_back(Jet_Data(ii,jj,jetmom,jetet));

  }
}


void Primitive_Detector::Evaluate(const Blob_List &,double value,int ncount) 
{
  PROFILE_HERE;
  Particle_List * pl=p_ana->GetParticleList("FinalState");
  Evaluate(*pl,value,ncount);
}

void Primitive_Detector::Evaluate(const Particle_List & pl,double value, int ncount) 
{
  PROFILE_HERE;
  // create jetlist
  Particle_List * pl_jets = p_ana->GetParticleList(m_name);
  if (pl_jets==0) {
    Init();
    FillCalorimeter(pl);
    CalcJets();
    //    Print();

    pl_jets = new Particle_List();
    int i=1;
    for (std::vector<Jet_Data>::iterator it=m_jets.begin();it!=m_jets.end();++it,++i) {
      pl_jets->push_back(new Particle(i,Flavour(kf::jet),it->mom));
    }

    p_ana->AddParticleList(m_name,pl_jets);
  }
}

void Primitive_Detector::Print(std::ostream & s)
{
  PROFILE_HERE;
  s<<" Primitive_Detector "<<std::endl;
  s<<" neta="<<m_neta<<" ("<<m_mineta<<","<<m_maxeta<<")  nphi="<<m_nphi<<std::endl;
  s<<" deta="<<m_delta_eta<<"       dphi="<<m_delta_phi<<std::endl;

  if (p_et_cal && p_jetno) {
    double maxet=0.;
    for (int i=0; i<m_neta;++i) 
      for (int j=0; j<m_nphi;++j)
	if (p_et_cal[i][j]>maxet) maxet=p_et_cal[i][j];
    if (maxet==0.) {
      s<<" --- no entries in detector!!! --- "<<std::endl;
      return;
    }
    maxet*=.11;
    for (int i=0; i<m_neta;++i) {
      for (int j=0; j<m_nphi;++j) {
	int et=int(p_et_cal[i][j]/maxet);
	if (p_et_cal[i][j]==0.)
	  s<<".";
	else
	  s<<et;
      }
      s<<"          ";
      for (int j=0; j<m_nphi;++j) {
	if (p_jetno[i][j]<1)
	  s<<".";
	else if (p_jetno[i][j]<10)
	  s<<p_jetno[i][j];
	else 
	  s<<"#";
      }
      s<<std::endl;
    }
  }
  
}

Primitive_Observable_Base * Primitive_Detector::Copy() const 
{
  //  std::cerr<<" ERROR in Primitive_Detector::Copy() : not splittable "<<std::endl;
  return new Primitive_Detector(m_dR,m_etjet_cut,m_mode,m_name);
}

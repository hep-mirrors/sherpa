#include "AddOns/Analysis/Triggers/Primitive_Calorimeter.H"
#include "ATOOLS/Phys/Particle_List.H"
#include "AddOns/Analysis/Tools/Particle_Qualifier.H"
#include "ATOOLS/Math/MathTools.H"

using namespace ANALYSIS;
using namespace ATOOLS;

Primitive_Calorimeter::Primitive_Calorimeter(const double mineta,const double maxeta,
					     const int neta,const int nphi,
					     const std::string qualifier) :
  Primitive_Detector_Element(neta,nphi,std::string("Hadronic Calorimeter")),
  m_mineta(mineta), m_maxeta(maxeta),
  p_qualifier(NULL)
{
  p_qualifier = ATOOLS::Particle_Qualifier_Getter::GetObject(qualifier,qualifier);
  m_delta_eta     = (m_maxeta-m_mineta)/double(m_nx);
  m_delta_phi     = 2.*M_PI/double(m_ny);

  p_costheta      = new double[m_nx];
  p_sintheta      = new double[m_nx];  
  for (int i=0; i<m_nx; ++i) {
    double eta    = m_mineta +m_delta_eta*(i+0.5);
    double theta  = 2.*atan(exp(-eta));
    p_costheta[i] = cos(theta);
    p_sintheta[i] = sin(theta);
  }
  p_cosphi        = new double[m_ny]; 
  p_sinphi        = new double[m_ny]; 
  for (int j=0; j<m_ny; ++j) {
    double phi    = m_delta_phi*(j+0.5);
    p_cosphi[j]   = cos(phi);
    p_sinphi[j]   = sin(phi);
  }
}

Primitive_Calorimeter::~Primitive_Calorimeter() 
{
  if (p_costheta) { delete [] p_costheta; p_costheta=0; }
  if (p_sintheta) { delete [] p_sintheta; p_sintheta=0; }
  if (p_cosphi)   { delete [] p_cosphi;   p_cosphi=0; }
  if (p_sinphi)   { delete [] p_sinphi;   p_sinphi=0; }
}


double Primitive_Calorimeter::PseudoRapidityNAzimuthalAngle(const Vec4D & p, double & phi)
{
  double pt2 = sqr(p[1])+sqr(p[2]);
  double pp  = sqrt(pt2+sqr(p[3]));
  double pz  = dabs(p[3]);
  double sn  = p[3]/pz;
  if (pt2<1.e-10*pp*pp) { phi = 0.; return sn*20.; }
  if (dabs(p[1])<1.e-10*dabs(p[2])) {
    if (p[2]>0) phi = 0.5*M_PI;
           else phi = 1.5*M_PI;
  }
  else {
    phi = atan2(p[1],p[2]);
    if (phi<0) phi+=2.*M_PI;
  }
  /*
    std::cout<<"Check : "<<p<<" : "<<std::endl
    <<sn*0.5*log(sqr(pp+pz)/pt2)<<"  "
    <<pt2<<" "<<pz<<" ... "<<atan2(sqrt(pt2),pz)<<"->"
    <<tan(atan2(sqrt(pt2),pz))<<"->"
    <<log(tan(atan2(sqrt(pt2),pz)/2.))<<std::endl;
  */
  return sn*0.5*log(sqr(pp+pz)/pt2);
}

void Primitive_Calorimeter::SmearEnergy(const Flavour & fl, double & E)
{
  // perform gauss smearing if necessary
}

void Primitive_Calorimeter::Fill(const Particle_List * pl)
{
  Reset();
  double maxet=0.;
  Vec4D  mom=Vec4D(0.,0.,0.,0.);
  int ii,jj;
  for (Particle_List::const_iterator it=pl->begin(); it!=pl->end();++it) {
    if (p_qualifier!=NULL && !(*p_qualifier)(*it)) continue;
    if (!((*it)->Flav().IsLepton())) {
      double phi = 0;
      double y   = PseudoRapidityNAzimuthalAngle((*it)->Momentum(),phi);
      
      double E   = (*it)->Momentum()[0];
      SmearEnergy((*it)->Flav(),E);

      int i = int((y-m_mineta)/m_delta_eta);
      int j = int(phi/m_delta_phi);
      if (i>=0&&i<m_nx) {
	p_cells[i][j] += p_sintheta[i]*E;
	if (maxet<p_sintheta[i]*E) { maxet = p_sintheta[i]*E; mom=(*it)->Momentum(); ii=i; jj=j; }
      }
    }
  }  
}

void Primitive_Calorimeter::Extract(Particle_List * pl)
{
  for (int i=0;i<m_nx;++i) {
    for (int j=0;j<m_ny;++j) {
      if (p_cells[i][j]!=0.0) {
	Vec4D mom(1.0,p_sintheta[i]*p_cosphi[j],p_sintheta[i]*p_sinphi[j],p_costheta[i]);
	mom*=p_cells[i][j]/dabs(p_sintheta[i]);
	Particle *track = new Particle(1,kf_jet,mom);
	pl->push_back(track);
      }
    }
  }
}

Primitive_Detector_Element * Primitive_Calorimeter::Copy() const 
{
  return new Primitive_Calorimeter(m_mineta,m_maxeta,m_nx,m_ny);
}


void Primitive_Calorimeter::Print(std::ostream & s)
{
  s<<" Primitive_Calorimeter "<<std::endl;
  s<<" neta="<<m_nx<<" ("<<m_mineta<<","<<m_maxeta<<")  nphi="<<m_ny<<std::endl;
  s<<" deta="<<m_delta_eta<<"       dphi="<<m_delta_phi<<std::endl;

  if (p_cells) {
    double maxet = 0.;
    for (int i=0;i<m_nx;++i) { 
      for (int j=0;j<m_ny;++j) {
	if (p_cells[i][j]>maxet) maxet=p_cells[i][j];
      }
    }
    if (maxet==0.) {
      s<<" --- no entries in detector!!! --- "<<std::endl;
      return;
    }
    else {
      for (int i=0;i<m_nx;++i) { 
	for (int j=0;j<m_ny;++j) {
	  if (p_cells[i][j]>0.) s<<i<<" "<<j<<" : "<<p_cells[i][j]<<std::endl;
	}
      }
    }
  }
}

void Primitive_Calorimeter::Reset()
{
  for (int i=0; i<m_nx; ++i) {
    for (int j=0; j<m_ny; ++j) {
      p_cells[i][j]=0.;
    }
  }
}


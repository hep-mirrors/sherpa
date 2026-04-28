#include "BEAM/Ions/Ion_Base.H"
#include "BEAM/Ions/Ion_Propagation.H"
#include "BEAM/Ions/Ion_Parameters.H"
#include "ATOOLS/Org/Message.H"

using namespace BEAM;
using namespace ATOOLS;
using namespace std;

Ion_Base::Ion_Base(beamspectrum type,const ATOOLS::Flavour & flav,
		   const double energy,const double polarisation,
		   const int dir, int mode) :
  Beam_Base(type,flav,energy,polarisation,dir,mode)
{
  ionpars = new Ion_Parameters();
  ionpars->Initialise();
  m_R      = m_beam.Radius();
  m_Rform  = ionpars->GetInt("R_Form");
  m_Pform  = ionpars->GetInt("P_Form");
  m_kFermi = ionpars->Get("K_Fermi");
  m_rho0   = ionpars->Get("Rho_0");
  m_rmin2  = ionpars->Get("R_min^2(Pauli)");
  m_dmin2  = 0.1;
  msg_Info()<<METHOD<<"("<<flav<<", A = "<<m_A<<", Z = "<<m_Z<<"):\n"
	    <<"R               = "<<m_R<<" fm,\n"
	    <<"K_Fermi         = "<<m_kFermi<<" GeV,\n"
	    <<"rho_0           = "<<m_rho0<<" fm^{-3},\n"
	    <<"r_min(Blocking) = "<<sqrt(m_rmin2)<<" fm.\n";
  Initialise();
}

void Ion_Base::Initialise() {
  long int comb = m_beam.Kfcode()-1000000000;
  m_Np   = m_Z = int(comb/10000.);
  m_A    = int((comb-m_Z*10000)/10.);
  m_Nn   = m_A-m_Z;
  p_prop = new Ion_Propagation(this);

  for (size_t i=0;i<m_Np;i++)
    m_constituents.push_back(new Ion_Constituent(Flavour(kf_p_plus)));
  for (size_t i=0;i<m_Nn;i++)
    m_constituents.push_back(new Ion_Constituent(Flavour(kf_n)));
  MakeInitialPositions(true);
  MakeInitialMomenta();
  PrintConstituents();
  p_prop->Evolve();
  exit(1);
}

void Ion_Base::MakeInitialPositions(bool withBlocking) {
  Vec4D summedpos = (m_Rform==0 ?
		     MakeInitialPositionsHardSphere(withBlocking) :
		     MakeInitialPositionsGaussian(withBlocking) );
  /////////////////////////////////////////////////////////////////////
  // Collective shift to center the nucelons such that the sum of
  // their positions is at the origin (which, assuming degeneracy of
  // neutron and proton mass is also the centre-of mass, and hence the
  // position of the ion).
  /////////////////////////////////////////////////////////////////////
  summedpos = summedpos/double(m_A);
  double r2 = 0;
  for (size_t i=0;i<m_A;i++) {
    m_constituents[i]->m_r -= summedpos;
    r2 += Vec3D(m_constituents[i]->m_r).Sqr();
  }
  double meanR2 = r2/double(m_A);
  double scale  = m_R/sqrt(meanR2);
  for (size_t i=0;i<m_A;i++)
    m_constituents[i]->m_r *= scale;
  CheckRadius();
}

Vec4D Ion_Base::MakeInitialPositionsHardSphere(bool withBlocking) {
  /////////////////////////////////////////////////////////////////////
  // Assume that the "radius per nucleon" is inverse to the nucleon 
  // density rho_0 in the nucleus, then use the scaling of the radius 
  // with the number of nucleons A^(1/3) - however, as the nucleons are 
  // smeared out and we only place their centres inside the hard sphere 
  // we have to subtract "one half of a nucleon layer" - this is the
  // weird-looking term [A^(1/3)-1]^3.
  /////////////////////////////////////////////////////////////////////
  double radius =  ( pow(3./(4.*M_PI*m_rho0),1./3.) * 
		     pow(1./2.*(m_A+pow(pow(m_A,1./3.)-1.,3.)),1./3.) );
  Vec4D  sum    = Vec4D(0.,0.,0.,0.), pos;
  for (size_t i=0;i<m_A;i++) {
    do {
      ////////////////////////////////////////////////////////////////
      // Assuming a distribution in r given by dr^3 = r^2dr d^2Omega
      // Rescaling the <r^2> of the hard sphere to the measured
      // radius through a factor of 5/3.
      ////////////////////////////////////////////////////////////////
      double r      = (m_A==2 ?
		       m_R :
		       sqrt(5./3.)*pow(ran->Get(),1./3.)*radius);
      double ctheta = 1.-2.*ran->Get(), stheta = sqrt(1.-ctheta*ctheta);
      double phi    = 2.*M_PI*ran->Get();
      pos           = r*Vec4D(0.,stheta*cos(phi),stheta*sin(phi),ctheta);
    } while (withBlocking && IsRBlocked(m_constituents[i]->m_flav,pos));
    sum += m_constituents[i]->m_r = pos;
  }
  return sum;
}

Vec4D Ion_Base::MakeInitialPositionsGaussian(bool withBlocking) {
  Vec4D  sum = Vec4D(0.,0.,0.,0.), pos;
  for (size_t i=0;i<m_A;i++) {
    do {
      double r      = (m_A==2 ?
		       m_R :
		       m_R*dabs(ran->GetGaussian()));
      double ctheta = 1.-2.*ran->Get(), stheta = sqrt(1.-ctheta*ctheta);
      double phi    = 2.*M_PI*ran->Get();
      pos           = r*Vec4D(0.,stheta*cos(phi),stheta*sin(phi),ctheta);
    } while (withBlocking && IsRBlocked(m_constituents[i]->m_flav,pos));
    sum += m_constituents[i]->m_r = pos;
  }
  return sum;
}

void Ion_Base::CheckRadius() {
  double R2 = 0.;
  for (size_t i=0;i<m_A;i++) 
    R2  += Vec3D(m_constituents[i]->m_r).Sqr();
  msg_Out()<<METHOD<<": <R> = "<<sqrt(R2)/double(m_A)<<" "
	   <<"vs. R = "<<m_R<<".\n";
}

bool Ion_Base::IsRBlocked(const Flavour & flav,const Vec4D & pos) {
  /////////////////////////////////////////////////////////////////////
  // Making sure that nucleons of the same kind do not "overlap", i.e.
  // have a minimal distance.
  /////////////////////////////////////////////////////////////////////
  for (size_t i=0;i<m_A;i++) {
    if (m_constituents[i]->m_flav==flav &&
	Vec3D(pos-m_constituents[i]->m_r).Sqr() < m_rmin2) {
      return true;
    }
  }
  return false;
}

void Ion_Base::MakeInitialMomenta(bool withBlocking) {
  withBlocking = false;
  Vec4D sum = Vec4D(0.,0.,0.,0.);
  if (m_Pform==1) { p_prop->Update(); }
  for (size_t i=0;i<m_A;i++) {
    msg_Out()<<METHOD<<"["<<i<<"]: "<<p_prop->Pmax(i)<<"\n";
    Vec4D  mom;
    do {
      double p      = (m_Pform==0 ?
		       pow(ran->Get(),1./3.)*m_kFermi :
		       ran->Get()*p_prop->Pmax(i));
      double ctheta = 1.-2.*ran->Get(), stheta = sqrt(1.-ctheta*ctheta);
      double phi    = 2.*M_PI*ran->Get();
      mom           = p*Vec4D(0.,stheta*cos(phi),stheta*sin(phi),ctheta);
    } while (withBlocking && IsPRBlocked(i,m_constituents[i]->m_flav,
					 m_constituents[i]->m_r,mom));
    sum          += m_constituents[i]->m_p = mom;
  }
  sum = sum/double(m_A);
  for (size_t i=0;i<m_A;i++) { m_constituents[i]->m_p -= sum; }
  FixEnergies();
}

bool Ion_Base::IsPRBlocked(const size_t ref,const ATOOLS::Flavour & flav,
			   const ATOOLS::Vec4D & pos,const ATOOLS::Vec4D & mom) {
  for (size_t i=0;i<m_A;i++) {
    if (i==ref) continue;
    if (m_constituents[i]->m_flav==flav &&
	(Vec3D(pos-m_constituents[i]->m_r).Sqr() *
	 Vec3D(mom-m_constituents[i]->m_p).Sqr() ) < m_dmin2) {
      msg_Out()<<"blocking: "
	       <<Vec3D(pos-m_constituents[i]->m_r).Sqr()<<" * "
	       <<Vec3D(mom-m_constituents[i]->m_p).Sqr()<<" = "
	       <<(Vec3D(pos-m_constituents[i]->m_r).Sqr() *
		  Vec3D(mom-m_constituents[i]->m_p).Sqr())<<" < "<<m_dmin2<<"\n";
      return true;
    }
  }
  return false;
}


void Ion_Base::FixEnergies() {
  ///////////////////////////////////////////////////////////////////
  //
  // Logic is that we define a global scale factor kappa such that
  // E = M = sum_i sqrt(p_i^2+kappa m_i^2)
  // with M the mass of the ion, p_i the three-momenta of the
  // constituents i and m_i their masses.
  // 
  ///////////////////////////////////////////////////////////////////
  double kappa = 1, f = kernel(kappa), M = m_beam.Mass();
  while (dabs(f)>M*1.e-6) {
    kappa -= f/dkernel(kappa);
    f      = kernel(kappa);
  }
  for (size_t i=0;i<m_A;i++) {
    Ion_Constituent * ic = m_constituents[i];
    double m   = ic->m_flav.Mass();
    ic->m_p[0] = sqrt(Vec3D(ic->m_p).Sqr()+kappa*sqr(m));
    ic->m_E    = ic->m_p[0]-m;
    ic->m_q2   = ic->m_p.Abs2()-sqr(m);
  }
}

void Ion_Base::PrintConstituents() {
  msg_Out()<<METHOD<<" for "<<m_beam<<", P = "<<m_lab<<", R = "<<m_position<<":\n";
  Vec4D R = Vec4D(0.,0.,0.,0.), P = Vec4D(0.,0.,0.,0.);
  for (size_t i=0;i<m_A;i++) {
    Ion_Constituent * ic = m_constituents[i];
    msg_Out()<<"* "<<std::setw(4)<<ic->m_flav<<": "
	     <<"E = "<<std::setprecision(4)<<std::setw(5)<<ic->m_E<<", "
	     <<"q2 = "<<std::setprecision(4)<<std::setw(5)<<ic->m_q2<<", "
	     <<"r = "<<ic->m_r<<", "
	     <<"p = "<<ic->m_p<<"\n";
    R += ic->m_r;
    P += ic->m_p; 
  }
  msg_Out()<<"*** Check the sums: P = "<<P<<", R = "<<R<<".\n";
}

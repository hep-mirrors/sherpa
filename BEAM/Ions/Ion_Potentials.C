#include "BEAM/Ions/Ion_Potentials.H"
#include "BEAM/Ions/Ion_Parameters.H"
#include "ATOOLS/Org/Message.H"

using namespace BEAM;
using namespace ATOOLS;
using namespace std;

Ion_Potentials::Ion_Potentials(size_t A) : m_A(A) {
  // Units of the constants below:
  // [m_rho0]  = 1/fm^3, [m_sigma] = fm,  [m_sigma2] = fm^2
  // [m_V2loc] = GeV,    [m_V3loc] = GeV, [m_nu] = 1
  // [m_aYuk]  = fm,     [m_Vyuk]  = GeV/fm^3
  m_rho0     = ionpars->Get("Rho_0");
  m_sigma    = ionpars->Get("Sigma_Psi");
  m_sigma2   = sqr(m_sigma);
  m_norm     = pow(4.*m_sigma2,-3./2.);
  m_rhonorm  = pow(2.*m_sigma2,-3./2.);  
  m_V2loc    = ionpars->Get("V^(loc)_2");
  m_V3loc    = ionpars->Get("V^(loc)_3");
  m_maxRloc  = ionpars->Get("R^(loc)_max");
  if (m_maxRloc<0.) m_maxRloc = 1.e99;
  m_nu       = ionpars->Get("nu^(loc)");
  m_U2pref   = m_V2loc;
  m_U3pref   = m_V3loc/(pow(m_nu+1.,3./2.)*pow(2.*M_PI*m_sigma2, 3.*m_nu/2.));
  m_aYuk     = ionpars->Get("a^(Yuk)");
  m_VYuk     = ionpars->Get("V^(Yuk)");
  m_maxRYuk  = ionpars->Get("R^(Yuk)_max");
  if (m_maxRYuk<0.) m_maxRYuk = 1.e99;
  m_VYukpref = m_VYuk*exp(m_sigma2/sqr(m_aYuk));
  m_UYukpref = m_VYuk*exp(m_sigma2/sqr(m_aYuk))/sqrt(2.*M_PI*m_sigma2);

  m_tilderho.resize(m_A);
  m_rho.resize(m_A);
  m_U.resize(m_A);
}

void Ion_Potentials::UpdateDensities() {
  //////////////////////////////////////////////////////////////////////
  // rho is the nucleon/baryon density, which comes from the product of
  // Gaussian-smeared single-nucleon wave functions. 
  // tilderho is the interaction density
  //////////////////////////////////////////////////////////////////////
  for (size_t i=0;i<m_A;i++) {
    double tilderho = 0., rho = 0.;
    for (size_t j=0;j<m_A;j++) {
      if (j==i) continue;
      tilderho += exp(   (*p_deltar2norm)[0][i][j]);
      rho      += exp(2.*(*p_deltar2norm)[0][i][j]);
    }
    m_tilderho[i] = m_norm*tilderho;
    m_rho[i]      = m_rhonorm*rho;
  }
}

void Ion_Potentials::UpdatePotentials() {
  for (size_t i=0;i<m_A;i++) {
    m_U[i] = LocalPotential2(i)+LocalPotential3(i)+YukawaPotential(i);
    msg_Out()<<"U("<<i<<") = "
	     <<LocalPotential2(i)<<" + "<<LocalPotential3(i)<<" + "
	     <<YukawaPotential(i)<<" = "<<m_U[i]<<" from rho = "<<m_tilderho[i]<<"\n";
    //msg_Out()<<"U("<<i<<") = "
    //	     <<(-0.38*m_rho[i]/0.17)<<" + "<<(0.303*sqr(m_rho[i]/0.17))<<" = "
    //	     <<(-0.38*m_rho[i]/0.17+0.303*sqr(m_rho[i]/0.17))
    //	     <<" from rho = "<<m_rho[i]<<"\n\n";
  }
}

const double Ion_Potentials::LocalPotential2(const size_t & i) {
  double U = 0.;
  if (dabs(m_U2pref)>1.e-10) {
    for (size_t j=0;j<m_A;j++) {
      if (i==j || (*p_deltar)[0][i][j]>m_maxRloc) continue;
      U += m_U2pref * m_norm * exp(-(*p_deltar2norm)[0][i][j]);
    }
  }
  return U;
}

ATOOLS::Vec3D
Ion_Potentials::NablaLocalPotential2(const size_t & level,const size_t & i) {
  Vec3D nabla(0.,0.,0.);
  if (dabs(m_U2pref)>1.e-10) {
    for (size_t j=0;j<m_A;j++) {
      if (i==j || (*p_deltar)[0][i][j]>m_maxRloc) continue;
      nabla -= ( m_U2pref * 
		 (*p_deltarvec)[level][i][j]/(2.*m_sigma2) *
		 m_norm * exp(-(*p_deltar2norm)[level][i][j]) );
    }
  }
  return nabla;
}

const double Ion_Potentials::LocalPotential3(const size_t & i) {
  double U = 0.;
  if (dabs(m_U3pref)>1.e-10 || m_A<3) {
    for (size_t j=0;j<m_A;j++) {
      for (size_t k=j;k<m_A;k++) {
	if (i==j || i==k || j==k ||
	    (*p_deltar)[0][i][j]>m_maxRloc || (*p_deltar)[0][i][k]>m_maxRloc) continue;
	U += ( m_U3pref *
	       m_norm *  exp(-(*p_deltar2norm)[0][i][j]) *
	       m_norm *  exp(-(*p_deltar2norm)[0][i][k]) );
      }
    }
  }
  return U;
}

ATOOLS::Vec3D
Ion_Potentials::NablaLocalPotential3(const size_t & level,const size_t & i) {
  Vec3D nabla(0.,0.,0.);
  if (dabs(m_U3pref)>1.e-10) {
    for (size_t j=0;j<m_A;j++) {
      for (size_t k=j;k<m_A;k++) {
	if (i==j || i==k || j==k ||
	    (*p_deltar)[0][i][j]>m_maxRloc || (*p_deltar)[0][i][k]>m_maxRloc) continue;
	nabla -= ( m_U3pref *
		   ((*p_deltarvec)[level][i][j]+(*p_deltarvec)[level][i][k])/(3.*m_sigma2) *
		   m_norm * exp(-(*p_deltar2norm)[level][i][j]) *
		   m_norm * exp(-(*p_deltar2norm)[level][i][k]) );
      }
    }
  }
  return nabla;
}

const double Ion_Potentials::YukawaPotential(const size_t & i) {
  double U = 0.;
  if (dabs(m_VYukpref)>1.e-10) {
    for (size_t j=0;j<m_A;j++) {
      if (i==j) continue;
      double rij_by_aY = (*p_deltar)[0][i][j]/m_aYuk;
      if (rij_by_aY>m_maxRYuk/m_aYuk) continue;
      double arg1      = m_sigma/m_aYuk-(*p_deltar)[0][i][j]/(2.*m_sigma);
      double arg2      = m_sigma/m_aYuk+(*p_deltar)[0][i][j]/(2.*m_sigma);
      U += ( m_VYukpref/rij_by_aY *
	     exp(-rij_by_aY) * (1.-erf(arg1)) - exp(+rij_by_aY) * (1.-erf(arg2)) );
    }
  }
  return U;
}

ATOOLS::Vec3D
Ion_Potentials::NablaYukawa(const size_t & level,const size_t & i) {
  Vec3D nabla(0.,0.,0.);
  if (dabs(m_UYukpref)>1.e-10) {
    for (size_t j=0;j<m_A;j++) {
      if (i==j) continue;
      double rij_by_aY = (*p_deltar)[level][i][j]/m_aYuk;
      if (rij_by_aY>m_maxRYuk/m_aYuk) continue;
      double exppref   = 1./(M_PI*m_sigma*rij_by_aY);
      double erfpref1  = 1.+1./rij_by_aY;
      double erfpref2  = 1.-1./rij_by_aY;
      double arg1      = m_sigma/m_aYuk-(*p_deltar)[level][i][j]/(2.*m_sigma);
      double arg2      = m_sigma/m_aYuk+(*p_deltar)[level][i][j]/(2.*m_sigma);
      nabla += ( m_UYukpref * 
		 ((*p_deltarvec)[level][i][j]/(*p_deltar)[level][i][j]) * 
		 (exp(-rij_by_aY)*( exppref*exp(-sqr(arg1))-erfpref1*(1.-erf(arg1))) -
		  exp(+rij_by_aY)*( exppref*exp(-sqr(arg2))+erfpref2*(1.-erf(arg2))) ) );
    }
  }
  return nabla;
}

ATOOLS::Vec3D Ion_Potentials::Nabla(const size_t & level,const size_t & i) {
  return ( NablaLocalPotential2(level,i) +
	   NablaLocalPotential3(level,i) +  
	   NablaYukawa(level,i) );
}

#include <iomanip>
#include <array>
#include <vector>
#include "ATOOLS/Math/Vector.H"
#include "METOOLS/Currents/C_Spinor.H"
#include "METOOLS/Main/Polarization_Tools.H"
#include "ATOOLS/Org/Message.H"
#include "ATOOLS/Org/Exception.H"
#include "assert.h"

#ifndef SQRT_05
#define SQRT_05 0.70710678118654757
#endif

#ifndef sqrttwo
#define sqrttwo 1.4142135623730950488
#endif

using namespace METOOLS;
using namespace ATOOLS;
using namespace std;
using SpinorBasis = std::vector<Polarization_Spinor>;


Polarization_Vector::Polarization_Vector(Vec4D p, bool anti,
                                         bool out) : std::vector<Vec4C>()
{
  Init(p);
}

Polarization_Vector::Polarization_Vector(Vec4D p, double m2, bool anti,
                                         bool out) : std::vector<Vec4C>()
{
  if(!IsZero((p.Abs2()-m2)/(fabs(p[0])+fabs(m2)), 1e-6)) {
    PRINT_INFO("Undefined for p.Abs2()="<<p.Abs2()<<" and m2="<<m2);
  }
  Init(p);
}

Polarization_Vector::Polarization_Vector(Vec4D p, Vec4D ref_mom) : std::vector<Vec4C>()
{
  Init(p, ref_mom);
  // delete auxiliary polarization vector
  pop_back();
}
Polarization_Spinor_Basis::Polarization_Spinor_Basis(const ATOOLS::Vec4D& p, double M2,
    const ATOOLS::Vec4D& ref_mom, bool anti, bool out) {
  m_s.push_back(Polarization_Spinor(p, M2, +1, ref_mom, anti, out));
  m_s.push_back(Polarization_Spinor(p, M2, -1, ref_mom, anti, out));
  }

  void Polarization_Spinor_Basis::BoostBack(const ATOOLS::Vec4D& beta) {
    for (size_t i(0); i<m_s.size(); ++i) METOOLS::BoostBack_Spinor(m_s[i], beta);
}

Polarization_Spinor::Polarization_Spinor(Vec4D p, double M2, int h, Vec4D ref_mom, bool anti, bool out)
  : m_hel(h), m_p(p), m_ref(ref_mom), m_M2(M2), m_anti(anti), m_out(out)
{
  if (h != 1 && h != -1) {
    THROW(fatal_error, "helicity must be +1 or -1")
  }
  Init(p, ref_mom);     // setzt m_k, m_kp, m_km (die Zeilen davor sind doppelt, kannst du löschen)
  FillComponents();     // schreibt m_u[0..3]
  std::cout << "u = " << (*this)[0] << " " << (*this)[1] << " " << (*this)[2] << " " << (*this)[3] << "\n";
  double nrm(0.0); for (size_t a(0); a<4; ++a) nrm += std::norm((*this)[a]);
  std::cout << "|u|^2 = " << nrm << "  (erwartet 2E = " << 2.0*p[0] << ")\n";
}

Polarization_Spinor::Polarization_Spinor(Vec4D p, int h, Vec4D ref_mom, bool anti, bool out)
  : m_hel(h), m_p(p), m_ref(ref_mom), m_M2(0.0), m_anti(anti), m_out(out)
{ 
  if (h != 1 && h != -1) {
    THROW(fatal_error, "helicity must be +1 or -1")
  }
  Init(p, ref_mom);
  FillComponents();
}

void Polarization_Spinor::FillComponents()
{
  typedef std::complex<double> C;
  const double m = sqrt(std::max(m_M2, 0.0));

  // masseloser Spinor zur Referenz k mit entgegengesetzter Helizität
  METOOLS::CSpinor<double> uk(1, 1, -m_hel, m_k);

  const double pk = m_p * m_k;            // Minkowski-Produkt
  if (pk < 1e-12) {
    THROW(fatal_error, "p and ref_mom are (anti)parallel: p.k = 0")
  }
  const double N = sqrt(2.0 * pk);

  const double p0 = m_p[0], p1 = m_p[1], p2 = m_p[2], p3 = m_p[3];
  const C i(0.0, 1.0);

  // p·sigma = ((p3, p1 - i p2), (p1 + i p2, -p3))
  const C s00 = p3,          s01 = p1 - i*p2;
  const C s10 = p1 + i*p2,   s11 = -p3;

  // (pslash uk): oben = (p0 - p·sigma) * uk[2..3], unten = (p0 + p·sigma) * uk[0..1]
  C t[4];
  t[0] = (p0 - s00)*uk[2] - s01*uk[3];
  t[1] = -s10*uk[2] + (p0 - s11)*uk[3];
  t[2] = (p0 + s00)*uk[0] + s01*uk[1];
  t[3] = s10*uk[0] + (p0 + s11)*uk[1];

  for (size_t a(0); a<4; ++a) m_u[a] = (t[a] + m*uk[a]) / N;
}

void Polarization_Vector::Init(Vec4D p, Vec4D ref_mom)
{
  if (ref_mom.Abs2() > 1e-8) {
    THROW(fatal_error, "new reference vector has to be light-like")
  }
  m_k=ref_mom;
  m_kp=SpinorType(1,m_k);
  m_km=SpinorType(-1,m_k);

  push_back(EMP(p));
  push_back(EMM(p));
  push_back(EML(p));
  push_back(Vec4C(Complex(0.,0.), Complex(0.,0.), Complex(0.,0.),
                  Complex(0.,0.)));
}

Vec4C Polarization_Vector::VT(const SpinorType &a,const SpinorType &b)
{
  Vec4C e;
  e[0]=a.U1()*b.U1()+a.U2()*b.U2();
  e[SpinorType::R3()]=a.U1()*b.U1()-a.U2()*b.U2();
  e[SpinorType::R1()]=a.U1()*b.U2()+a.U2()*b.U1();
  e[SpinorType::R2()]=Complex(0.0,1.0)*(a.U1()*b.U2()-a.U2()*b.U1());
  return e;
} //Das transformiert den vector dann im Endeffekt?

Vec4C Polarization_Vector::EM(const Vec4D &p)
{
  SpinorType pp(1,p);
  Vec4C e(VT(pp,m_km));
  return e/(sqrttwo*std::conj(m_kp*pp));
}

Vec4C Polarization_Vector::EP(const Vec4D &p)
{
  SpinorType pm(-1,p);
  Vec4C e(VT(m_kp,pm));
  return e/(sqrttwo*std::conj(m_km*pm));
}

Vec4C Polarization_Vector::EMM(const Vec4D &p)
{
  return EM(p-p.Abs2()/(2.0*m_k*p)*m_k);
}

Vec4C Polarization_Vector::EMP(const Vec4D &p)
{
  return EP(p-p.Abs2()/(2.0*m_k*p)*m_k);
}

Vec4C Polarization_Vector::EML(const Vec4D &p)
{
  double a(p.Abs2()/(2.0*m_k*p));
  Vec4D b(p-a*m_k);
  SpinorType bm(-1,b), bp(1,b), am(-1,m_k), ap(1,m_k);
  Vec4C e(VT(bp,bm)-a*VT(ap,am));
  return e/(2.0*p.Mass());
}


Polarization_Tensor::Polarization_Tensor(Vec4D p, double m2, bool anti,
                                         bool out) : std::vector<CMatrix>()
{
  Polarization_Vector eps(p,m2);
  
  CMatrix tensor(4);
  for(int mu = 0; mu < 4; mu++){
    for(int nu = 0; nu < 4; nu++){
      Complex eps_munu = (eps[0][mu]*eps[0][nu]);
      tensor[mu][nu] = eps_munu;
    }
  }
  push_back(tensor);

  for(int mu = 0; mu < 4; mu++){
    for(int nu = 0; nu < 4; nu++){
      Complex eps_munu = (eps[0][mu]*eps[2][nu]+eps[2][mu]*
                          eps[0][nu])/sqrt(2.0);
      tensor[mu][nu] = eps_munu;
    }
  }
  push_back(tensor);

  for(int mu = 0; mu < 4; mu++){
    for(int nu = 0; nu < 4; nu++){
      Complex eps_munu = (eps[0][mu]*eps[1][nu]+eps[1][mu]*
                          eps[0][nu]
                          -2.0*eps[2][mu]*eps[2][nu])/sqrt(6.0);
      tensor[mu][nu] = eps_munu;
    }
  }
  push_back(tensor);

  for(int mu = 0; mu < 4; mu++){
    for(int nu = 0; nu < 4; nu++){
      Complex eps_munu = (eps[1][mu]*eps[2][nu]+eps[2][mu]*
                          eps[1][nu])/sqrt(2.0);
      tensor[mu][nu] = eps_munu;
    }
  }
  push_back(tensor);

  for(int mu = 0; mu < 4; mu++){
    for(int nu = 0; nu < 4; nu++){
      Complex eps_munu = (eps[1][mu]*eps[1][nu]);
      tensor[mu][nu] = eps_munu;
    }
  }
  push_back(tensor);
}


double g(int mu, int nu)
{
  if (mu<0 || mu>3 || nu<0 || nu>3)
    std::cout<<"wrong indices in g(mu, nu)."<<std::endl;
  if (mu!=nu) return 0.0;
  if (mu==0) return 1.0;
  if (mu>0 && mu<4) return -1.0;
  return 0.0;
}

double S(int mu, int nu, Vec4D p)
{
  return p[mu]*p[nu]/p.Abs2()-g(mu, nu);
}

void Polarization_Vector::Test(Vec4D p)
{
  bool massless(IsZero(p.Mass(),1e-6));

  std::cout<<METHOD<<": Testing transversality..."<<std::endl;
  bool success=true;
  for(int s=0;s<(massless?2:3);s++) {
    Complex d = (*this)[s]*p;
    if(!IsZero(d)) {
      success=false;
      msg_Out()<<"s="<<s<<std::endl;
      msg_Out()<<"  d="<<d<<" should be: 0.0"<<std::endl;
    }
  }
  std::cout<<METHOD<<"... "<<(success?"passed":"failed")<<std::endl<<std::endl;

  std::cout<<METHOD<<": Testing orthogonality..."<<std::endl;
  success=true;
  for(int s=0;s<(massless?2:3);s++) {
    for(int sp=0;sp<(massless?2:3);sp++) {
      Complex d = conj((*this)[s])*(*this)[sp];
      if( (s==sp && !IsEqual(d,Complex(-1.0,0.0))) ||
          (s!=sp && !IsZero(d))) {
        success=false;
        msg_Out()<<"s="<<s<<" sp="<<sp<<std::endl;
        msg_Out()<<"  d="<<d<<" should be: "<<(s==sp?-1.0:0.0)<<std::endl;
      }
    }
  }
  std::cout<<METHOD<<"... "<<(success?"passed":"failed")<<std::endl<<std::endl;

  std::cout<<METHOD<<": Testing completeness..."<<std::endl;
  success=true;
  for(int mu=0;mu<4;mu++) {
    for(int nu=0;nu<4;nu++) {
      Complex d(0.0,0.0);
      Complex ref;
      if (massless) {
        for(int s=0;s<2;s++) {
          d+=conj((*this)[s])[nu]*(*this)[s][mu];
        }
        Vec4D q(p[0], -p[1], -p[2], -p[3]);
        ref=(p[mu]*q[nu]+p[nu]*q[mu])/(p*q)-g(mu, nu);
      }
      else {
        for(int s=0;s<3;s++) {
          d+=conj((*this)[s])[nu]*(*this)[s][mu];
        }
        ref=p[mu]*p[nu]/p.Abs2()-g(mu, nu);
      }
      if(!IsEqual(d,ref)) {
        success=false;
        msg_Out()<<"  mu="<<mu<<std::endl;
        msg_Out()<<"    nu="<<nu<<std::endl;
        msg_Out()<<"      d="<<d<<std::endl;
        msg_Out()<<"      r="<<ref<<std::endl;
      }
    }
  }
  std::cout<<METHOD<<"... "<<(success?"passed":"failed")<<std::endl<<std::endl;
}

void Polarization_Spinor::BoostBack(const ATOOLS::Vec4D& beta) {
  if (beta == ATOOLS::Vec4D(0.,0.,0.,0.)) return;
  METOOLS::BoostBack_Spinor(*this, beta);   // boostet m_u von *this
  m_beta = beta;                            // m_p und m_ref bleiben CM-Werte
}

void METOOLS::BoostBack_Spinor(METOOLS::CSpinor<double> &s, const ATOOLS::Vec4D &l) 
{
  //TODO: Spinor boost back (spinor is not a 4vector - another representation of Lorentz group)
  // need rapidity
  double m(sqrt(l.Abs2()));
  //if (s.B()<0) {
    //METOOLS::CSpinor<double> k(s.Bar());
    //BoostBack(k);
    //s=k.Bar();
    //return;
  //}
  double n(1.0/sqrt(2.0*m*(l[0]+m)));
  // cosh(eta)
  double c(n*(l[0]+m));
  // sinh * vec(n)
  double sx(n*l[1]), sy(n*l[2]), sz(n*l[3]);
  double sgn(-1.0);
  for (int i(0);i<4;i+=2,sgn=-sgn) {
    Complex a(s[i]), b(s[i+1]);
    s[i] = c*a+sgn*(sz*a+Complex(sx,-sy)*b);
    s[i+1] = c*b+sgn*(Complex(sx,sy)*a-sz*b);
  }

}

std::vector<std::vector<Complex>>
Polarization_Spinor::BasisTrafo_Spinor(const SpinorBasis& old_basis, bool pol_checks) const {
  if (old_basis.size() != 2) {
    THROW(fatal_error, "BasisTrafo_Spinor expects two spinors in old basis")
  }
  
  int m_nhel = old_basis.size();
  Polarization_Spinor partner(m_p,m_M2,-m_hel,m_ref,m_anti,m_out);
  if (m_beta != ATOOLS::Vec4D(0.,0.,0.,0.)){
    METOOLS::BoostBack_Spinor(partner, m_beta);
  }
  std::array<const Polarization_Spinor*,2> new_spinor;
  new_spinor[0] = (m_hel > 0) ? this : &partner;
  new_spinor[1] = (m_hel > 0) ? &partner : this;
  std::cout << "p = " << m_p << ", hel = " << m_hel << ", ref = " << m_ref << "\n";
  std::cout << "m_kp = " << m_kp << "\n m_km = " << m_km << "\n";
  std::cout << "this[0..3] = " << (*this)[0] << ", " << (*this)[1] << ", "
          << (*this)[2] << ", " << (*this)[3] << "\n";
  // erstelle koeffizientenmatrix
  vector<vector<Complex>> coeff(m_nhel, vector<Complex>(m_nhel, Complex(0, 0)));

  // erstelle arrays, um dinge besser zu speichern
  std::array<std::array<Complex,2>, 2> array_old_basis;
  std::array<std::array<Complex,2>, 2> array_old_basis_inv;
  std::array<std::array<Complex,2>, 2> array_new_basis;
  for (size_t i(0); i<2; i++)
  std::cout << "old[" << i << "] = ( "
            << old_basis[i][0] << ", " << old_basis[i][1] << ", "
            << old_basis[i][2] << ", " << old_basis[i][3] << " )\n";
  //determinante
  Complex det = old_basis[0][0] * old_basis[1][1] - old_basis[1][0] * old_basis[0][1];
  if (abs(det.real()) < 1e-8 && abs(det.imag()) < 1e-8){
    std::cout << "determinante is" << det << std::endl;
    THROW(fatal_error, "matrix of polarisation vectors in old spin basis is not invertible")
  }

  //inverse
  array_old_basis_inv[0][0]=(Complex(1, 0) / det)*(old_basis[1][1]);
  array_old_basis_inv[0][1]=(Complex(1, 0) / det)*(-old_basis[0][1]);
  array_old_basis_inv[1][0]=(Complex(1, 0) / det)*(-old_basis[1][0]);
  array_old_basis_inv[1][1]=(Complex(1, 0) / det)*(old_basis[0][0]);

  for (size_t i(0); i<m_nhel; i++){
    for (size_t j(0); j<m_nhel; j++){
      array_new_basis[i][j] = (*new_spinor[i])[j+1];
      array_old_basis[i][j] = old_basis[i][j+1];
    }
  }
  
  // calculate basis transformation matrix: coeff = new_basis * (old_basis)^-1; each row contains the coeff for one
  // polarisation vector in new basis
  Complex tmp(0, 0);
  Complex test(0, 0);
  for (size_t i(0); i<m_nhel; i++){
    Complex zero_comp(0, 0), one_comp(0, 0), two_comp(0,0), three_comp(0,0);
    for (size_t j(0); j<m_nhel; j++){
      tmp = Complex(0, 0);
      test = Complex(0, 0);
      for (size_t k(0); k<m_nhel; k++){
        tmp += (array_new_basis[i][k] * array_old_basis_inv[k][j]);
        test += (array_old_basis_inv[i][k] * array_old_basis[k][j]);
      }

      coeff[i][j] = tmp;

      // test, whether inverse matrix is calculated correctly
      if ((i==j && (abs(test.real()-1)>1e-8 || abs(test.imag())>1e-8)) || (i!=j &&
                                                                           (abs(test.real())>1e-8 || abs(test.imag())>1e-8))){
        THROW(fatal_error, "multiplication of matrix and its inverse does not result in identity matrix")
      }
      // test, whether transformation also works for the zeroth component of the polarisation vectors
      zero_comp += coeff[i][j] * old_basis[j][0];
      // test system of equation's solution
      one_comp += coeff[i][j] * old_basis[j][1];
      two_comp   += coeff[i][j] * old_basis[j][2];
      three_comp += coeff[i][j] * old_basis[j][3];
      //std::cout << "i=" << i << " j=" << j << " test=" << test << " det=" << det << std::endl;
      
    }
    if (pol_checks) {
      if (((fabs(zero_comp.real()) > 1e-8 || fabs(array_new_basis[i][0].real()) > 1e-8) &&
           ((array_new_basis[i][0] - zero_comp).real() > fabs(zero_comp.real()) * 1e-8)) ||
          ((fabs(zero_comp.imag()) > 1e-8 || fabs(array_new_basis[i][0].imag()) > 1e-8) &&
           ((array_new_basis[i][0] - zero_comp).imag() > fabs(zero_comp.imag()) * 1e-8))) {
        std::cout << "Polarization_Warning in " << METHOD <<
                  ": Testing polarisation spinor transformation's system of equations - "
                  "zeroth component of polarization spinor failed..." << std::endl;
        msg_Out() << "zeroth component of the new polarization spinr calculated from linear superposition "
                  << std::setprecision(20) << array_new_basis[i][0] << std::endl;
        msg_Out() << "has to be the same as the zeroth component in new "
                     "polarization spinor" << std::setprecision(20) <<
                  zero_comp << std::endl;
      }
      if (((fabs(one_comp.real()) > 1e-8 || fabs(array_new_basis[i][1].real()) > 1e-8) &&
           ((array_new_basis[i][1] - one_comp).real() > fabs(one_comp.real()) * 1e-8)) ||
          ((fabs(one_comp.imag()) > 1e-8 || fabs(array_new_basis[i][1].imag()) > 1e-8) &&
           ((array_new_basis[i][1] - one_comp).imag() > fabs(one_comp.imag()) * 1e-8))) {
        std::cout << "Polarization_Warning in " << METHOD <<
                  ": Testing polarisation spinor transformation's system of equations - "
                  "first component of polarization spinor failed..." << std::endl;
        msg_Out() << "first component of the new polarization spinor calculated from linear superposition "
                  << std::setprecision(20) << array_new_basis[i][1] << std::endl;
        msg_Out() << "has to be the same as the first component in new polarization spinor" <<
                  std::setprecision(20) << one_comp << std::endl;
      }
    }
    if (pol_checks) {
      if (std::abs(array_new_basis[i][2] - two_comp)   > 1e-8 * std::max(1.0, std::abs(array_new_basis[i][2])) ||
      std::abs(array_new_basis[i][3] - three_comp) > 1e-8 * std::max(1.0, std::abs(array_new_basis[i][3]))) {
      msg_Out() << "Polarization_Warning in " << METHOD << ": components 2/3 of the transformed spinor do not match" << std::endl;
      }
    } 
  }
 
  return coeff;
}

vector<vector<Complex>> Polarization_Vector::BasisTrafo(const Polarization_Vector& old_basis, bool pol_checks) const {
  // TODO: Generalization to other particles than massive spin 1 particles
  // TODO: Avoid matrix formulation
  int m_nhel = old_basis.size();
  vector<vector<Complex>> coeff(m_nhel, vector<Complex>(m_nhel, Complex(0, 0)));
  if (m_nhel != this->size()){
    THROW(fatal_error, "polarisation vectors in new and old basis describes particles with different helicity "
                       "degrees of freedom!")
  }
  if (m_nhel != 3){
    THROW(not_implemented, "basis trafo not implemented for others then massive vector bosons")
  }
  std::array<std::array<Complex,3>, 3> array_old_basis;
  std::array<std::array<Complex,3>, 3> array_old_basis_inv;
  std::array<std::array<Complex,3>, 3> array_new_basis;
  // Inverse of matrix of old basis, using explizit form for 3x3 matrix inverse (polarisation components as rows)
  // determinante of 3*3 matrix (Sarrus rule)
  Complex det = old_basis[0][1] * old_basis[1][2] * old_basis[2][3]
                + old_basis[0][2] * old_basis[1][3] * old_basis[2][1]
                + old_basis[0][3] * old_basis[1][1] * old_basis[2][2]
                - old_basis[0][3] * old_basis[1][2] * old_basis[2][1]
                - old_basis[0][1] * old_basis[1][3] * old_basis[2][2]
                - old_basis[0][2] * old_basis[1][1] * old_basis[2][3];
  if (abs(det.real()) < 1e-8 && abs(det.imag()) < 1e-8){
    std::cout << "determinante is" << det << std::endl;
    THROW(fatal_error, "matrix of polarisation vectors in old spin basis is not invertible")
  }
  // inverse
  array_old_basis_inv[0][0]= (Complex(1, 0) / det)*(old_basis[1][2]*old_basis[2][3]-old_basis[1][3]*old_basis[2][2]);
  array_old_basis_inv[0][1]= (Complex(1, 0) / det)*(old_basis[0][3]*old_basis[2][2]-old_basis[0][2]*old_basis[2][3]);
  array_old_basis_inv[0][2]= (Complex(1, 0) / det)*(old_basis[0][2]*old_basis[1][3]-old_basis[0][3]*old_basis[1][2]);
  array_old_basis_inv[1][0]= (Complex(1, 0) / det)*(old_basis[1][3]*old_basis[2][1]-old_basis[1][1]*old_basis[2][3]);
  array_old_basis_inv[1][1]= (Complex(1, 0) / det)*(old_basis[0][1]*old_basis[2][3]-old_basis[0][3]*old_basis[2][1]);
  array_old_basis_inv[1][2]= (Complex(1, 0) / det)*(old_basis[0][3]*old_basis[1][1]-old_basis[0][1]*old_basis[1][3]);
  array_old_basis_inv[2][0]= (Complex(1, 0) / det)*(old_basis[1][1]*old_basis[2][2]-old_basis[1][2]*old_basis[2][1]);
  array_old_basis_inv[2][1]= (Complex(1, 0) / det)*(old_basis[0][2]*old_basis[2][1]-old_basis[0][1]*old_basis[2][2]);
  array_old_basis_inv[2][2]= (Complex(1, 0) / det)*(old_basis[0][1]*old_basis[1][2]-old_basis[0][2]*old_basis[1][1]);

  // formulate polarisation vector in new and old basis as array for easier matrix multiplikation afterwards
  for (size_t i(0); i<m_nhel; i++){
    for (size_t j(0); j<m_nhel; j++){
      array_new_basis[i][j] = (*this)[i][j+1];
      array_old_basis[i][j] = old_basis[i][j+1];
    }
  }
  // calculate basis transformation matrix: coeff = new_basis * (old_basis)^-1; each row contains the coeff for one
  // polarisation vector in new basis
  Complex tmp(0, 0);
  Complex test(0, 0);
  for (size_t i(0); i<m_nhel; i++){
    Complex zero_comp(0, 0), one_comp(0, 0), two_comp(0, 0), three_comp(0, 0);
    for (size_t j(0); j<m_nhel; j++){
      tmp = Complex(0, 0);
      test = Complex(0, 0);
      for (size_t k(0); k<m_nhel; k++){
        tmp += (array_new_basis[i][k] * array_old_basis_inv[k][j]);
        test += (array_old_basis_inv[i][k] * array_old_basis[k][j]);
      }

      coeff[i][j] = tmp;

      // test, whether inverse matrix is calculated correctly
      if ((i==j && (abs(test.real()-1)>1e-8 || abs(test.imag())>1e-8)) || (i!=j &&
                                                                           (abs(test.real())>1e-8 || abs(test.imag())>1e-8))){
        THROW(fatal_error, "multiplication of matrix and its inverse does not result in identity matrix")
      }
      // test, whether transformation also works for the zeroth component of the polarisation vectors
      zero_comp += coeff[i][j] * old_basis[j][0];
      // test system of equation's solution
      one_comp += coeff[i][j] * old_basis[j][1];
      two_comp += coeff[i][j] * old_basis[j][2];
      three_comp += coeff[i][j] * old_basis[j][3];
    }
    if (pol_checks) {
      if (((fabs(zero_comp.real()) > 1e-8 || fabs((*this)[i][0].real()) > 1e-8) &&
           (((*this)[i][0] - zero_comp).real() > fabs(zero_comp.real()) * 1e-8)) ||
          ((fabs(zero_comp.imag()) > 1e-8 || fabs((*this)[i][0].imag()) > 1e-8) &&
           (((*this)[i][0] - zero_comp).imag() > fabs(zero_comp.imag()) * 1e-8))) {
        std::cout << "Polarization_Warning in " << METHOD <<
                  ": Testing polarisation vector transformation's system of equations - "
                  "zeroth component of polarization vector failed..." << std::endl;
        msg_Out() << "zeroth component of the new polarization vector calculated from linear superposition "
                  << std::setprecision(20) << (*this)[i][0] << std::endl;
        msg_Out() << "has to be the same as the zeroth component in new "
                     "polarization vector" << std::setprecision(20) <<
                  zero_comp << std::endl;
      }
      if (((fabs(one_comp.real()) > 1e-8 || fabs((*this)[i][1].real()) > 1e-8) &&
           (((*this)[i][1] - one_comp).real() > fabs(one_comp.real()) * 1e-8)) ||
          ((fabs(one_comp.imag()) > 1e-8 || fabs((*this)[i][1].imag()) > 1e-8) &&
           (((*this)[i][1] - one_comp).imag() > fabs(one_comp.imag()) * 1e-8))) {
        std::cout << "Polarization_Warning in " << METHOD <<
                  ": Testing polarisation vector transformation's system of equations - "
                  "first component of polarization vector failed..." << std::endl;
        msg_Out() << "first component of the new polarization vector calculated from linear superposition "
                  << std::setprecision(20) << (*this)[i][1] << std::endl;
        msg_Out() << "has to be the same as the first component in new polarization vector" <<
                  std::setprecision(20) << one_comp << std::endl;
      }
      if (((fabs(two_comp.real()) > 1e-8 || fabs((*this)[i][2].real()) > 1e-8) &&
           (((*this)[i][2] - two_comp).real() > fabs(two_comp.real()) * 1e-8)) ||
          ((fabs(two_comp.imag()) > 1e-8 || fabs((*this)[i][2].imag()) > 1e-8) &&
           (((*this)[i][2] - two_comp).imag() > fabs(two_comp.imag()) * 1e-8))) {
        std::cout << "Polarization_Warning in " << METHOD <<
                  ": Testing polarisation vector transformation's system of equations - "
                  "second component of polarization vector failed..." << std::endl;
        msg_Out() << "second component of the new polarization vector calculated from linear superposition "
                  << std::setprecision(20) << (*this)[i][2] << std::endl;
        msg_Out() << "has to be the same as the second component in new polarization vector" <<
                  std::setprecision(20) << two_comp << std::endl;
      }
      if (((fabs(three_comp.real()) > 1e-8 || fabs((*this)[i][3].real()) > 1e-8) &&
           (((*this)[i][3] - three_comp).real() > fabs(three_comp.real()) * 1e-8)) ||
          ((fabs(three_comp.imag()) > 1e-8 || fabs((*this)[i][3].imag()) > 1e-8) &&
           (((*this)[i][3] - three_comp).imag() > fabs(three_comp.imag()) * 1e-8))) {
        std::cout << "Polarization_Warning in " << METHOD <<
                  ": Testing polarisation vector transformation's system of equations - "
                  "third component of polarization vector failed... " << std::endl;
        msg_Out() << "third component of the new polarization vector calculated from linear superposition "
                  << std::setprecision(20) << (*this)[i][3] << std::endl;
        msg_Out() << "has to be the same as the third component in new polarization vector" <<
                  std::setprecision(20) << three_comp << std::endl;
      }
    }
  }

  // test, whether transformation matrix is unitary
  if (pol_checks){
    for (size_t i(0); i<m_nhel; i++){
      for (size_t j(0); j<m_nhel; j++){
        tmp = Complex(0,0);
        for (size_t k(0); k<m_nhel; k++){
          tmp += coeff[i][k] * conj(coeff[j][k]);
        }
        if ((i==j && (abs((tmp-Complex(1, 0)).real()) > 1e-8 || abs(tmp.imag()) > 1e-8)) ||
            (i!=j && (abs(tmp.real()) > 1e-8 || abs(tmp.imag()) > 1e-8))) {
          std::cout<<"Polarization_Warning in "<< METHOD <<
                   ": Testing unitarity of transformation matrix failed..." << std::endl;
          msg_Out() << "failed for matrix entry" << i << j << "of unit matrix" << std::endl;
          msg_Out() << std::setprecision(20) << "actual value of this entry is " << tmp << std::endl;
        }
      }
    }
  }
  return coeff;
}

void Polarization_Tensor::Test(Vec4D p)
{
  bool massless(IsZero(p.Mass(),1e-6));
  if (massless) {
    THROW(not_implemented, "Massless polarisation tensors not implemented.");
  }
  
  bool success = true;
  std::cout<<METHOD<<": Testing symmetry..."<<std::endl;
  for(int s=0;s<5;s++) {
    for(int mu=0;mu<4;mu++) {
      for(int nu=mu+1;nu<4;nu++) {
        if(!IsEqual((*this)[s][mu][nu],(*this)[s][nu][mu])) {
        success=false;
        msg_Out()<<"s="<<s<<std::endl;
        msg_Out()<<"  eps["<<mu<<"]["<<nu<<"]="<<(*this)[s][mu][nu]<<std::endl
                 <<"  eps["<<nu<<"]["<<mu<<"]="<<(*this)[s][nu][mu]<<std::endl;
        }
      }
    }
  }
  std::cout<<"... "<<(success?"passed":"failed")<<std::endl<<std::endl;

  success=true;
  std::cout<<METHOD<<": Testing transversality..."<<std::endl;
  for(int s=0;s<5;s++) {
    Vec4C test1;
    for(int nu=0;nu<4;nu++) {
      for(int mu=0;mu<4;mu++) {
        double pmu = mu==0 ? p[mu] : -p[mu];
        test1[nu] += (*this)[s][mu][nu]*pmu;
      }
      if(!IsZero(test1[nu])) {
        success=false;
        msg_Out()<<"s="<<s<<std::endl;
        msg_Out()<<"  test1["<<nu<<"]="<<test1[nu]<<std::endl;
      }
    }
  }
  std::cout<<"... "<<(success?"passed":"failed")<<std::endl<<std::endl;

  success=true;
  std::cout<<METHOD<<": Testing tracelessness..."<<std::endl;
  for(int s=0;s<5;s++) {
    Complex d(0.0,0.0);
    for(int mu=0;mu<4;mu++) {
      double facmu = mu==0 ? 1.0 : -1.0;
      // should the metric go here, or is it really a trace?
      d += facmu*(*this)[s][mu][mu];
    }
    if(!IsZero(d)) {
      success=false;
      msg_Out()<<"s="<<s<<std::endl;
      msg_Out()<<"  d="<<d<<std::endl;
    }
  }
  std::cout<<"... "<<(success?"passed":"failed")<<std::endl<<std::endl;

  success=true;
  std::cout<<METHOD<<": Testing orthonormality..."<<std::endl;
  for(int s=0;s<5;s++) {
    for(int sp=0;sp<5;sp++) {
      Complex d(0.0,0.0);
      for(int nu=0;nu<4;nu++) {
        for(int mu=0;mu<4;mu++) {
          double facmu = mu==0 ? 1.0 : -1.0;
          double facnu = nu==0 ? 1.0 : -1.0;
          d += facmu*facnu*(*this)[s][mu][nu]*(*this)[sp].Conjugate()[mu][nu];
        }
      }
      if( (s==sp && !IsEqual(d,Complex(1.0,0.0))) ||
          (s!=sp && !IsZero(d))) {
        success=false;
        msg_Out()<<"s="<<s<<" sp="<<sp<<std::endl;
        msg_Out()<<"  d="<<d<<" should be: "<<(s==sp?1.0:0.0)<<std::endl;
      }
    }
  }
  std::cout<<"... "<<(success?"passed":"failed")<<std::endl<<std::endl;

  success=true;
  std::cout<<METHOD<<": Testing completeness..."<<std::endl;
  for(int mu=0;mu<4;mu++) {
    for(int nu=0;nu<4;nu++) {
      for(int rho=0;rho<4;rho++) {
        for(int sigma=0;sigma<4;sigma++) {
          Complex d(0.0,0.0);
          for(int s=0;s<5;s++) {
            d+=(*this)[s][mu][nu]*(*this)[s].Conjugate()[rho][sigma];
          }
          Complex ref=S(mu,rho,p)*S(nu,sigma,p)/2.0+
                      S(mu,sigma,p)*S(nu,rho,p)/2.0-
                      S(mu,nu,p)*S(rho,sigma,p)/3.0;
          if(!IsEqual(d,ref)) {
            success=false;
            msg_Out()<<"  mu="<<mu<<std::endl;
            msg_Out()<<"    nu="<<nu<<std::endl;
            msg_Out()<<"      rho="<<rho<<std::endl;
            msg_Out()<<"        sigma="<<sigma<<std::endl;
            msg_Out()<<"          d="<<d<<std::endl;
            msg_Out()<<"          r="<<ref<<std::endl;
          }
        }
      }
    }
  }
  std::cout<<"... "<<(success?"passed":"failed")<<std::endl<<std::endl;
}


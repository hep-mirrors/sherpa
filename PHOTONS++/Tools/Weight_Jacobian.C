#include "PHOTONS++/Tools/Weight_Jacobian.H"
#include <math.h>

using namespace PHOTONS;
using namespace std;


// class Weight_Jacobian
Weight_Jacobian::Weight_Jacobian() {
}

Weight_Jacobian::~Weight_Jacobian() {
}

Vec4D Weight_Jacobian::CalculateMomentumSum(Particle_Vector pv) {
  Vec4D sum = Vec4D(0.,0.,0.,0.);
  for (unsigned int i=0; i<pv.size(); i++) {
    sum = sum + pv.at(i)->Momentum();
  }
  return sum;
}


// class Weight_Jacobian_Lorentz
Weight_Jacobian_Lorentz::Weight_Jacobian_Lorentz(Particle_Vector newdip, Particle_Vector newspec, Particle_Vector olddip, Particle_Vector oldspec, Particle_Vector phot, Dipole_Type::code dtype) {
  m_dtype   = dtype;
  if (m_dtype == Dipole_Type::ff) {
    // only final state dipole constituents
    m_QC0      = CalculateMomentumSum(olddip)[0];
    m_PC0      = CalculateMomentumSum(newdip)[0];
    m_mMQ   = m_QC0;
    m_mMP   = m_PC0;
  }
  else if (m_dtype == Dipole_Type::fi) {
    // initial state dipole constituent at position 0
    m_QC0      = (CalculateMomentumSum(olddip) - olddip.at(0)->Momentum())[0];
    m_PC0      = (CalculateMomentumSum(newdip) - newdip.at(0)->Momentum())[0];
    m_mMQ   = CalculateMomentumSum(olddip)[0];
    m_mMP   = CalculateMomentumSum(newdip)[0];
  }
  m_QN0     = CalculateMomentumSum(oldspec)[0];
  m_PN0     = CalculateMomentumSum(newspec)[0];
  m_K0      = CalculateMomentumSum(phot)[0];
  CalculateWeight();
  CalculateMax();
}

Weight_Jacobian_Lorentz::~Weight_Jacobian_Lorentz() {
}

void Weight_Jacobian_Lorentz::CalculateWeight() {
  m_weight = pow(m_mMP,3)/pow(m_mMQ,3) * (m_QC0+m_QN0)/(m_PC0 + m_PN0 + m_K0);
}

void Weight_Jacobian_Lorentz::CalculateMax() {
  m_maxweight = 1;
}



// class Weight_Jacobian_Mapping
Weight_Jacobian_Mapping::Weight_Jacobian_Mapping(Particle_Vector newdip, Particle_Vector newspec, Particle_Vector olddip, Particle_Vector oldspec, Vec4D phot, double M, double u, Dipole_Type::code dtype) {
  m_dtype         = dtype;
  m_newdipole     = newdip;
  m_olddipole     = olddip;
  m_newspectator  = newspec;
  m_oldspectator  = oldspec;
  m_K             = Vec3D(phot);
  m_M             = M;
  m_u             = u;
  CalculateWeight();
  CalculateMax();
}

Weight_Jacobian_Mapping::~Weight_Jacobian_Mapping() {
}

void Weight_Jacobian_Mapping::CalculateWeight() {
  Vec3D Q         = Vec3D(0.,0.,0.);
  Vec3D P         = Vec3D(0.,0.,0.);
  Vec3D QN        = Vec3D(0.,0.,0.);
  Vec3D PN        = Vec3D(0.,0.,0.);
  double sumq     = 0;
  double sump     = 0;
  double prod     = 1;
  int N           = m_newdipole.size() + m_newspectator.size();
  for (unsigned int i=0; i<m_newdipole.size(); i++) {
    // exclude initial state particle (in initial-final dipoles at position 0)
    if (!((m_dtype == Dipole_Type::fi) && (i == 0))) {
      sumq = sumq -  Vec3D(m_olddipole.at(i)->Momentum())*Vec3D(m_olddipole.at(i)->Momentum())/m_olddipole.at(i)->Momentum()[0];
      sump = sump -  Vec3D(m_newdipole.at(i)->Momentum())*Vec3D(m_olddipole.at(i)->Momentum())/m_newdipole.at(i)->Momentum()[0];
      prod = prod * m_olddipole.at(i)->Momentum()[0]/m_newdipole.at(i)->Momentum()[0];
      Q = Q + Vec3D(m_olddipole.at(i)->Momentum());
      P = P + Vec3D(m_newdipole.at(i)->Momentum());
    }
  }
  for (unsigned int i=0; i<m_newspectator.size(); i++) {
    sumq = sumq -  Vec3D(m_oldspectator.at(i)->Momentum())*Vec3D(m_oldspectator.at(i)->Momentum())/m_oldspectator.at(i)->Momentum()[0];
    sump = sump -  Vec3D(m_newspectator.at(i)->Momentum())*Vec3D(m_oldspectator.at(i)->Momentum())/m_newspectator.at(i)->Momentum()[0];
    prod = prod * m_oldspectator.at(i)->Momentum()[0]/m_newspectator.at(i)->Momentum()[0];
    QN = QN + Vec3D(m_oldspectator.at(i)->Momentum());
    PN = PN + Vec3D(m_newspectator.at(i)->Momentum());
  }
  if (m_dtype == Dipole_Type::ff) {
    sumq = sumq + QN*QN/sqrt(m_M*m_M + QN*QN);
    sump = sump + (PN+m_K)*QN/sqrt(m_M*m_M + (PN+m_K)*(PN+m_K));
  }
  else if (m_dtype == Dipole_Type::fi) {
    sumq = sumq + Q*Q/sqrt(m_M*m_M + Q*Q);
    sump = sump + P*Q/sqrt(m_M*m_M + P*P);
    // subtract initial state particle (position 0 in old/newdipole)
    N = N-1;
  }
  else {
    prod = 0;
  }
  m_weight = pow(m_u,3*N-4) * prod * sumq/sump;
}

void Weight_Jacobian_Mapping::CalculateMax() {
  m_maxweight = 1;
}

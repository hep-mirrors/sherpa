#include "Poincare.H"
#include "MathTools.H"

#ifdef DEBUG__Poincare
#include <iostream>
#endif

using namespace ATOOLS;

Poincare::Poincare(): m_status(0),m_beta(1.,0.,0.,0.),m_rsq(1.)        // standard constructor Unity
{
  for (int i=0;i<4;i++)       // unitary
    for (int j=0;j<4;j++) {
      if (i==j) m_mat[i][j]=1.; else m_mat[i][j]=0.;
    }
}

Poincare::Poincare(Vec4D v1, Vec4D v2)   // rotations constructor
{
  /*
    Rotation of a 4-vector:  v2 -> v1 !!!!                                    
                             p1 = mat*p2                          
  */
  //intitialisation of rotation matrix
  m_status = 2;    
  short int i,k,l;
  double r[2][3][3],pm[2],sp[2],cp[2],st[2],ct[2];
  int n0,n1,n2;
  Vec4D pp[2];
  // normalisation
  pm[0] = Vec3D(v1).Abs();
  pm[1] = Vec3D(v2).Abs();
  pp[0] = (1./pm[0])*v1;
  pp[1] = (1./pm[1])*v2;
  // look for smaller cs to improve numeric
  if ((dabs(pp[0][3])<0.7)&&(dabs(pp[1][3])<0.7)) {
    n2 = 2;       // 3rd component
    n1 = 1;       // 2nd component
    n0 = 0;       // 1st component
  } 
  else if ((dabs(pp[1][2])<0.7)&&(dabs(pp[1][2])<0.7)) {
    n2 = 1;      
    n1 = 0;      
    n0 = 2;      
  } 
  else {
    n2 = 0;
    n1 = 2;
    n0 = 1;
  }
  for (i=0;i<2;i++) {
    // set angles
    ct[i] = pp[i][n2+1];
    st[i] = sqrt(1.-sqr(ct[i]));
    cp[i] = pp[i][n1+1]/st[i];
    sp[i] = pp[i][n0+1]/st[i];
#ifdef DEBUG__Poincare
    std::cout<<"Poincare::Poincare("<<v1<<","<<v2<<"): \n"
	     <<"cos(\\theta)="<<ct[i]<<" \n"
	     <<"sin(\\theta)="<<st[i]<<" \n"
	     <<"cos(\\phi)="<<cp[i]<<" \n"
	     <<"sin(\\phi)="<<sp[i]<<std::endl
	     <<"\\phi   -> "<<(sqr(sp[i])+sqr(cp[i])-1.0)<<std::endl
	     <<"\\theta -> "<<(sqr(st[i])+sqr(ct[i])-1.0)<<std::endl;
#endif
    // calculate matrix
    r[i][n0][n0] = cp[i]; 
    r[i][n0][n1] = sp[i]*ct[i]; 
    r[i][n0][n2] = st[i]*sp[i]; 
    r[i][n1][n0] = -sp[i]; 
    r[i][n1][n1] = ct[i]*cp[i];   
    r[i][n1][n2] = cp[i]*st[i];
    r[i][n2][n0] = 0.;
    r[i][n2][n1] = -st[i]; 
    r[i][n2][n2] = ct[i];
  }
  for (i=0;i<3;i++) {
    for (l=0;l<3;l++) {
      m_mat[i+1][l+1] = 0.;
      for (k=0;k<3;k++) 
	m_mat[i+1][l+1] += r[0][i][k]*r[1][l][k]; 
    }
  } 
}

Poincare::Poincare(Vec4D v)              // boost constuctor 
{
  m_status = 1;
  m_beta   = v; 
  m_rsq    = sqrt(m_beta.Abs2());
}
 
void Poincare::Boost(Vec4D& v)          // boosts vectors in CMS
{
  double v0  = (m_beta[0]*v[0]-Vec3D(m_beta)*Vec3D(v))/m_rsq;   // was m_beta*p/rsq before.
  double c1  = (v[0]+v0)/(m_rsq+m_beta[0]);
  v          = Vec4D(v0,Vec3D(v)-c1*Vec3D(m_beta)); 
}

void Poincare::BoostBack(Vec4D& v)      // boost back to LAB Frame
{
  double v0  = (m_beta[0]*v[0]+Vec3D(m_beta)*Vec3D(v))/m_rsq;
  double c1  = (v[0]+v0)/(m_rsq+m_beta[0]);
  v          = Vec4D(v0,Vec3D(v)+c1*Vec3D(m_beta));  
}


void Poincare::Rotate(Vec4D& v)           // rotate 
{
  Vec4D p2  = v;
  for (int i=1;i<4;i++) {
    v[i]    = 0.;
    for (int j=1;j<4;j++) 
      v[i] += m_mat[j][i]*p2[j];
  }
}

void Poincare::RotateBack(Vec4D& v)         // rotate back
{
  Vec4D p2  = v;
  for (int i=1;i<4;i++) {
    v[i]    = 0.;
    for (int j=1;j<4;j++) 
      v[i] += m_mat[i][j]*p2[j];
  }
}

bool Poincare::CheckBoost() 
{
  if (m_beta.Abs2()<=0.) return false;
  return true;
}

bool Poincare::CheckRotation() 
{
  for (short unsigned i=0;i<3;++i) {
    for (short unsigned j=0;j<3;++j) {
      double test=dabs(m_mat[i+1][j+1]);
      if (test>1.0 || !(test>=0.0)) return false; 
    }
  } 
  return true;
}

void Poincare::Invert() 
{
  for (int i=1;i<4;++i) m_beta[i]=-m_beta[i];

  for (int i=1;i<4;i++) {
    for (int j=i+1;j<4;j++) {
      double help=m_mat[i][j];
      m_mat[i][j]=m_mat[j][i];
      m_mat[j][i]=help;
    }
  }
}

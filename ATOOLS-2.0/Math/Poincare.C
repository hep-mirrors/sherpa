#include "Poincare.H"
#include "MathTools.H"

using namespace AMATOOLS;

Poincare::Poincare()                     // standard constructor Unity
{
  status=0;
  beta=vec4d(1.,0.,0.,0.);
  for (int i=0;i<4;i++)       // unitary
    for (int j=0;j<4;j++) {
      if (i==j) mat[i][j]=1.; else mat[i][j]=0.;
    }
}

Poincare::Poincare(vec4d v1, vec4d v2)   // rotations constructor
{
  /*
    Rotation of a 4-vector:  v2 -> v1 !!!!                                    
                             p1 = mat*p2                          
  */
  //intitialisation of rotation matrix
  //Poincare();
  status=2;    
  short int i,k,l;
  double r[2][3][3],pm[2],sp[2],cp[2],st[2],ct[2];
  int n0,n1,n2;
  vec4d pp[2];
  // normalisation
  pm[0] = vec3d(v1).abs();
  pm[1] = vec3d(v2).abs();
  pp[0] = (1./pm[0])*v1;
  pp[1] = (1./pm[1])*v2;
  // look for smaller cs to improve numeric
  if ((pp[0][3]<0.7)&&(pp[1][3]<0.7)) {
    n2=2;       // 3rd component
    n1=1;       // 2nd component
    n0=0;       // 1st component
  } 
  else if ((pp[1][2]<0.7)&&(pp[1][2]<0.7)) {
    n2=1;      
    n1=0;      
    n0=2;      
  } 
  else {
    n2=0;
    n1=2;
    n0=1;
  }
  for (i=0;i<2;i++) {
    // set angles
    ct[i] = pp[i][n2+1];
    st[i] = sqrt(1.-sqr(ct[i]));
    cp[i] = pp[i][n1+1]/st[i];
    sp[i] = pp[i][n0+1]/st[i];
    // calculate matrix
    r[i][n0][n0]=  cp[i]; 
    r[i][n0][n1]=  sp[i]*ct[i]; 
    r[i][n0][n2]=  st[i]*sp[i]; 
    r[i][n1][n0]=  -sp[i]; 
    r[i][n1][n1]=  ct[i]*cp[i];   
    r[i][n1][n2]=  cp[i]*st[i];
    r[i][n2][n0]=  0.;
    r[i][n2][n1]=  -st[i]; 
    r[i][n2][n2]=  ct[i];
  }
  for (i=0;i<3;i++) {
    for (l=0;l<3;l++) {
      mat[i+1][l+1] = 0.;
      for (k=0;k<3;k++) 
	mat[i+1][l+1] += r[0][i][k]*r[1][l][k]; 
    }
  } 
}

Poincare::Poincare(vec4d v)              // boost constuctor 
{
  status = 1;
  beta   = v; 
}
 
void Poincare::Boost(vec4d& v)          // boosts vectors in CMS
{
  vec4d ph   = v;
  double rsq = sqrt(beta.abs2());
  v[0]       = (beta[0]*ph[0]-vec3d(beta)*vec3d(ph))/rsq;   // was beta*p/rsq before.
  double c1  = (ph[0]+v[0])/(rsq+beta[0]);
  v          = vec4d(v[0],vec3d(ph)-c1*vec3d(beta));  
}

void Poincare::BoostBack(vec4d& v)      // boost back to LAB Frame
{
  vec4d ph   = v;
  double rsq = sqrt(beta.abs2());
  v[0]       = (beta[0]*ph[0]+vec3d(beta)*vec3d(ph))/rsq;
  double c1  = (ph[0]+v[0])/(rsq+beta[0]);
  v          = vec4d(v[0],vec3d(ph)+c1*vec3d(beta));  
}


void Poincare::Rotate(vec4d& v)           // rotate 
{
  vec4d p2  = v;
  for (int i=1;i<4;i++) {
    v[i]    = 0.;
    for (int j=1;j<4;j++) 
      v[i] += mat[j][i]*p2[j];
  }
}

void Poincare::RotateBack(vec4d& v)         // rotate back
{
  vec4d p2  = v;
  for (int i=1;i<4;i++) {
    v[i]    = 0.;
    for (int j=1;j<4;j++) 
      v[i] += mat[i][j]*p2[j];
  }
}

#include "PHASIC++/Channels/Sarge.H"
#include "ATOOLS/Org/Run_Parameter.H"
#include "ATOOLS/Org/Message.H"
#include "ATOOLS/Math/Random.H"

using namespace PHASIC;
using namespace ATOOLS;

Sarge::Sarge(int _nin,int _nout) : 
  nin(_nin), nout(_nout), x(NULL), rot(NULL)
{
  n_xi  = 2*nout-4;
  x     = new double[n_xi+1];
  rot   = new double*[3];
  for (short int i=0;i<3;i++) rot[i] = new double[3];
}

Sarge::~Sarge()
{
  delete[] x;
  for (short int i=0;i<3;i++) delete[] rot[i];
  delete[] rot;
}

void Sarge::GenerateWeight(Vec4D * p,Cut_Data * cuts)
{
  Vec4D sump(0.,0.,0.,0.);
  CalculateS0(cuts);  
  for (short int i=0;i<nin;i++) sump += p[i];
  double xi_min = sump.Abs2()/s0-((double)(nout)+1.)*((double)(nout)-2.)/2.;

  weight = 2.*pow(M_PI,nout-1)*pow(log(xi_min),2*nout-4)*(2.*nout-3.)/sqr(sump.Abs2());
  for (short int j=nin;j<nin+nout-1;j++) weight *= p[j]*p[j+1]; // denominator, product of 
  weight *= p[nin+nout-1]*p[nin];                              // momenta of Eq. 10
  weight *= 1./(pow(2.*M_PI,3*nout-4));
}

void Sarge::GeneratePoint(Vec4D* p,Cut_Data * cuts,double* ran)
{
  GeneratePoint(p,cuts);
}

void Sarge::GeneratePoint(Vec4D* p,Cut_Data * cuts)
{
  Vec4D sump(0.,0.,0.,0.);
  CalculateS0(cuts);  
  for (short int i=0;i<nin;i++) sump += p[i];
  double ET     = sqrt(sump.Abs2());
  double xi_min = sqr(ET)/s0-((double)(nout)+1.)*((double)(nout)-2.)/2.;

  Vec4D* q  = &p[nin];
  double costheta = 2.*ran->Get()-1.;
  double sintheta = sqrt(1.-costheta*costheta);
  double phi      = 2.*M_PI*ran->Get();
  q[0]      = ET/2.*Vec4D(1.,sintheta*::sin(phi),sintheta*cos(phi),costheta);
  q[nout-1] = ET/2.*Vec4D(1.,-sintheta*::sin(phi),-sintheta*cos(phi),-costheta);
  QcdAntenna(q,xi_min);
}

void Sarge::CalculateS0(Cut_Data * cuts) {
  s0 = 0.;
  for (short int i=0;i<cuts->ncut;i++) {
    for (short int j=i+1;j<cuts->ncut;j++) {
      if (s0<cuts->scut[i][j]) s0 = cuts->scut[i][j];
    }
  }
}

void Sarge::QcdAntenna(Vec4D* &p,double xi_min)
{
  Polytope(n_xi,x);
  double phi;
  double xi[2];
  Vec4D lab = p[0]+p[nout-1];
  double logm = log(xi_min); 
  for (short int j=0;j<nout-2;j++) {
    phi = 2.*M_PI*ran->Get();
    xi[0] = exp((x[2*j+1]-x[2*j])*logm);
    xi[1] = exp((x[2*j+2]-x[2*j])*logm);
    BasicAntenna(p[j],p[nout-1],p[j+1],xi,phi);
    lab += p[j+1];
  }
  
  double Elab  = sqrt(lab.Abs2());
  Vec3D  B     = (-1.)*Vec3D(lab)/Elab;
  double G     = lab[0]/Elab;
  double A     = 1./(1.+G);
  double scale = rpa->gen.Ecms()/Elab;
  double e,BQ;
  for (short int j=0;j<nout;j++) {
    e     = p[j][0];
    BQ    = B*Vec3D(p[j]);
    p[j]  = scale*Vec4D((G*e+BQ),Vec3D(p[j])+B*(e+A*BQ));
  }
}

void Sarge::BasicAntenna(Vec4D pin0,Vec4D pin1, Vec4D &k, double* xi, double phi)
{
  /* Basic antenna of SARGE, see hep-ph/0004047 */
  Vec4D  cms     = pin0+pin1;

  Vec4D ref(1.,0.,0.,1.);
  Vec4D pin0_cms;
  boost(1,cms,pin0_cms,pin0);
  double E1 = pin0_cms[0];

  rotat(0,pin0_cms,ref,rot);

  double k0_cms = E1*(xi[0]+xi[1]);
  double cos_cm = (xi[1]-xi[0])/(xi[1]+xi[0]);

  double sin_cm = sqrt(1.-sqr(cos_cm));
  Vec4D k_help  = k0_cms*Vec4D(1.,sin_cm*::sin(phi),sin_cm*cos(phi),cos_cm);
  Vec4D k_cms;
  rotat(1,k_cms,k_help,rot);
  boost(0,cms,k_cms,k);
}

void Sarge::PermP(int number,int* &perm)
{
  /* P-algorithm of Knuth "The Art of Computer Programming",
     shuffles a number of integers contained in perm  */
  int j,k,dummy;
  for (j=number-1;j>1;j--) {
    k = (int)(j*ran->Get())+1;
    dummy   = perm[j];
    perm[j] = perm[k];
    perm[k] = dummy;
  } 
}

void Sarge::Polytope(int m,double* &x)
{
  /* Produces a uniform random distribution inside a polytope with
     | x_k | < 1, | x_k-x_l | < 1, see physics/0003078 */
  int* perm = new int[m+1];
  short int i,k;
  for (i=1;i<m+1;i++) perm[i] = i;
  PermP(m+1,perm);
  // number of negative values
  k = (int)((m+1)*ran->Get());
  x[0] = 0.;
  if (k==0) {
    for (i=1;i<m+1;i++) x[perm[i]] = ran->Get();
    delete[] perm;
    return;
  }
  if (k==m) {
    for (i=1;i<m+1;i++) x[perm[i]] = -ran->Get();
    delete[] perm;
    return;
  }
  double v1,v2,prod,y1;
  prod = 1.;
  for (i=1;i<k+1;i++) prod *= ran->Get();
  v1 = -log(prod);
  prod = 1.;
  for (i=1;i<m-k+2;i++) prod *= ran->Get();
  v2 = -log(prod);                           // -> b(k)
  y1 = v1/(v1+v2);
  x[perm[1]]   = -y1;
  x[perm[m]] = (1-y1)*pow(ran->Get(),double(1./(m-k)));
  for (i=2;i<k+1;i++) x[perm[i]] = x[perm[1]]*ran->Get();
  for (i=k+1;i<m;i++) x[perm[i]] = x[perm[m]]*ran->Get();
  delete[] perm;
}

void Sarge::boost(int lflag,Vec4D q,Vec4D& ph,Vec4D& p)
{
  /*
                                        _                       
   Boost of a 4-vector ( relative speed q/q(0) ):               
                                                                
   ph is the 4-vector in the rest frame of q                    
   p is the corresponding 4-vector in the lab frame             
                                                                  
              INPUT                               OUTPUT         
                                                            
  lflag= 0:   q, ph                               p             
                                                                
  lflag= 1:   q, p                                ph            

  */

  double rsq = sqrt(q.Abs2());
  if (lflag==0) {
    p[0] = (q[0]*ph[0]+Vec3D(q)*Vec3D(ph))/rsq;
    double c1 = (ph[0]+p[0])/(rsq+q[0]);
    p = Vec4D(p[0],Vec3D(ph)+c1*Vec3D(q));  
  }
  else {
    ph[0] = q*p/rsq;
    double c1 = (p[0]+ph[0])/(rsq+q[0]);
    ph = Vec4D(ph[0],Vec3D(p)-c1*Vec3D(q));  
  }
}

void Sarge::rotat(int lflag,Vec4D& p1 ,Vec4D p2,double** _rot)
{
/*
   Rotation of a 4-vector:                                      
                                                                
                            p1= rot*p2                          
                                                               
              INPUT                               OUTPUT        
                                                                
  lflag= 0:   p1, p2  (even if |p1|.ne.|p2|)      rot           
                                                                
  lflag= 1:   p2, rot                             p1            
*/

  if (lflag==0) {
    short int i,k,l;
    double r[2][3][3],pm[2],sp[2],cp[2],st[2],ct[2];
    Vec4D pp[2];
    pm[0] = Vec3D(p1).Abs();
    pm[1] = Vec3D(p2).Abs();
    pp[0] = (1./pm[0])*p1;
    pp[1] = (1./pm[1])*p2;
    for (i=0;i<2;i++) {
      ct[i] = pp[i][3];
      st[i] = sqrt(1.-sqr(ct[i]));
      if (ATOOLS::IsEqual(dabs(ct[i]),1.)) {
	cp[i] = 1.;
	sp[i] = 0.;
      }
      else {
	cp[i] = pp[i][2]/st[i];
	sp[i] = pp[i][1]/st[i];
      }
      r[i][0][0]=  cp[i]; 
      r[i][0][1]=  sp[i]*ct[i]; 
      r[i][0][2]=  st[i]*sp[i]; 
      r[i][1][0]=  -sp[i]; 
      r[i][1][1]=  ct[i]*cp[i];   
      r[i][1][2]=  cp[i]*st[i];
      r[i][2][0]=  0.;
      r[i][2][1]=  -st[i]; 
      r[i][2][2]=  ct[i];
    }
    for (i=0;i<3;i++) {
      for (l=0;l<3;l++) {
	_rot[i][l] = 0.;
	for (k=0;k<3;k++) 
	  _rot[i][l] += r[0][i][k]*r[1][l][k];
      }
    }
  }
  else {
    short int i,j;
    p1[0] = p2[0];
    for (i=0;i<3;i++) {
      p1[i+1] = 0.;
      for (j=0;j<3;j++) {
	p1[i+1] += _rot[i][j]*p2[j+1];
      }
    }
  }
}



#include "AMEGIC++/Main/Helicity.H"
#include "ATOOLS/Math/MyComplex.H"
#include "ATOOLS/Math/Vector.H"
#include "ATOOLS/Math/MathTools.H"
#include "ATOOLS/Org/Message.H"
#include "AMEGIC++/Main/Pol_Info.H"
#include "ATOOLS/Org/Run_Parameter.H"
#include <algorithm>

using namespace ATOOLS;
using namespace AMEGIC;
using namespace std;

Helicity::Helicity(int Nin,int Nout,Flavour* fl,Pol_Info* pl) : 
  m_flavours(Nin+Nout, kf_none), m_nPols(Nin+Nout),
  m_allowTrafo(true), m_needsTrafo(false), m_spos(-1)  
{
 
  int N=Nin+Nout;
  p_pol_types = new char[N+1];
  p_angles    = new double[N];

  /* Store the informations. Look up massive polarized spinors.  */
  for(int i=0;i<N;i++){
    m_flavours[i]  = fl[i];
    p_pol_types[i] = pl[i].pol_type;
    p_angles[i]    = pl[i].angle;
    m_nPols[i]     = pl[i].num;
    if (p_pol_types[i]=='e') m_spos=i;
    
    // Treatment for massive polarized spinors
    if ( (m_flavours[i].IsMassive() && m_flavours[i].IsFermion())
	 && !(m_nPols[i]==2 && pl[i].factor[0] == pl[i].factor[1]) )  {
      m_trafoList.push_back(i);
      m_nPols[i] = 2;
      msg_Out() <<"Flavour "<<m_flavours[i]<<" has been added to the list of "
		<<"polarized massive fermions."<<endl;
    }
    // Treatment for flavours using only one polarisation amplitude
    else if ((pl[i].pol_type=='h' || pl[i].pol_type=='c') && (m_nPols[i]==1)) {
      if (pl[i].type[0]==mt::p_p) p_pol_types[i] = '+';
      if (pl[i].type[0]==mt::p_m) p_pol_types[i] = '-';
    }
  }
  m_trafoMatrices.resize(m_trafoList.size(), ATOOLS::CMatrix(2));
  p_pol_types[N] = 0;  // end of char-string


  // The total number of different polarisation combinations
  m_nsign = 1; for(int pNum=0; pNum<N; ++pNum) m_nsign *= m_nPols[pNum];


  /* Create the list of all nessecary polarisation combinations that need to be 
     evaluated to obtain the cross section for the specified process. Also calculate the
     respective weights of the physical helicity amplitudes */
  p_slist = new Sign_List[m_nsign];
  for (size_t i=0;i<m_nsign;i++) {
    p_slist[i].s = new int[N];
    for (int flav=0;flav<N;++flav){

      int div = 1;
      for (int k=0;k<flav;k++) div *= m_nPols[k]; 
      int l = (i/div) % m_nPols[flav];
      
      // flavour in question is a polarized massive fermion
      if (find(m_trafoList.begin(), m_trafoList.end(), flav) != m_trafoList.end()) {
        l==0 ? p_slist[i].s[flav] = -1 :  p_slist[i].s[flav] = +1;
	if (pl[flav].num==1) {
	  if (pl[flav].type[0] != p_slist[i].s[flav]) p_slist[i].polfactor=0.;
	  if(flav<Nin) p_slist[i].polfactor *= pl[flav].factor[0];
	} 
	else if(flav<Nin) p_slist[i].polfactor *= pl[flav].factor[l];  
      }
      // flavour in question does not cause a transformation of the amplitudes
      else {
	p_slist[i].s[flav] = pl[flav].type[l];
	if (pl[flav].num==1 && l==1) p_slist[i].polfactor = 0.;
	if(flav<Nin) p_slist[i].polfactor *= pl[flav].factor[l];  
	Tensor_Struc ts;
	p_slist[i].polfactor*=ts.GetTfactor(pl[flav].type[l]);  
      }          
     
    } 
  }
}

Helicity::~Helicity() 
{
  if (p_slist) {
    for (size_t i=0;i<m_nsign;i++) delete [] p_slist[i].s;
    delete [] p_slist;
  }
  delete [] p_pol_types;
  delete [] p_angles;
}

int Helicity::GetPartnerPol(const int heli,const int flav, int& lambda) // inline
{// only works for spinors (two possible polariations)!
  // get the "div"
  int div=1; 
  for (int i=0; i<flav; ++i) div *= m_nPols[i];   // -> to a list
  
  // get the sign and the oppsite sign of the polarisation
  lambda = (heli/div) % m_nPols[flav];
  int lambda2= 1-lambda;
  
  // obtain the corresponding helComb
  return heli + div*(lambda2 - lambda);
}

bool Helicity::IsContrib(int i,int* pm,int length)
{
  if (!pm) return 1;
  for (int j=0;j<length;j++) if(pm[j]<99 && p_slist[i].s[j]!=pm[j]) return 0;
  return 1;
}

int Helicity::Compare(Helicity* h_cmp, int N)
{
  if (MaxHel() != h_cmp->MaxHel()) return 0;
  for (size_t i=0;i<MaxHel();i++) {
    for (size_t j=0;j<size_t(N);j++) 
      if (p_slist[i].s[j] != (*h_cmp)[i][j]) return 0;
  }
  return 1;
}

void Helicity::InitializeSpinorTransformation(Basic_Sfuncs * BS)  
{
  if ((m_needsTrafo = m_trafoList.size()) > 0) {

    std::vector<int>::const_iterator flNum(m_trafoList.begin());
    for (int i=0; flNum != m_trafoList.end(); flNum++, ++i) {
   
      // initializing several phyiscal and non-physical values
      ATOOLS::Vec4D  k0(BS->Getk0());
      ATOOLS::Vec4D  p(BS->Momentum(*flNum));
      double         mu(BS->Mu(*flNum).real());
      double         m(m_flavours[*flNum].Mass());
      double         antiSign;
      m_flavours[*flNum].IsAnti() ? antiSign = -1. : antiSign = 1.;
      
      /* Constructing the polarization vector s. */
      double         p3abs(ATOOLS::Vec3D(p).Abs());
      ATOOLS::Vec4D  s = Vec4D(p3abs, p[0]/p3abs*Vec3D(p))/m;

      // Determine the common scaling factor c.
      double denominator((p + antiSign*m*s)*k0);
      if (ATOOLS::IsZero(denominator)) {
	msg_Error()<<"Warning: Encountered a zero-denominator while trying to "
	           <<"construct the matrices for the polarisation transformation."<<endl
	           <<"No transformation will occur."<<endl;
	m_needsTrafo=false;
	return; }
      Complex  c( csqrt( 2. * k0 * p / denominator ));    

      /* From the scaling factor, compute the transformation matrices using the S-funcs 
	 S(+,s,p) and S(-,s,p) from Basic_Sfuncs */
      std::pair<Complex, Complex> S = BS->GetS(s, *flNum);
      Complex eta_s(csqrt(2. * s * k0));

      Complex A( (k0 * p) / (k0 * s) );
      
      m_trafoMatrices[i][0][0] = 
	.5*c + (antiSign*c*0.25/m) * (m*m/A + A + S.first*S.second);
      m_trafoMatrices[i][0][1] = -0.5*antiSign*c*mu*eta_s*S.second/m;    
      m_trafoMatrices[i][1][0] = -0.5*antiSign*c*mu*eta_s*S.first/m; 
      m_trafoMatrices[i][1][1] = m_trafoMatrices[i][0][0];
    }
  }
}

void Helicity::SpinorTransformation(std::vector<Complex>& A)
{ 
  if (UseTransformation()) {
    std::vector<int>::const_iterator flav(m_trafoList.begin());
    for(int i=0; flav != m_trafoList.end(); flav++, ++i) {
      std::vector<Complex> Atemp(A);
      for (size_t helComb=0; helComb<MaxHel(); ++helComb) {
      	int lambda;
	int helComb2 = GetPartnerPol(helComb, *flav, lambda);  
	A[helComb] =    m_trafoMatrices[i][lambda][lambda  ] * Atemp[helComb]
	              + m_trafoMatrices[i][lambda][1-lambda] * Atemp[helComb2];
      }
    }
  }
}

size_t Helicity::MaxHel(size_t i) {
  return (size_t) m_nPols[i];
}

size_t Helicity::GetAmplitudeNumber(std::vector<int> *Helis)
{
  int mult(1);
  size_t num(0);
  for (size_t i=0; i < Helis->size(); ++i) {
    num += mult * (*Helis)[i];
    mult *= m_nPols[i];
  }
  return num;
}

int Helicity::GetPol(const int& flav, const int& hNumber)
{ 
  /* Return the polarisation state of flavour "flav" in helicity combination "hNumber" */
  return p_slist[hNumber].s[flav];
}

#include "METOOLS/Main/Partial_Amplitude_Base.H"
#include "ATOOLS/Org/Exception.H"
#include "ATOOLS/Org/MyStrStream.H"
#include <stdarg.h>

#include "METOOLS/Main/Three_Particle_Amplitudes.H"
#include "METOOLS/Main/Four_Particle_Amplitudes.H"

using namespace METOOLS;
using namespace ATOOLS;
using namespace std;

string METOOLS::GetName(Flavour* flavs, int n)
{
  string name=ToString(flavs[0])+" --> ";
  for(int i=1; i<n; ++i)
    name+=" "+ToString(flavs[i]);
  return name;
}

Partial_Amplitude_Base::Partial_Amplitude_Base(Flavour* flavs,size_t size,
                                               int* i,bool* out) :
  Spin_Structure<Complex>(flavs,size,i), 
  p_flavs(flavs), p_i(i), p_out(out)
{
}

Partial_Amplitude_Base::~Partial_Amplitude_Base()
{
  if(p_out) { delete [] p_out; }
  if(p_i)   { delete [] p_i; }
}

void Partial_Amplitude_Base::AssertSpins(int spin, ...)
{
  va_list ap;
  va_start(ap, spin);
  int sp=spin;
  for (size_t i(0); i<m_spins.size(); ++i) {
    if(p_flavs[p_i[i]].IntSpin()!=sp)
      THROW(fatal_error, ToString(p_flavs[p_i[i]])+" does not have spin "+
            ToString(sp)+" in "+GetName(p_flavs,m_spins.size())+".");
    sp=va_arg(ap,int);
  }
  va_end(ap);
}

void Partial_Amplitude_Base::AssertIn(int nin)
{
  int isin(0);
  for (size_t i(0);i<m_spins.size();i++) if (p_out[i]==false) isin++;
  if (isin!=nin) {
    THROW(fatal_error, "Expected "+ToString(nin)+" incomings, but got "
          +ToString(isin)+" in "+GetName(p_flavs,m_spins.size())+".");
  }
}

#define SELECT_ISOTROPIC \
  msg_Debugging()<<METHOD<<": Generic hadron decay ME for "<<GetName(flavs,n) \
  <<" not implemented. Using Isotropic."<<endl; \
  me=new Isotropic(flavs, n, inds, out)

Partial_Amplitude_Base* Partial_Amplitude_Base::Select(Flavour* flavs, int n)
{
  Partial_Amplitude_Base* me(NULL);
  int* inds = new int[n];
  for (int i=0; i<n; ++i) inds[i]=i;
  bool* out = new bool[n];
  out[0]=false;
  for(int i(1);i<n;++i) out[i]=true;

  for (int i=0; i<n; ++i) {
    if (flavs[i].IsVector() && IsZero(flavs[i].Mass())) {
      SELECT_ISOTROPIC;
      return me;
    }
  }

  if(flavs[0].IsScalar()) {
    if(n==3) {
      if(flavs[1].IsScalar() && flavs[2].IsScalar())
        me=new SSS(flavs, inds, out);
      else if(flavs[1].IsFermion() && flavs[2].IsFermion()) {
        inds[0]=0;
        inds[1]=1;
        inds[2]=2;
        if((flavs[inds[1]].IsLepton() || flavs[inds[2]].IsLepton()) &&
           (  flavs[inds[0]].Kfcode()==kf_pi_plus ||
              flavs[inds[0]].Kfcode()==kf_K_plus ||
              flavs[inds[0]].Kfcode()==kf_D_plus ||
              flavs[inds[0]].Kfcode()==kf_D_s_plus ||
              flavs[inds[0]].Kfcode()==kf_B_plus ||
              flavs[inds[0]].Kfcode()==kf_B_c
              ))
          me=new SFF_FPI(flavs, inds, out);
        else
          me=new SFF(flavs, inds, out, Complex(1.0,0.0), Complex(1.0,0.0));
      }
      else if(flavs[1].IsScalar() && flavs[2].IsVector())
        me=new SSV(flavs, inds, out);
      else if(flavs[2].IsScalar() && flavs[1].IsVector()) {
        inds[0]=0;
        inds[1]=2;
        inds[2]=1;
        me=new SSV(flavs, inds, out);
      }
      else if(flavs[1].IsVector() && flavs[2].IsVector())
        me=new SVV(flavs, inds, out);
      else if(flavs[1].IsScalar() && flavs[2].IsTensor()) {
        inds[0]=2;
        inds[1]=1;
        inds[2]=0;
        me=new TSS(flavs, inds, out);
      }
      else if(flavs[2].IsScalar() && flavs[1].IsTensor()) {
        inds[0]=1;
        inds[1]=2;
        inds[2]=0;
        me=new TSS(flavs, inds, out);
      }
      else if(flavs[1].IsVector() && flavs[2].IsTensor()) {
        inds[0]=2;
        inds[1]=1;
        inds[2]=0;
        me=new TVS(flavs, inds, out);
      }
      else if(flavs[2].IsVector() && flavs[1].IsTensor()) {
        inds[0]=1;
        inds[1]=2;
        inds[2]=0;
        me=new TVS(flavs, inds, out);
      }
      else {
        SELECT_ISOTROPIC;
      }
    }
    else {
      SELECT_ISOTROPIC;
    }
  }
  else if(flavs[0].IsFermion()) {
    if(n==3) {
      if(flavs[1].IsFermion() && flavs[2].IsScalar()) {
        inds[0]=2;
        inds[1]=1;
        inds[2]=0;
        if((flavs[inds[1]].IsLepton() || flavs[inds[2]].IsLepton()) &&
           (  flavs[inds[0]].Kfcode()==kf_pi_plus ||
           flavs[inds[0]].Kfcode()==kf_K_plus ||
           flavs[inds[0]].Kfcode()==kf_D_plus ||
           flavs[inds[0]].Kfcode()==kf_D_s_plus ||
           flavs[inds[0]].Kfcode()==kf_B_plus ||
              flavs[inds[0]].Kfcode()==kf_B_c
              ))
          me=new SFF_FPI(flavs, inds, out);
        else
          me=new SFF(flavs, inds, out, Complex(1.0,0.0), Complex(1.0,0.0));
      }
      else if(flavs[1].IsScalar() && flavs[2].IsFermion()) {
        inds[0]=1;
        inds[1]=2;
        inds[2]=0;
        if((flavs[inds[1]].IsLepton() || flavs[inds[2]].IsLepton()) &&
           (  flavs[inds[0]].Kfcode()==kf_pi_plus ||
           flavs[inds[0]].Kfcode()==kf_K_plus ||
           flavs[inds[0]].Kfcode()==kf_D_plus ||
           flavs[inds[0]].Kfcode()==kf_D_s_plus ||
           flavs[inds[0]].Kfcode()==kf_B_plus ||
              flavs[inds[0]].Kfcode()==kf_B_c
              ))
          me=new SFF_FPI(flavs, inds, out);
        else
          me=new SFF(flavs, inds, out, Complex(1.0,0.0), Complex(1.0,0.0));
      }
      else if(flavs[1].IsFermion() && flavs[2].IsVector()) {
        inds[0]=2;
        inds[1]=1;
        inds[2]=0;
        me=new VFF(flavs, inds, out, Complex(1.0,0.0), Complex(1.0,0.0));
      }
      else if(flavs[2].IsFermion() && flavs[1].IsVector()) {
        inds[0]=1;
        inds[1]=2;
        inds[2]=0;
        me=new VFF(flavs, inds, out, Complex(1.0,0.0), Complex(1.0,0.0));
      }
      else {
        SELECT_ISOTROPIC;
      }
    }
    else {
      SELECT_ISOTROPIC;
    }
  }
  else if(flavs[0].IsVector()) {
    if(n==3) {
      if(flavs[1].IsScalar() && flavs[2].IsScalar()) {
        inds[0]=1;
        inds[1]=2;
        inds[2]=0;
        me=new SSV(flavs, inds, out);
      }
      else if(flavs[1].IsFermion() && flavs[2].IsFermion())
        me=new VFF(flavs, inds, out, Complex(1.0,0.0), Complex(0.0,0.0));
      else if(flavs[1].IsVector() && flavs[2].IsVector())
        me=new VVV(flavs, inds, out);
      else if(flavs[2].IsScalar() && flavs[1].IsVector()) {
        inds[0]=2;
        inds[1]=0;
        inds[2]=1;
        me=new SVV(flavs, inds, out);
      }
      else if(flavs[1].IsScalar() && flavs[2].IsVector()) {
        inds[0]=1;
        inds[1]=0;
        inds[2]=2;
        me=new SVV(flavs, inds, out);
      }
      else if(flavs[1].IsVector() && flavs[2].IsTensor()) {
        inds[0]=2;
        inds[1]=1;
        inds[2]=0;
        me=new TVV(flavs, inds, out);
      }
      else if(flavs[2].IsVector() && flavs[1].IsTensor()) {
        inds[0]=1;
        inds[1]=2;
        inds[2]=0;
        me=new TVV(flavs, inds, out);
      }
      else {
        SELECT_ISOTROPIC;
      }
    }
    else if(n==4) {
      if(flavs[1].IsScalar() && flavs[2].IsScalar() && flavs[3].IsScalar()) {
        me=new VSSS(flavs, inds, out);
      }
      else {
        SELECT_ISOTROPIC;
      }
    }
    else {
      SELECT_ISOTROPIC;
    }
  }
  else if(flavs[0].IsRaritaSchwinger()) {
    SELECT_ISOTROPIC;
  }
  else if(flavs[0].IsTensor()) {
    if(n==3) {
      if(flavs[1].IsScalar() && flavs[2].IsScalar())
        me=new TSS(flavs, inds, out);
      else if(flavs[2].IsScalar() && flavs[1].IsVector())
        me=new TVS(flavs, inds, out);
      else if(flavs[1].IsScalar() && flavs[2].IsVector()) {
        inds[0]=0;
        inds[1]=2;
        inds[2]=1;
        me=new TVS(flavs, inds, out);
      }
      else if(flavs[1].IsVector() && flavs[2].IsVector())
        me=new TVV(flavs, inds, out);
      else {
        SELECT_ISOTROPIC;
      }
    }
    else {
      SELECT_ISOTROPIC;
    }
  }
  else {
    SELECT_ISOTROPIC;
  }
  return me;
}


Isotropic::Isotropic(Flavour *fl,int n,int* i,bool *out) :
  Partial_Amplitude_Base(fl,n,i,out)
{
  AssertIn(1);
}

void Isotropic::operator()(const Vec4D * moms,const bool anti)
{
  CreateTrivial(Complex(1.0,0.0));
}

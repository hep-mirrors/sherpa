#include "METOOLS/Main/Partial_Amplitude_Base.H"
#include "ATOOLS/Org/Exception.H"
#include "ATOOLS/Org/MyStrStream.H"
#include <stdarg.h>

#include "METOOLS/Main/Three_Particle_Amplitudes.H"
#include "METOOLS/Main/Four_Particle_Amplitudes.H"

using namespace METOOLS;
using namespace ATOOLS;
using namespace std;

namespace METOOLS {
  string GetName(const Flavour_Vector& flavs)
  {
    string name=ToString(flavs[0])+" --> ";
    for(size_t i=1; i<flavs.size(); ++i)
      name+=" "+ToString(flavs[i]);
    return name;
  } 
}

Partial_Amplitude_Base::Partial_Amplitude_Base(const Flavour_Vector& flavs,
                                               const vector<int>& i,
                                               const vector<bool>& out) :
  Spin_Amplitudes(flavs,i), 
  p_flavs(flavs), p_out(out)
{
  p_i.resize(i.size());
  for (size_t j=0; j<i.size(); ++j) p_i[j]=i[j];
  p_out.resize(out.size());
  for (size_t j=0; j<out.size(); ++j) p_out[j]=out[j];
}

Partial_Amplitude_Base::~Partial_Amplitude_Base()
{
}

void Partial_Amplitude_Base::AssertSpins(int spin, ...)
{
  va_list ap;
  va_start(ap, spin);
  int sp=spin;
  for (size_t i(0); i<m_spins.size(); ++i) {
    if(p_flavs[p_i[i]].IntSpin()!=sp)
      THROW(fatal_error, ToString(p_flavs[p_i[i]])+" does not have spin "+
            ToString(sp)+" in "+GetName(p_flavs)+".");
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
          +ToString(isin)+" in "+GetName(p_flavs)+".");
  }
}

#define SELECT_ISOTROPIC \
  msg_Debugging()<<METHOD<<": Generic hadron decay ME for "<<GetName(flavs) \
  <<" not implemented. Using Isotropic."<<endl; \
  me=new Isotropic(flavs, inds, out)

Partial_Amplitude_Base*
Partial_Amplitude_Base::Select(const Flavour_Vector& flavs)
{
  int n=flavs.size();
  Partial_Amplitude_Base* me(NULL);
  vector<int> inds(n);
  for (int i=0; i<n; ++i) inds[i]=i;
  vector<bool> out(n);
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


Isotropic::Isotropic(const Flavour_Vector& fl,
                     const vector<int>& i, const vector<bool>& out) :
  Partial_Amplitude_Base(fl,i,out)
{
  AssertIn(1);
}

void Isotropic::Calculate(const ATOOLS::Vec4D_Vector& moms, bool anti)
{
  CreateTrivial(Complex(1.0,0.0));
}

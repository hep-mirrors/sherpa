#include "METOOLS/Main/Three_Particle_Amplitudes.H"
#include "ATOOLS/Org/Message.H"
#include "ATOOLS/Org/Exception.H"

using namespace METOOLS;
using namespace ATOOLS;
using namespace std;

////////////////////////////////////////////////////////////////////////////////
////////////////////////////////////////////////////////////////////////////////
// 0 -> SSS
////////////////////////////////////////////////////////////////////////////////
////////////////////////////////////////////////////////////////////////////////

SSS::SSS(Flavour *fl,int* i,bool *out) :
  Partial_Amplitude_Base(fl,3,i,out)
{
  AssertIn(1);
  AssertSpins(0,0,0);
}

void SSS::operator()(const Vec4D * moms,const bool anti)
{
  Insert(Complex(1.0,0.0),0);
}


////////////////////////////////////////////////////////////////////////////////
////////////////////////////////////////////////////////////////////////////////
// 0 -> SFF
////////////////////////////////////////////////////////////////////////////////
////////////////////////////////////////////////////////////////////////////////


SFF::SFF(Flavour *fl,int* i,bool *out,const Complex cL,const Complex cR) :
  Partial_Amplitude_Base(fl,3,i,out), m_cL(cL), m_cR(cR), p_xyz(NULL) 
{
  AssertIn(1);
  AssertSpins(0,1,1);

  Flavour f1(p_flavs[p_i[1]]), f2(p_flavs[p_i[2]]);
  if ((f1.IsAnti()^p_out[p_i[1]])&&(f2.IsAnti()==p_out[p_i[2]])) {
    m_bar=1;
    m_nonbar=2;
  }
  else if ((f2.IsAnti()^p_out[p_i[2]])&&(f1.IsAnti()==p_out[p_i[1]])) {
    m_bar=2;
    m_nonbar=1;
  }
  else {
    msg_Error()<<METHOD<<": you provided an impossible combination of incoming/"
               <<"outgoing particles/antiparticles: "<<p_flavs[p_i[0]]
               <<" + "<<p_flavs[p_i[1]]<<" + "<<p_flavs[p_i[2]]<<endl;
  }
  p_xyz = new XYZFunc(3,fl,1,p_i);
}

void SFF::operator()(const Vec4D * moms,const bool anti)
{
  Complex amp(0.,0.);
  p_xyz->Prepare(moms,anti);
  for (int hel1(0);hel1<2;hel1++) {
    for (int hel2(0);hel2<2;hel2++) {
      amp = p_xyz->Y(m_bar,hel1,m_nonbar,hel2,m_cL,m_cR);
      vector<pair<int,int> > spins;
      spins.push_back(make_pair(0,0));
      spins.push_back(make_pair(1,hel1));
      spins.push_back(make_pair(2,hel2));
      Insert(amp,spins);
    }
  }
}


////////////////////////////////////////////////////////////////////////////////
////////////////////////////////////////////////////////////////////////////////
// 0 -> SFF_FPI
////////////////////////////////////////////////////////////////////////////////
////////////////////////////////////////////////////////////////////////////////


SFF_FPI::SFF_FPI(Flavour *fl,int *i,bool *out,const Complex c) :
  Partial_Amplitude_Base(fl,3,i,out), m_c(c), p_xyz(NULL) 
{
  AssertIn(1);
  AssertSpins(0,1,1);

  Flavour f1(p_flavs[p_i[1]]), f2(p_flavs[p_i[2]]);
  if ((f1.IsAnti()^p_out[p_i[1]])&&(f2.IsAnti()==p_out[p_i[2]])) {
    m_bar=1;
    m_nonbar=2;
  }
  else if ((f2.IsAnti()^p_out[p_i[2]])&&(f1.IsAnti()==p_out[p_i[1]])) {
    m_bar=2;
    m_nonbar=1;
  }
  else {
    msg_Error()<<METHOD<<": you provided an impossible combination of incoming/"
               <<"outgoing particles/antiparticles: "<<p_flavs[p_i[0]]
               <<" + "<<p_flavs[p_i[1]]<<" + "<<p_flavs[p_i[2]]<<endl;
  }
  p_xyz = new XYZFunc(3,fl,1,p_i);
}

void SFF_FPI::operator()(const Vec4D * moms,const bool anti)
{
  Vec4D pS(p_i[0]==0?moms[p_i[0]]:-moms[p_i[0]]);
  Complex amp(0.,0.);
  p_xyz->Prepare(moms,anti);
  for (int hel1(0);hel1<2;hel1++) {
    for (int hel2(0);hel2<2;hel2++) {
      amp = p_xyz->X(m_bar,hel1,pS,m_nonbar,hel2,m_c,0.);
      vector<pair<int,int> > spins;
      spins.push_back(make_pair(0,0));
      spins.push_back(make_pair(m_bar,hel1));
      spins.push_back(make_pair(m_nonbar,hel2));
      Insert(amp,spins);
    }
  }
}


////////////////////////////////////////////////////////////////////////////////
////////////////////////////////////////////////////////////////////////////////
// 0 -> SSV
////////////////////////////////////////////////////////////////////////////////
////////////////////////////////////////////////////////////////////////////////

SSV::SSV(Flavour *fl,int *i,bool *out) :
  Partial_Amplitude_Base(fl,3,i,out)
{
  AssertIn(1);
  AssertSpins(0,0,2);
}

void SSV::operator()(const Vec4D * moms,const bool anti)
{
  Vec4D pS1(p_i[0]==0?moms[p_i[0]]:-moms[p_i[0]]);
  Vec4D pS2(p_i[1]==0?moms[p_i[1]]:-moms[p_i[1]]);
  Vec4D pV(p_i[2]==0?moms[p_i[2]]:-moms[p_i[2]]);
  Flavour flV(p_flavs[p_i[2]]);
  Polarization_Vector eps(pV,sqr(flV.HadMass()),flV.IsAnti()^anti,p_out[2]);
  int npol=IsZero(flV.HadMass())?2:3;
  for (int Vpol(0);Vpol<npol;Vpol++) Insert(eps[Vpol]*(pS1-pS2),Vpol);
}


////////////////////////////////////////////////////////////////////////////////
////////////////////////////////////////////////////////////////////////////////
// 0 -> SVV
////////////////////////////////////////////////////////////////////////////////
////////////////////////////////////////////////////////////////////////////////

SVV::SVV(Flavour *fl,int *i,bool *out) :
  Partial_Amplitude_Base(fl,3,i,out)
{
  AssertIn(1);
  AssertSpins(0,2,2);
}

void SVV::operator()(const Vec4D * moms,const bool anti)
{
  Vec4D pV1(p_i[1]==0?moms[p_i[1]]:-moms[p_i[1]]);
  Vec4D pV2(p_i[2]==0?moms[p_i[2]]:-moms[p_i[2]]);
  Flavour flV1(p_flavs[p_i[1]]);
  Flavour flV2(p_flavs[p_i[2]]);
  Polarization_Vector eps1(pV1,sqr(flV1.HadMass()),flV1.IsAnti()^anti,p_out[1]);
  int npol1=IsZero(flV1.HadMass())?2:3;
  Polarization_Vector eps2(pV2,sqr(flV2.HadMass()),flV2.IsAnti()^anti,p_out[2]);
  int npol2=IsZero(flV2.HadMass())?2:3;
  for (int V1pol(0);V1pol<npol1;V1pol++) {
    for (int V2pol(0);V2pol<npol2;V2pol++) {
      vector<pair<int,int> > spins;
      spins.push_back(make_pair(0,0));
      spins.push_back(make_pair(1,V1pol));
      spins.push_back(make_pair(2,V2pol));
      Insert((eps1[V1pol]*eps2[V2pol])*(pV1*pV2)-
	     (eps1[V1pol]*pV2)*(pV1*eps2[V2pol]),spins);
    }
  }
}


////////////////////////////////////////////////////////////////////////////////
////////////////////////////////////////////////////////////////////////////////
// 0 -> VFF
////////////////////////////////////////////////////////////////////////////////
////////////////////////////////////////////////////////////////////////////////

VFF::VFF(Flavour *fl,int *i,bool *out,const Complex cL,const Complex cR) :
  Partial_Amplitude_Base(fl,3,i,out), m_cL(cL), m_cR(cR), p_xyz(NULL)
{
  AssertIn(1);
  AssertSpins(2,1,1);

  Flavour f1(p_flavs[p_i[1]]), f2(p_flavs[p_i[2]]);
  if ((f1.IsAnti()^p_out[p_i[1]])&&(f2.IsAnti()==p_out[p_i[2]])) {
    m_bar=1;
    m_nonbar=2;
  }
  else if ((f2.IsAnti()^p_out[p_i[2]])&&(f1.IsAnti()==p_out[p_i[1]])) {
    m_bar=2;
    m_nonbar=1;
  }
  else {
    msg_Error()<<METHOD<<": you provided an impossible combination of incoming/"
               <<"outgoing particles/antiparticles: "<<p_flavs[p_i[0]]
               <<" + "<<p_flavs[p_i[1]]<<" + "<<p_flavs[p_i[2]]<<endl;
  }
  p_xyz = new XYZFunc(3,fl,1,p_i);
}

void VFF::operator()(const Vec4D * moms,const bool anti)
{
  Vec4D pV(p_i[0]==0?moms[p_i[0]]:-moms[p_i[0]]);
  Flavour flV(p_flavs[p_i[0]]);
  Complex amp(0.,0.);
  p_xyz->Prepare(moms,anti);
  Polarization_Vector eps(pV,sqr(flV.HadMass()),flV.IsAnti()^anti,p_out[0]);
  int npol=IsZero(flV.HadMass())?2:3;
  for (int Vpol(0);Vpol<npol;Vpol++) {
    for (int hel1(0);hel1<2;hel1++) {
      for (int hel2(0);hel2<2;hel2++) {
	amp = p_xyz->X(m_bar,hel1,eps[Vpol],m_nonbar,hel2,m_cL,m_cR);
	vector<pair<int,int> > spins;
	spins.push_back(make_pair(0,Vpol));
	spins.push_back(make_pair(m_bar,hel1));
	spins.push_back(make_pair(m_nonbar,hel2));
	Insert(amp,spins);
      }
    }
  }
}


////////////////////////////////////////////////////////////////////////////////
////////////////////////////////////////////////////////////////////////////////
// 0 -> VVV
////////////////////////////////////////////////////////////////////////////////
////////////////////////////////////////////////////////////////////////////////

VVV::VVV(Flavour *fl,int *i,bool *out) :
  Partial_Amplitude_Base(fl,3,i,out)
{
  AssertIn(1);
  AssertSpins(2,2,2);
}

void VVV::operator()(const Vec4D * moms,const bool anti)
{
  Vec4D pV0(p_i[0]==0?moms[p_i[0]]:-moms[p_i[0]]);
  Vec4D pV1(p_i[1]==0?moms[p_i[1]]:-moms[p_i[1]]);
  Vec4D pV2(p_i[2]==0?moms[p_i[2]]:-moms[p_i[2]]);
  Flavour flV0(p_flavs[p_i[0]]);
  Flavour flV1(p_flavs[p_i[1]]);
  Flavour flV2(p_flavs[p_i[2]]);
  Polarization_Vector eps0(pV0,sqr(flV0.HadMass()),flV0.IsAnti()^anti,p_out[0]);
  Polarization_Vector eps1(pV1,sqr(flV1.HadMass()),flV1.IsAnti()^anti,p_out[1]);
  Polarization_Vector eps2(pV2,sqr(flV2.HadMass()),flV2.IsAnti()^anti,p_out[2]);
  int npol0=IsZero(flV0.HadMass())?2:3;
  int npol1=IsZero(flV1.HadMass())?2:3;
  int npol2=IsZero(flV2.HadMass())?2:3;
  Complex amp;
  for (int V0pol(0);V0pol<npol0;V0pol++) {
    for (int V1pol(0);V1pol<npol1;V1pol++) {
      for (int V2pol(0);V2pol<npol2;V2pol++) {
	vector<pair<int,int> > spins;
	spins.push_back(make_pair(0,V0pol));
	spins.push_back(make_pair(1,V1pol));
	spins.push_back(make_pair(2,V2pol));
        amp=(eps0[V0pol]*eps1[V1pol])*((pV0-pV1)*eps2[V2pol])+
            (eps1[V1pol]*eps2[V2pol])*((pV1-pV2)*eps0[V0pol])+
            (eps0[V0pol]*eps2[V2pol])*((pV2-pV0)*eps1[V1pol]);
	Insert(amp,spins);
      }
    }
  }
}


////////////////////////////////////////////////////////////////////////////////
////////////////////////////////////////////////////////////////////////////////
// 0 -> TSS
////////////////////////////////////////////////////////////////////////////////
////////////////////////////////////////////////////////////////////////////////

TSS::TSS(Flavour *fl,int *i,bool *out) :
  Partial_Amplitude_Base(fl,3,i,out)
{
  AssertIn(1);
  AssertSpins(4,0,0);
}

void TSS::operator()(const Vec4D * moms,const bool anti)
{
  Vec4D pS1(p_i[1]==0?moms[p_i[1]]:-moms[p_i[1]]);
  Vec4D pS2(p_i[2]==0?moms[p_i[2]]:-moms[p_i[2]]);
  Vec4D pT(p_i[0]==0?moms[p_i[0]]:-moms[p_i[0]]);
  Flavour flT(p_flavs[p_i[0]]);
  Polarization_Tensor eps(pT,sqr(flT.HadMass()),flT.IsAnti()^anti,p_out[0]);
  if(IsZero(flT.HadMass()))
    THROW(fatal_error, "Zero mass tensors not implemented yet.");
  for (int Tpol(0);Tpol<5;Tpol++)
    Insert(pS1*(eps[Tpol]*pS2),Tpol);
  /*!
    \f[
    {\cal M} = g\epsilon_{\mu\nu}p_1\mu p_2^\nu
    \f]
  */
}



////////////////////////////////////////////////////////////////////////////////
////////////////////////////////////////////////////////////////////////////////
// 0 -> TVS
////////////////////////////////////////////////////////////////////////////////
////////////////////////////////////////////////////////////////////////////////

TVS::TVS(Flavour *fl,int *i,bool *out) :
  Partial_Amplitude_Base(fl,3,i,out)
{
  AssertIn(1);
  AssertSpins(4,2,0);
}

void TVS::operator()(const Vec4D * moms,const bool anti)
{
  Vec4D pT(p_i[0]==0?moms[p_i[0]]:-moms[p_i[0]]);
  Vec4D pV(p_i[1]==0?moms[p_i[1]]:-moms[p_i[1]]);
  Vec4D pS(p_i[2]==0?moms[p_i[2]]:-moms[p_i[2]]);
  Flavour flT(p_flavs[p_i[0]]);
  Polarization_Tensor epsT(pT,sqr(flT.HadMass()),flT.IsAnti()^anti,p_out[0]);
  Flavour flV(p_flavs[p_i[1]]);
  Polarization_Vector epsV(pV,sqr(flV.HadMass()),flV.IsAnti()^anti,p_out[1]);
  int npolV=IsZero(flV.HadMass())?2:3;
  if(IsZero(flT.HadMass()))
    THROW(fatal_error, "Zero mass tensors not implemented yet.");
  for (int Tpol(0);Tpol<5;Tpol++) {
    for (int Vpol(0);Vpol<npolV;Vpol++) {
      vector<pair<int,int> > spins;
      spins.push_back(make_pair(0,Tpol));
      spins.push_back(make_pair(1,Vpol));
      spins.push_back(make_pair(2,0));
      Insert((epsT[Tpol]*pS)*cross(pV,epsV[Vpol],pS),spins);
    }
  }
  /*!
    \f[
    {\cal M} = g\epsilon_T^{\mu\nu}p_{P,\mu}\epsilon_{\nu\rho\sigma\kappa}
    p_V^\nu\epsilon_V^\rho p_P^\kappa
    \f]
  */
}


////////////////////////////////////////////////////////////////////////////////
////////////////////////////////////////////////////////////////////////////////
// 0 -> TVV
////////////////////////////////////////////////////////////////////////////////
////////////////////////////////////////////////////////////////////////////////

TVV::TVV(Flavour *fl,int *i,bool *out) :
  Partial_Amplitude_Base(fl,3,i,out)
{
  AssertIn(1);
  AssertSpins(4,2,2);
}

void TVV::operator()(const Vec4D * moms,const bool anti)
{
  Vec4D pT(p_i[0]==0?moms[p_i[0]]:-moms[p_i[0]]);
  Vec4D pV1(p_i[1]==0?moms[p_i[1]]:-moms[p_i[1]]);
  Vec4D pV2(p_i[2]==0?moms[p_i[2]]:-moms[p_i[2]]);
  Flavour flT(p_flavs[p_i[0]]);
  Polarization_Tensor epsT(pT,sqr(flT.HadMass()),flT.IsAnti()^anti,p_out[0]);
  Flavour flV1(p_flavs[p_i[1]]);
  Polarization_Vector epsV1(pV1,sqr(flV1.HadMass()),flV1.IsAnti()^anti,p_out[1]);
  int npolV1=IsZero(flV1.HadMass())?2:3;
  Flavour flV2(p_flavs[p_i[2]]);
  Polarization_Vector epsV2(pV2,sqr(flV2.HadMass()),flV2.IsAnti()^anti,p_out[2]);
  int npolV2=IsZero(flV2.HadMass())?2:3;
  if(IsZero(flT.HadMass()))
    THROW(fatal_error, "Zero mass tensors not implemented yet.");
  for (int T(0);T<5;T++) {
    for (int V1(0);V1<npolV1;V1++) {
      for (int V2(0);V2<npolV2;V2++) {
        vector<pair<int,int> > spins;
        spins.push_back(make_pair(0,T));
        spins.push_back(make_pair(1,V1));
        spins.push_back(make_pair(2,V2));
        Complex p1p2=pV1*pV2;
        Complex amp=((epsT[T]*pV1)*pV2-(epsT[T]*pV2)*pV1)*(epsV1[V1]*epsV2[V2])-
          ((epsT[T]*epsV1[V1])*pV2-(epsT[T]*pV2)*epsV1[V1])*(pV1*epsV2[V2])-
          ((epsT[T]*pV1)*epsV2[V2]-(epsT[T]*epsV2[V2])*pV1)*(pV2*epsV1[V1])-
          ((epsT[T]*epsV1[V1])*epsV2[V2]-(epsT[T]*epsV2[V2])*epsV1[V1])*p1p2;
        Insert(amp,spins);
      }
    }
  }
  /*!
    \f[
    {\cal M} = g\epsilon_T^{\mu\nu}[
    \left(p_\mu\epsilon_\rho-p_\rho\epsilon_\mu\right)
    \left(p'_\nu{\epsilon'}^\rho-{p'}^\rho\epsilon'_\nu\right)
    -\left(\mu\leftrightarrow\nu\right)]
    \f]
  */
}

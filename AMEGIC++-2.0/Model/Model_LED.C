#include "Model_LED.H"
#include "Run_Parameter.H"
#include "Running_AlphaS.H"
#include "Running_AlphaQED.H"
#include "MathTools.H"

using namespace AMEGIC;
using namespace AORGTOOLS;

void Model_LED::Init()
{
  moqcd.Init();
  moew.Init();
  moqcdgrav.Init();
  moewgrav.Init();
  
}

void Model_LED::c_FFV(Single_Vertex* v,int& vanz)
{
  moqcd.c_FFV(v,vanz);
  moew.c_FFV(v,vanz);
}

void Model_LED::c_VVV(Single_Vertex* v,int& vanz)
{
  moqcd.c_VVV(v,vanz);
  moew.c_VVV(v,vanz);
}
void Model_LED::c_VVVV(Single_Vertex* v,int& vanz)
{
  moqcd.c_VVVV(v,vanz);
  moew.c_VVVV(v,vanz);
}

void Model_LED::c_FFS(Single_Vertex* v,int& vanz) {moew.c_FFS(v,vanz);}
void Model_LED::c_VVS(Single_Vertex* v,int& vanz) {moew.c_VVS(v,vanz);}
void Model_LED::c_SSS(Single_Vertex* v,int& vanz) {moew.c_SSS(v,vanz);}

void Model_LED::c_FFT(Single_Vertex* v,int& vanz)
{
  moewgrav.c_FFT(v,vanz);
}
void Model_LED::c_VVT(Single_Vertex* v,int& vanz)
{
  moqcdgrav.c_VVT(v,vanz);
  moewgrav.c_VVT(v,vanz);
}
void Model_LED::c_SST(Single_Vertex* v,int& vanz)
{
  moewgrav.c_SST(v,vanz);
}
void Model_LED::c_VVVT(Single_Vertex* v,int& vanz)
{
  moqcdgrav.c_VVVT(v,vanz);
  moewgrav.c_VVVT(v,vanz);
}
void Model_LED::c_SSST(Single_Vertex* v,int& vanz)
{
  moewgrav.c_SSST(v,vanz);
}
void Model_LED::c_FFVT(Single_Vertex* v,int& vanz)
{
  moqcdgrav.c_FFVT(v,vanz);
  moewgrav.c_FFVT(v,vanz);
}


inline double Model_LED::Aqcd(double t)    {return APHYTOOLS::as->AlphaS(t);} // switch dependent running
inline double Model_LED::Aqcd()            {return APHYTOOLS::as->AsFixed();} // alpha_S _eff (read in)
inline double Model_LED::Aqed(double t)    {return APHYTOOLS::aqed->Aqed(t);}
inline double Model_LED::Aqed()            {return APHYTOOLS::aqed->AqedFixed();}

/*
inline double Model_SM::Aqcd(double t) {return (*APHYTOOLS::as)(t);}
inline double Model_SM::Aqcd()         {return (*APHYTOOLS::as)(sqr(rpa.gen.Ecms()));}
inline double Model_SM::Aqed(double t) {return (*APHYTOOLS::aqed)(t);}
inline double Model_SM::Aqed()         {return (*APHYTOOLS::aqed)(sqr(rpa.gen.Ecms()));}
*/



















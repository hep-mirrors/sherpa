#include "Model_THDM.H"
#include "Run_Parameter.H"
#include "Running_AlphaS.H"
#include "Running_AlphaQED.H"

using namespace AMEGIC;
using namespace AMATOOLS;

void Model_THDM::Init()
{
  mosm.Init();
  moHiggs.Init(&mosm,isa);
}

void Model_THDM::c_FFV(Single_Vertex* v,int& vanz)
{
  mosm.c_FFV(v,vanz);
}

void Model_THDM::c_VVV(Single_Vertex* v,int& vanz)
{
  mosm.c_VVV(v,vanz);
}

void Model_THDM::c_FFS(Single_Vertex* v,int& vanz)
{
  moHiggs.c_FFS(v,vanz);
}

void Model_THDM::c_VVS(Single_Vertex* v,int& vanz) 
{
  moHiggs.c_VVS(v,vanz);
}

void Model_THDM::c_SSS(Single_Vertex* v,int& vanz) 
{
  moHiggs.c_SSS(v,vanz);
}


void Model_THDM::c_SSV(Single_Vertex* v,int& vanz)
{ 
  moHiggs.c_SSV(v,vanz);
}

void Model_THDM::c_SSVV(Single_Vertex* v,int& vanz)
{ 
  moHiggs.c_SSVV(v,vanz);
}

void Model_THDM::c_VVVV(Single_Vertex* v,int& vanz)
{
  mosm.c_VVVV(v,vanz);
}

inline double Model_THDM::Aqcd(double t)    {return APHYTOOLS::as->AlphaS(t);} // switch dependent running
inline double Model_THDM::Aqcd()            {return APHYTOOLS::as->AsFixed();} // alpha_S _eff (read in)
inline double Model_THDM::Aqed(double t)    {return APHYTOOLS::aqed->Aqed(t);}
inline double Model_THDM::Aqed()            {return APHYTOOLS::aqed->AqedFixed();}











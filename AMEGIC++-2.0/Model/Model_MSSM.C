#include "Model_MSSM.H"
#include "Run_Parameter.H"
#include "Running_AlphaS.H"
#include "Running_AlphaQED.H"

using namespace AMEGIC;
using namespace AMATOOLS;

void Model_MSSM::Init()
{
  mosm.Init();
  moHiggs.Init(&mosm,isa);
  moch.Init(&moHiggs,isa);
  mosL.Init(&moHiggs,isa);
  mosLChi.Init(&mosm,&moch,&mosL);  
  mosq.Init(&moHiggs,isa);    
  mosqChi.Init(&mosm,&moch,&mosq);  
    
}

void Model_MSSM::c_FFV(Single_Vertex* v,int& vanz)
{
  mosm.c_FFV(v,vanz);
  moch.c_FFV(v,vanz);
  mosq.c_FFV(v,vanz);
}

void Model_MSSM::c_VVV(Single_Vertex* v,int& vanz)
{
  mosm.c_VVV(v,vanz);
}

void Model_MSSM::c_FFS(Single_Vertex* v,int& vanz)
{
  mosm.c_FFS(v,vanz);
  moHiggs.c_FFS(v,vanz);
  moch.c_FFS(v,vanz);
  mosLChi.c_FFS(v,vanz);
  mosq.c_FFS(v,vanz);
  mosqChi.c_FFS(v,vanz);
  
}

void Model_MSSM::c_VVS(Single_Vertex* v,int& vanz) 
{
  mosm.c_VVS(v,vanz);
  moHiggs.c_VVS(v,vanz);
}


void Model_MSSM::c_SSS(Single_Vertex* v,int& vanz) 
{
  mosm.c_SSS(v,vanz);
  moHiggs.c_SSS(v,vanz);
  mosL.c_SSS(v,vanz);
  mosq.c_SSS(v,vanz);
}


void Model_MSSM::c_SSV(Single_Vertex* v,int& vanz)
{ 
  moHiggs.c_SSV(v,vanz);
  mosL.c_SSV(v,vanz);
  mosq.c_SSV(v,vanz);
}

inline double Model_MSSM::Aqcd(double t)    {return APHYTOOLS::as->AlphaS(t);} // switch dependent running
inline double Model_MSSM::Aqcd()            {return APHYTOOLS::as->AsFixed();} // alpha_S _eff (read in)
inline double Model_MSSM::Aqed(double t)    {return APHYTOOLS::aqed->Aqed(t);}
inline double Model_MSSM::Aqed()            {return APHYTOOLS::aqed->AqedFixed();}
/*
inline double Model_MSSM::Aqcd(double t)    {return (*APHYTOOLS::as)(t);}
inline double Model_MSSM::Aqcd()            {return (*APHYTOOLS::as)(sqr(AORGTOOLS::rpa.gen.Ecms()));}
inline double Model_MSSM::Aqed(double t)    {return (*APHYTOOLS::aqed)(t);}
inline double Model_MSSM::Aqed()            {return (*APHYTOOLS::aqed)(sqr(AORGTOOLS::rpa.gen.Ecms()));}
*/
inline double Model_MSSM::SinTW()           {return mosm.SinTW();}
inline double Model_MSSM::CosTW()           {return mosm.CosTW();}
inline double Model_MSSM::TanB()            {return moHiggs.TanB();}
inline double Model_MSSM::CosA()            {return moHiggs.CosA();}
inline double Model_MSSM::SinA()            {return moHiggs.SinA();}










#include <stdio.h>
#include "Model_SM.H"
#include "Run_Parameter.H"
#include "Running_AlphaS.H"
#include "Running_AlphaQED.H"
#include "MathTools.H"

using namespace AMEGIC;
using namespace AORGTOOLS;
using namespace AMATOOLS;

void Model_SM::Init()
{
  moqcd.Init();
  moew.Init();
}

void Model_SM::c_FFV(Single_Vertex* v,int& vanz)
{
  moqcd.c_FFV(v,vanz);
  moew.c_FFV(v,vanz);
}

void Model_SM::c_VVV(Single_Vertex* v,int& vanz)
{
  moqcd.c_VVV(v,vanz);
  moew.c_VVV(v,vanz);
}
void Model_SM::c_VVVV(Single_Vertex* v,int& vanz)
{
  moqcd.c_VVVV(v,vanz);
  moew.c_VVVV(v,vanz);
}

void Model_SM::c_FFS(Single_Vertex* v,int& vanz) {moew.c_FFS(v,vanz);}
void Model_SM::c_VVS(Single_Vertex* v,int& vanz) {moew.c_VVS(v,vanz);}
void Model_SM::c_SSS(Single_Vertex* v,int& vanz) {moew.c_SSS(v,vanz);}

inline double Model_SM::Aqcd(double t)    {return APHYTOOLS::as->AlphaS(t);} // switch dependent running
inline double Model_SM::Aqcd()            {return APHYTOOLS::as->AsFixed();} // alpha_S _eff (read in)
inline double Model_SM::Aqed(double t)    {return APHYTOOLS::aqed->Aqed(t);}
inline double Model_SM::Aqed()            {return APHYTOOLS::aqed->AqedFixed();}

/*
inline double Model_SM::Aqcd(double t) {return (*APHYTOOLS::as)(t);}
inline double Model_SM::Aqcd()         {return (*APHYTOOLS::as)(sqr(rpa.gen.Ecms()));}
inline double Model_SM::Aqed(double t) {return (*APHYTOOLS::aqed)(t);}
inline double Model_SM::Aqed()         {return (*APHYTOOLS::aqed)(sqr(rpa.gen.Ecms()));}
*/

#include "Interaction_Model_EHC_S.H"
#include "MathTools.H"
#include "Message.H"
#include "Run_Parameter.H"
#include <stdio.h>


using namespace MODEL;
using namespace ATOOLS;
using namespace std;

Interaction_Model_EHC_S::Interaction_Model_EHC_S(MODEL::Model_Base * _model,
						 std::string _cplscheme,std::string _yukscheme) :
  Interaction_Model_Base("",_model,_cplscheme,_yukscheme)
{ 
 



void Interaction_Model_EHC_S::c_VVS(std::vector<Single_Vertex>& vertex,int& vanz)
{
}
 

void Interaction_Model_EHC_S::c_VVVV(std::vector<Single_Vertex>& vertex,int& vanz)
{
}

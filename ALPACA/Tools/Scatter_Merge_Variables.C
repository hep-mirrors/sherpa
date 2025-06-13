#include "ALPACA/Tools/Scatter_Merge_Variables.H"


using namespace ALPACA;

Scatter_Merge_Variables::Scatter_Merge_Variables():
  p_partons(0.)
{};

Scatter_Merge_Variables::Scatter_Merge_Variables(double parton):
  p_partons(parton)
{};


Scatter_Merge_Variables::~Scatter_Merge_Variables() {
};
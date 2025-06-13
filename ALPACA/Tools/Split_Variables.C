#include "ALPACA/Tools/Split_Variables.H"


using namespace ALPACA;

Split_Variables::Split_Variables():
  p_part_in(nullptr), p_part_recoil(nullptr)
{};

Split_Variables::Split_Variables(std::shared_ptr<Parton> ptr_part_in):
  p_part_in(ptr_part_in), p_part_recoil(nullptr)
{};


Split_Variables::~Split_Variables() {
};
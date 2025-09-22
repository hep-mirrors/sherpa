#include "EXTRA_XS/One2Three/H_to_bbg_NLO_Complete.H"

#include "EXTRA_XS/One2Three/H_to_bbg_Real.H"
#include "EXTRA_XS/One2Three/H_to_bb_Virtual.H"

#include "ATOOLS/Org/Message.H"

using namespace EXTRAXS;
using namespace ATOOLS;

H_to_bbg_NLO_Complete::H_to_bbg_NLO_Complete(
    const std::vector<Flavour>& flavs,
    const Flavour& prop,
    size_t non_prop, size_t gluon, size_t propj)
  : Spin_Amplitudes(flavs, Complex(0.0, 0.0)),
    p_real(nullptr), p_virtual(nullptr)
{/*
  try {
    p_real = new H_to_bbg_Real(flavs, prop, non_prop, gluon, propj);
    p_virtual = new H_to_bb_Virtual(flavs, prop);
  } catch (...) {
    delete p_real;
    delete p_virtual;
    throw;
  }*/
}

H_to_bbg_NLO_Complete::~H_to_bbg_NLO_Complete()
{/*
  delete p_real;
  delete p_virtual;*/
}

void H_to_bbg_NLO_Complete::Calculate(
    const Vec4D_Vector& momenta, bool anti)
{}

double H_to_bbg_NLO_Complete::GetRealContribution() const
{
  //return p_real ? p_real->GetMatrixElement() : 0.0;
  return 2.0;
}

double H_to_bbg_NLO_Complete::GetVirtualContribution() const
{
  //return p_virtual ? p_virtual->GetMatrixElement() : 0.0;
  return 2.0;
}
#include "SHRiMPS/Eikonals/Eikonal_Contributor.H"
#include "SHRiMPS/Tools/MinBias_Parameters.H"
#include "ATOOLS/Org/Message.H"
#include "ATOOLS/Math/MathTools.H"

using namespace SHRIMPS;
using namespace ATOOLS;

Eikonal_Contributor::
Eikonal_Contributor(Form_Factor * ffi,Form_Factor * ffk,
		    const Eikonal_Parameters & params) :
  p_ffi(ffi), p_ffk(ffk), m_params(params), m_grid(Grid(params.Ymax))
{ }

void Eikonal_Contributor::
PrepareGrid(const int & ff1steps, const int & ff2steps)
{
  m_grid.Initialise(ff1steps,ff2steps,p_ffi->Maximum(),p_ffk->Maximum());
}

void Eikonal_Contributor::InsertValues(const size_t & i,const size_t & j,
				       const std::vector<double> & values) {
  m_grid.InsertValues(i,j,values);
}

void Eikonal_Contributor::SetB1B2(const double & b1,const double & b2) { 
  m_b1 = b1; m_b2 = b2; 
}

double Eikonal_Contributor::
operator()(const double & b1,const double & b2,const double & y) {
  double ff1 = p_ffi->FourierTransform(b1);
  double ff2 = p_ffk->FourierTransform(b2);
  return m_grid(ff1,ff2,y);
}

double Eikonal_Contributor::operator()(const double & y) {
  double ff1 = p_ffi->FourierTransform(m_b1);
  double ff2 = p_ffk->FourierTransform(m_b2);
  return m_grid(ff1,ff2,y);
}

bool Eikonal_Contributor::Valid(const double & y) const 
{
  if (IsNan(y)) return false;
  if (m_b1<0. || m_b1>=m_params.bmax || 
      m_b2<0. || m_b2>=m_params.bmax)   return false;
  if (p_ffi->FourierTransform(m_b1)<0. ||
      p_ffk->FourierTransform(m_b2)<0.) return false;
  return true;
}





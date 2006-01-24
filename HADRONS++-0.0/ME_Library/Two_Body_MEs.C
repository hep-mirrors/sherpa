#include "Two_Body_MEs.H"

using namespace HADRONS;
using namespace ATOOLS;
using namespace std;

P_2Gamma::P_2Gamma(int _nout,Flavour * _flavs) :
  HD_ME_Base(_nout,_flavs)
{ 
  m_metype = string("P_2Gamma");
}

void   P_2Gamma::operator()( 
    const ATOOLS::Vec4D  * _p, 
    std::vector<Complex> * _ampls_tensor, 
    std::vector<std::pair<int,int> > * _indices,
    int                    k0_n )
{
  _ampls_tensor->clear();
  _ampls_tensor->push_back( p_masses2[0]/2. );
  _indices->clear();
}

double P_2Gamma::operator()(const Vec4D * moms)
{
  return sqr(p_masses2[0])/4.;
}



#include "HELICITIES/Main/Full_Amplitude_Base.H"
#include "HELICITIES/Main/Partial_Amplitude_Base.H"

using namespace HELICITIES;
using namespace ATOOLS;
using namespace std;

Full_Amplitude_Base::Full_Amplitude_Base(Flavour* flavs,size_t size) :
  Spin_Structure<std::vector<Complex> >(flavs,size,std::vector<Complex>(0)), 
  p_flavs(flavs), p_colormatrix(NULL)
{}

Full_Amplitude_Base::~Full_Amplitude_Base() {
  for (size_t i=0;i<m_amplitudes.size();i++) {
    if (m_amplitudes[i]!=NULL) {
      delete m_amplitudes[i]; m_amplitudes[i]=NULL;
    }
  }
}

double Full_Amplitude_Base::SummedSquared(const Vec4D * moms,bool anti) {
  double result(0.0);
  for (size_t i(0);i<m_amplitudes.size();++i) {
    (*m_amplitudes[i])(moms,anti);
    for (size_t j(0);j<=i;++j) {
      for (size_t hels=0;hels<size();hels++) {
	result += 
	  abs(m_amplitudes[i]->Get(hels)*
	      conj(m_amplitudes[j]->Get(hels))*
	      (*p_colormatrix)[i][j]);
//   	std::cout<<METHOD<<": (ij) = ("<<i<<j<<"), hels = "<<hels
//   		 <<" --> "<<m_amplitudes[i]->Get(hels)
//   		 <<" * "<<conj(m_amplitudes[j]->Get(hels))
//   		 <<" with "<<(*p_colormatrix)[i][j]<<std::endl
//   		 <<"             --->  Result = "<<result<<std::endl;
	if (i!=j)
	  result += 
	    abs(m_amplitudes[j]->Get(hels)*
		conj(m_amplitudes[i]->Get(hels))*
		(*p_colormatrix)[j][i]);
      }
    }
  }
  return result;
}

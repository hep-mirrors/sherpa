#include "AddOns/Analysis/Observables/Momentum_Tensor.H"
#include <vector>

using namespace ANALYSIS;
using namespace ATOOLS;


Momentum_Tensor::Momentum_Tensor(double power) : m_power(power)
{
}

void Momentum_Tensor::Calculate(const ATOOLS::Particle_List & pl)
{
  std::vector<double> scale;
  double denom=0.;
  for (Particle_List::const_iterator pit=pl.begin();pit!=pl.end();++pit) {
    if (m_power==2.) {
      // quadratic (normed) momentum tensor (i.e. Sphericity, Aplanarity, ...)
      denom+=Vec3D((*pit)->Momentum()).Sqr();
      scale.push_back(1.);
    }
    else {
      // e.g. m_power==1 linear moment tensor (i.e. C-parameter, D-parameter, ...)
      double p=Vec3D((*pit)->Momentum()).Abs();
      denom+=pow(p,m_power);
      scale.push_back(pow(p,(m_power-2.)));
    }
  }

  // now calculate the tensor
  double s[3][3]={{0,0,0}, {0,0,0}, {0,0,0}};
  for (int i=0; i<3; i++)
    for (int j=0; j<=i; j++) {
      for (size_t k=0; k<pl.size(); ++k)  
	s[i][j]+=scale[k]*pl[k]->Momentum()[i+1]*pl[k]->Momentum()[j+1];
      s[i][j]/=denom; 
      s[j][i]=s[i][j];
    }

  Evaluate(Matrix<3>(s));
}

void Momentum_Tensor::Evaluate(const ATOOLS::Matrix<3> & m)
{
  Matrix<3> v;
  double r[3];
  m.DiagonalizeSort(r,v);

  ///  std::cout<<"R: "<<r[0]<<";"<<r[1]<<";"<<r[2]<<std::endl;
  
  // sphericity and co. are usually stored in descending order
  // Matrix returns them ascending!
  for (int i=2; i>=0; i--) {
    m_eigenvalues[2-i]=r[i];
    m_eigenvectors[2-i]=Vec3D(v[0][i],v[1][i],v[2][i]);
  }
}

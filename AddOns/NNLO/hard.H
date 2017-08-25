#ifndef NNLO__hard_H
#define NNLO__hard_H

#include "param.H"
#ifdef USING__VV
#include "VVhard.H"
#endif

namespace SHNNLO {

#ifndef USING__VV
  std::complex<double> spa(int i,int j);
  std::complex<double> spb(int i,int j);
  double sij(int i,int j);
#endif

  double hf0qq(double mu, double Q);
  double hf1qq(double mu, double Q);
  double hf2qq(double mu, double Q);

  void updatecH1qq(int i, int vm);
  void updatecH2qq(int i, int vm);

  double hf0gg(double mu, double Q);
  double hf1gg(double mu, double Q);
  double hf2gg(double mu, double Q);

  double hf0tt(double mu, double Q, double mt);
  double hf1tt(double mu, double Q, double mt);
  double hf2tt(double mu, double Q, double mt);
  double hf3tt(double mu, double Q, double mt);

}

#endif

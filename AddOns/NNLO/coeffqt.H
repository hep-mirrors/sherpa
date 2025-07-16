#include "param.H"

namespace SHNNLO {

  double PDF(int i, double x, double mu);

  double Cqq0qiqi(int i1, int i2, double x1, double x2, double z1, double z2, double qT, double mu, double Q);
  double Cqq1qg(int i1, int i2, double x1, double x2, double z1, double z2, double qT, double mu, double Q);
  double Cqq1qiqi(int i1, int i2, double x1, double x2, double z1, double z2, double qT, double mu, double Q);
  double Cqq2qiqj(int i1, int i2, double x1, double x2, double z1, double z2, double qT, double mu, double Q);
  double Cqq2qiqbi(int i1, int i2, double x1, double x2, double z1, double z2, double qT, double mu, double Q);
  double Cqq2qg(int i1, int i2, double x1, double x2, double z1, double z2, double qT, double mu, double Q);
  double Cqq2gg(int i1, int i2, double x1, double x2, double z1, double z2, double qT, double mu, double Q);
  double Cqq2qiqi(int i1, int i2, double x1, double x2, double z1, double z2, double qT, double mu, double Q);

  double Cgg0gg(double x1, double x2, double z1, double z2, double qT, double mu, double Q);
  double Cgg1gq(double x1, double x2, double z1, double z2, double qT, double mu, double Q);
  double Cgg1gg(double x1, double x2, double z1, double z2, double qT, double mu, double Q);
  double Cgg2gq(double x1, double x2, double z1, double z2, double qT, double mu, double Q);
  double Cgg2qq(double x1, double x2, double z1, double z2, double qT, double mu, double Q);
  double Cgg2gg(double x1, double x2, double z1, double z2, double qT, double mu, double Q);

}

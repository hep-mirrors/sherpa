#include "ATOOLS/Math/Tensor.H"
#include "ATOOLS/Math/MyComplex.H"
#include "ATOOLS/Math/MathTools.H"
#include "ATOOLS/Org/Exception.H"
#include "ATOOLS/Org/MyStrStream.H"

using namespace std;
using namespace ATOOLS;

// IsEqual methods

bool ATOOLS::IsEqual(const Lorentz_Ten2D& t1, const Lorentz_Ten2D& t2, const double crit)
{
  double max(0.);
  for(unsigned short int i=0;i<4;++i) {
    for(unsigned short int j=0;j<4;++j) {
      double absij = 0.5*dabs(t1.at(i,j)+t2.at(i,j));
      if (absij > max) max = absij;
    }
  }
  if (IsZero(max))  return true;
  for(unsigned short int i=0;i<4;++i) 
    for(unsigned short int j=0;j<4;++j) 
      if (dabs((t1.at(i,j)-t2.at(i,j))/max)>crit &&
         !(dabs(t1.at(i,j))<=crit && 
           dabs(t2.at(i,j))<=crit)) return false;
  return true;
}

bool ATOOLS::IsEqual(const Lorentz_Ten2C& t1, const Lorentz_Ten2C& t2, const double crit)
{
  double max(0.);
  for(unsigned short int i=0;i<4;++i) {
    for(unsigned short int j=0;j<4;++j) {
      double absij = 0.5*cabs(t1.at(i,j)+t2.at(i,j));
      if (absij > max) max = absij;
    }
  }
  if (IsZero(max))  return true;
  for(unsigned short int i=0;i<4;++i) 
    for(unsigned short int j=0;j<4;++j) 
      if (cabs((t1.at(i,j)-t2.at(i,j))/max)>crit &&
         !(cabs(t1.at(i,j))<=crit && 
           cabs(t2.at(i,j))<=crit)) return false;
  return true;
}

bool ATOOLS::IsEqual(const Lorentz_Ten3D& t1, const Lorentz_Ten3D& t2, const double crit)
{
  double max(0.);
  for(unsigned short int i=0;i<4;++i) {
    for(unsigned short int j=0;j<4;++j) {
      for(unsigned short int k=0;k<4;++k) {
        double absij = 0.5*dabs(t1.at(i,j,k)+t2.at(i,j,k));
        if (absij > max) max = absij;
      }
    }
  }
  if (IsZero(max))  return true;
  for(unsigned short int i=0;i<4;++i) 
    for(unsigned short int j=0;j<4;++j) 
      for(unsigned short int k=0;k<4;++k)
        if (dabs((t1.at(i,j,k)-t2.at(i,j,k))/max)>crit &&
          !(dabs(t1.at(i,j,k))<=crit && 
            dabs(t2.at(i,j,k))<=crit)) return false;
  return true;
}

bool ATOOLS::IsEqual(const Lorentz_Ten3C& t1, const Lorentz_Ten3C& t2, const double crit)
{
  double max(0.);
  for(unsigned short int i=0;i<4;++i) {
    for(unsigned short int j=0;j<4;++j) {
      for(unsigned short int k=0;k<4;++k) {
        double absij = 0.5*cabs(t1.at(i,j,k)+t2.at(i,j,k));
        if (absij > max) max = absij;
      }
    }
  }
  if (IsZero(max))  return true;
  for(unsigned short int i=0;i<4;++i) 
    for(unsigned short int j=0;j<4;++j) 
      for(unsigned short int k=0;k<4;++k) 
        if (cabs((t1.at(i,j,k)-t2.at(i,j,k))/max)>crit &&
          !(cabs(t1.at(i,j,k))<=crit && 
            cabs(t2.at(i,j,k))<=crit)) return false;
  return true;
}

bool ATOOLS::IsEqual(const Lorentz_Ten4D& t1, const Lorentz_Ten4D& t2, const double crit)
{
  double max(0.);
  for(unsigned short int i=0;i<4;++i) {
    for(unsigned short int j=0;j<4;++j) {
      for(unsigned short int k=0;k<4;++k) {
        for(unsigned short int l=0;l<4;++l) {
          double absij = 0.5*dabs(t1.at(i,j,k,l)+t2.at(i,j,k,l));
          if (absij > max) max = absij;
        }
      }
    }
  }
  if (IsZero(max))  return true;
  for(unsigned short int i=0;i<4;++i) 
    for(unsigned short int j=0;j<4;++j) 
      for(unsigned short int k=0;k<4;++k)
        for(unsigned short int l=0;l<4;++l)
          if (dabs((t1.at(i,j,k,l)-t2.at(i,j,k,l))/max)>crit &&
            !(dabs(t1.at(i,j,k,l))<=crit && 
              dabs(t2.at(i,j,k,l))<=crit)) return false;
  return true;
}

bool ATOOLS::IsEqual(const Lorentz_Ten4C& t1, const Lorentz_Ten4C& t2, const double crit)
{
  double max(0.);
  for(unsigned short int i=0;i<4;++i) {
    for(unsigned short int j=0;j<4;++j) {
      for(unsigned short int k=0;k<4;++k) {
        for(unsigned short int l=0;l<4;++l) {
          double absij = 0.5*cabs(t1.at(i,j,k,l)+t2.at(i,j,k,l));
          if (absij > max) max = absij;
        }
      }
    }
  }
  if (IsZero(max))  return true;
  for(unsigned short int i=0;i<4;++i) 
    for(unsigned short int j=0;j<4;++j) 
      for(unsigned short int k=0;k<4;++k) 
        for(unsigned short int l=0;l<4;++l) 
          if (cabs((t1.at(i,j,k,l)-t2.at(i,j,k,l))/max)>crit &&
            !(cabs(t1.at(i,j,k,l))<=crit && 
              cabs(t2.at(i,j,k,l))<=crit)) return false;
  return true;
}









#include "HELICITIES/Loops/Golem95_Wrapper.H"
#include "ATOOLS/Org/CXXFLAGS_PACKAGES.H"
#include "ATOOLS/Org/Message.H"

#ifdef USING__GOLEM95

#define f90(l,n) __ ## l ## _MOD_ ## n

extern "C" {
  void f90(matrice_s, allocation_s)(int*);
  void f90(cache,allocate_cache)(int*);
  void f90(matrice_s,init_invs)();
  void f90(cache,clear_cache)();
  void f90(matrice_s,deallocation_s)();
  void f90(golem95_wrap,set_s_mat)(int*, int*, double*);
  void f90(golem95_wrap,set_ref3)(int*, int*, int*);
  void f90(golem95_wrap,set_ref4)(int*, int*, int*, int*);
  void f90(golem95_wrap,a40w)(double*,double*,double*,double*,double*,double*);
  void f90(golem95_wrap,a32w)(int*, int*,
                              double*,double*,double*,double*,double*,double*);
  void f90(golem95_wrap,b32w)(double*,double*,double*,double*,double*,double*);

  void f90(golem95_wrap,a44w)(int*, int*, int*, int*,
                              double*,double*,double*,double*,double*,double*);
  void f90(golem95_wrap,b44w)(int*, int*,
                              double*,double*,double*,double*,double*,double*);
  void f90(golem95_wrap,c44w)(double*,double*,double*,double*,double*,double*);
}

using namespace HELICITIES;
using namespace std;

namespace golem95 {

void init(int npoint)
{
  f90(matrice_s,allocation_s)(&npoint);
  
  int i(1), j(2), k(3), l(4);
  switch (npoint) {
  case 3:
    f90(golem95_wrap, set_ref3)(&i, &j, &k);
    break;
  case 4:
    f90(golem95_wrap, set_ref4)(&i, &j, &k, &l);
    break;
  default:
    cout<<"Error: npoint="<<npoint
        <<" not implemented in golem95 wrapper."<<endl;
    abort();
  }
  
  f90(cache,allocate_cache)(&npoint);

  
}

void set_s_mat(const std::vector<std::vector<double> >& smat)
{
  DEBUG_FUNC("");
  for (int i=1; i<int(smat.size()+1); ++i) {
    for (int j=1; j<int(smat[i-1].size()+1); ++j) {
      double smatij(smat[i-1][j-1]);
      f90(golem95_wrap,set_s_mat)(&i, &j, &smatij);
      msg_Debugging()<<"s_mat("<<i<<","<<j<<")="<<smatij<<".d0"<<std::endl;
    }
  }
  f90(matrice_s,init_invs)();
}


void a40(DivArrC& result)
{
  f90(golem95_wrap,a40w)(&result.IR2().real(), &result.IR2().imag(),
                         &result.IR().real(), &result.IR().imag(),
                         &result.Finite().real(), &result.Finite().imag());
}

void a32(int l1, int l2, DivArrC& result)
{
  f90(golem95_wrap,a32w)(&l1, &l2,
                         &result.IR2().real(), &result.IR2().imag(),
                         &result.IR().real(), &result.IR().imag(),
                         &result.Finite().real(), &result.Finite().imag());
}

void b32(DivArrC& result)
{
  f90(golem95_wrap,b32w)(&result.IR2().real(), &result.IR2().imag(),
                         &result.IR().real(), &result.IR().imag(),
                         &result.Finite().real(), &result.Finite().imag());
}

void a44(int l1, int l2, int l3, int l4, DivArrC& result)
{
  f90(golem95_wrap,a44w)(&l1, &l2, &l3, &l4,
                         &result.IR2().real(), &result.IR2().imag(),
                         &result.IR().real(), &result.IR().imag(),
                         &result.Finite().real(), &result.Finite().imag());
}

void b44(int l1, int l2, DivArrC& result)
{
  f90(golem95_wrap,b44w)(&l1, &l2,
                         &result.IR2().real(), &result.IR2().imag(),
                         &result.IR().real(), &result.IR().imag(),
                         &result.Finite().real(), &result.Finite().imag());
}

void c44(DivArrC& result)
{
  f90(golem95_wrap,c44w)(&result.IR2().real(), &result.IR2().imag(),
                         &result.IR().real(), &result.IR().imag(),
                         &result.Finite().real(), &result.Finite().imag());
}

void finish()
{
  f90(cache,clear_cache)();
  f90(matrice_s,deallocation_s)();
}
}

#endif

#ifndef METOOLS_Loops_GOLEM95_WRAPPER_H
#define METOOLS_Loops_GOLEM95_WRAPPER_H

#include "ATOOLS/Org/CXXFLAGS_PACKAGES.H"
#ifdef USING__GOLEM95
#include <vector>
#include "METOOLS/Loops/Divergence_Array.H"


namespace golem95 {

void init(int npoint);
void set_s_mat(const std::vector<std::vector<double> >& smat);

void a40(METOOLS::DivArrC& result);
void a32(int l1, int l2, METOOLS::DivArrC& result);
void b32(METOOLS::DivArrC& result);

void a44(int l1, int l2, int l3, int l4, METOOLS::DivArrC& result);
void b44(int l1, int l2, METOOLS::DivArrC& result);
void c44(METOOLS::DivArrC& result);

void finish();

}

#endif
#endif

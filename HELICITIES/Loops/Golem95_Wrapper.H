#ifndef HELICITIES_Loops_GOLEM95_WRAPPER_H
#define HELICITIES_Loops_GOLEM95_WRAPPER_H

#include "ATOOLS/Org/CXXFLAGS_PACKAGES.H"
#ifdef USING__GOLEM95
#include <vector>
#include "HELICITIES/Loops/Divergence_Array.H"


namespace golem95 {

void init(int npoint);
void set_s_mat(const std::vector<std::vector<double> >& smat);

void a40(HELICITIES::DivArrC& result);
void a32(int l1, int l2, HELICITIES::DivArrC& result);
void b32(HELICITIES::DivArrC& result);

void a44(int l1, int l2, int l3, int l4, HELICITIES::DivArrC& result);
void b44(int l1, int l2, HELICITIES::DivArrC& result);
void c44(HELICITIES::DivArrC& result);

void finish();

}

#endif
#endif

#ifndef SHERPA_HerwigTools_MCatNLO_Wrapper_H
#define SHERPA_HerwigTools_MCatNLO_Wrapper_H

#include <ctype.h>
#include <cstring>

extern struct {
  char qqin[50];
} vvjin_;
#define vvjin vvjin_

extern struct {
  double tmplam;
} ctmplam_;
#define ctmplam ctmplam_
  
#define mccint mccint_ 
#define mccgse mccgse_
#define mcctrm mcctrm_

extern "C" {
  void mccint(int*);
  void mccgse();
  void mcctrm();
}

void mccinit(int _print) { mccint(&_print); }

//////////////////
// Doxygen part //
//////////////////

/*!
  \file MCatNLOWrapper.H
  \brief Contains definitions of \b MC@NLO common blocks
*/

#endif

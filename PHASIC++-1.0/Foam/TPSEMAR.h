#ifndef TPSEMAR_H
#define TPSEMAR_H
#include "ROOT_DEF.h"

#include <iostream>
/////////////////////////////////////////////////////////////////////////////
//                                                                         //
//       PseudoRandom number generator RANMAR used by Foam 2.x             //
//                                                                         //
/////////////////////////////////////////////////////////////////////////////
#ifdef ROOT_DEF
#include "TROOT.h"
class TPSEMAR : public TObject {
#else
class TPSEMAR{
#endif
 private:
  int       m_modcns;
  double    m_U[98], m_C,   m_CD,   m_CM,   m_twom24;
  long int  m_i97,   m_j97, m_ntot, m_ijkl, m_ntot2;
 public:
  TPSEMAR(void);
  TPSEMAR(long ,long ,long );
  ~TPSEMAR();
  void Initialize(long ,long ,long );
  void MakeVec(int, double *);
  void Out(long int &ijkl_out,long int &ntot_out,long int &ntot2_out);
 private:
  //inline functions
  long int max(long int x, long int y){ if (x>=y) return(x);  else return(y);}
  long int MOD(long int x, long int y){ return(x % y);}
/////////////////////////////////////////////////////////////////////////////
#ifdef ROOT_DEF
  ClassDef(TPSEMAR,1) //PseudoRandom number generator (RANMAR)
#endif
};
/////////////////////////////////////////////////////////////////////////////
//                                                                         //
/////////////////////////////////////////////////////////////////////////////
#endif


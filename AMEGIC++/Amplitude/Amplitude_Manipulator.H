#ifndef AMEGIC_Amplitude_Amplitude_Manipulator_H
#define AMEGIC_Amplitude_Amplitude_Manipulator_H

#include "ATOOLS/Phys/Flavour.H"
#include "AMEGIC++/Main/Topology.H"
#include "AMEGIC++/Amplitude/Single_Amplitude.H"

namespace AMEGIC {

  class Amplitude_Manipulator {
    int N;
    ATOOLS::Flavour* fl;
    int* b;
    int dec;

    void   SetPrev(Point*);
    int    Permutation(int*,int);
    void   GetPerm(int*,Single_Amplitude*,int&);
    Point* FindNext(Point*);
    void   GetFermionLine(Point*,Point*&,Point*&);
    Point* ForwardLine(Point*);
    Point* BackwardLine(Point*);
    
    int  SetPropOrientation(Point*,Point*);
    void ForwardLineOrientation(Point*,int&);
    void BackwardLineOrientation(Point*,int&);
    
    int  SetFermionNumberFlow(Point*,Point*);
    void SetForwardFNFlow(Point*,int,int&);
    void SetBackwardFNFlow(Point*,int,int&);
    
    void SetMajoFlowForward(Point*,int);
    void SetMajoFlowBackward(Point*,int);
    
  public:
    Amplitude_Manipulator(int _no,ATOOLS::Flavour* _fl,int* _b,int d=0): 
      N(_no), fl(_fl), b(_b), dec(d){}
    void   FixSign(Single_Amplitude*);  
  };

 /*! 
    \file
    \brief this file contains the class AMEGIC::Amplitude_Manipulator
  
  */

  /*! \class Amplitude_Manipulator
      \brief is capable for fixing the relative sign of interferring Feynman graphs 
      
      A basic feature of many extensions of the SM, {\e e.g.} the MSSM, is the 
      occurrence of Majorana fermion fields. Due to their self--conjugated nature 
      diagrams with ``clashing arrows'' enter and a proper determination of the 
      Relative Sign of Interfering Feynman graphs (RSIF) is far from being trivial.
      An approach easily applicable for an automated generation of Feynman graphs has 
      been provided in \cite{Denner:1992vz}. The Feynman rules defined there are close 
      to the rules for Dirac fermions. Since for Majorana fermions the fermion number 
      flow is violated the authors define a new kind of fermion chain orientation, called 
      fermion flow. This continuous orientation forces the introduction of two analytic 
      expressions for each vertex, one for fermion flow parallel and one for fermion 
      flow anti--parallel to the conventional fermion number flow. However, the derived 
      Feynman rules ask only for a minimal number of vertices and avoid the explicit 
      occurrence of charge conjugation matrices in the vertices and propagators. In this
      framework the RSIF is determined from the permutation of the external
      fermions in the usual way. The precise determination of the RSIF according to the 
      algorithm defined in \cite{Denner:1992vz} is provided by the class 
      AMEGIC::Amplitude_Manipulator      
  */
  
  /*!
    \fn void Amplitude_Manipulator::SetPrev(Point*);
    \brief Describe SetPrev(Point*) ...;
    
  */
  /*!  
    \fn int Amplitude_Manipulator::Permutation(int*,int);
     \brief Describe Permutation(int*,int);
    */
  /*!
   \fn void Amplitude_Manipulator::GetPerm(int*,Single_Amplitude*,int&);
     \brief Describe GetPerm(int*,Single_Amplitude*,int&);
    */
  /*!
    \fn Point* Amplitude_Manipulator::FindNext(Point*);
     \brief Describe FindNext(Point*);
    */    
  /*! 
    \fn void Amplitude_Manipulator::GetFermionLine(Point*,Point*&,Point*&);
       \brief Describe GetFermionLine(Point*,Point*&,Point*&);
      */
  /*! 
   \fn Point* Amplitude_Manipulator::ForwardLine(Point*);
       \brief Describe ForwardLine(Point*);
      */
  /*! 
      \fn Point* Amplitude_Manipulator::BackwardLine(Point*);
       \brief Describe BackwardLine(Point*);
      */    
  /*! 
      \fn int  Amplitude_Manipulator::SetPropOrientation(Point*,Point*);
      \brief Describe SetPropOrientation(Point*,Point*);
  */ 
  /*! 
    \fn void Amplitude_Manipulator::ForwardLineOrientation(Point*,int&);
     \brief Describe ForwardLineOrientation(Point*,int&);
    */
  /*! 
    \fn void Amplitude_Manipulator::BackwardLineOrientation(Point*,int&);
     \brief Describe BackwardLineOrientation(Point*,int&);
    */
  /*! 
    \fn void Amplitude_Manipulator::SetFermionNumberFlow(Point*,Point*);
   \brief Describe SetFermionNumberFlow(Point*,Point*);
  */
  /*! 
    \fn void Amplitude_Manipulator::SetForwardFNFlow(Point*,int);
    \brief Describe SetForwardFNFlow(Point*,int);
  */
  /*! 
    \fn void Amplitude_Manipulator::SetBackwardFNFlow(Point*,int);
   \brief Describe SetBackwardFNFlow(Point*,int);
  */
}
#endif

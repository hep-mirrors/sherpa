#ifndef Perturbative_Interface_h
#define Perturbative_Interface_h

#include "Matrix_Element_Handler.H"
#include "Shower_Handler.H"
#include "Blob.H"

namespace SHERPA {
  class Perturbative_Interface {
  protected :
    bool   m_ini,m_fin;
    int    m_nin,m_nout,m_maxjetnumber;
    double m_weight;

    Matrix_Element_Handler * p_me;
    Shower_Handler         * p_shower;        
    APHYTOOLS::Flavour     * p_fl;
    AMATOOLS::Vec4D        * p_moms;
  public :
    Perturbative_Interface(Matrix_Element_Handler * _me,Shower_Handler * _shower) :
      p_me(_me), p_shower(_shower), 
      m_ini(_shower->ISROn()), m_fin(_shower->FSROn()), m_maxjetnumber(_me->MaxJets()), m_weight(1.),
      p_fl(NULL), p_moms(NULL) {}

    virtual ~Perturbative_Interface() 
    {
      if (p_fl)   { delete p_fl;   p_fl   = NULL; }
      if (p_moms) { delete p_moms; p_moms = NULL; }
    }
    virtual bool   DefineInitialConditions(APHYTOOLS::Blob *) = 0;
    double GetWeight() { return m_weight; }
  };
}
#endif








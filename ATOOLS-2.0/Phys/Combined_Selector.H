#ifndef Combined_Selector_h
#define Combnied_Selector_h

#include "Selector.H"


namespace APHYTOOLS {
  class Combined_Selector : public Selector_Base {
  protected:
    long int       m_count;
    std::vector<Selector_Base *> m_sels;
  public:
    Combined_Selector() {};
    Combined_Selector(int,int,Flavour *,Selector_Data *);
    ~Combined_Selector();
    void Add(Selector_Base *);
    bool Trigger(const AMATOOLS::Vec4D *);
    void BuildCuts(Cut_Data *);
    void UpdateCuts(double,double,Cut_Data *);
    void Output();
  };


/*!
  \file contains APHYTOOLS::Combined_Selector
*/

/*!
  \class Combined_Selector
  \brief A list of selectors that are applied.

  This is the class to generate lists of selectors that come into
  play jointly. All methods have the same functionality as in the
  Selector_Base from which it is derived.
*/
}

#endif

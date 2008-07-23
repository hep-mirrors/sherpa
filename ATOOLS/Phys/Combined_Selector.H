#ifndef Combined_Selector_h
#define Combnied_Selector_h

#include "Selector.H"


namespace ATOOLS {
  class Combined_Selector : public Selector_Base {
  protected: 
    long int       m_count;
    int            m_update;
    std::vector<Selector_Base *> m_sels;
    std::vector<std::pair<std::string,double> > m_osc;
  public:
    Combined_Selector() {};
    Combined_Selector(int,int,Flavour *,Selector_Data *,
		      std::string ycut="-1.");
    ~Combined_Selector();
    void Add(Selector_Base *);
    bool Trigger(const Vec4D *);
    void BuildCuts(Cut_Data *);
    int  NeedUpdate();
    void UpdateCuts(double,double,Cut_Data *);
    void Output();
    size_t Size();
    Selector_Base * operator[](size_t);
    Selector_Base * GetSelector(std::string);
    void SetProcessName(const std::string &name);
    void AddOnshellCondition(std::string,double);
  };


/*!
  \file contains Combined_Selector
*/

/*!
  \class Combined_Selector
  \brief A list of selectors that are applied.

  This is the class to generate lists of selectors that come into
  play jointly. All methods have the same functionality as in the
  Selector_Base from which it is derived.
*/

  inline size_t Combined_Selector::Size() { return m_sels.size(); }
  inline Selector_Base * Combined_Selector::operator[](size_t i) { return m_sels[i]; }


}

#endif

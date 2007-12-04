#ifndef Variable_Selector_H
#define Variable_Selector_H

#include "Selector.H"

namespace ATOOLS {

  template <typename NType> class Variable_Base;

  class Order_Base;

  class Variable_Selector : public Selector_Base {

    Variable_Base<double> *p_variable;
    Order_Base            *p_order;

    std::vector<std::pair<double,double> > m_bounds;

    std::vector<std::vector<int> >   m_sels;
    std::vector<std::vector<Vec4D> > m_moms;
    std::vector<Flavour>             m_cfl;
    std::vector<size_t>              m_nfl;

    std::string m_omode;

    bool Trigger(const Vec4D *p,size_t &l,std::vector<Vec4D> &moms,
		 const size_t &f,const size_t &n,const size_t &m);

  public:

    // constructor
    Variable_Selector(const int &nin,const int &nout,
		      Flavour *const fl,const std::string &name);

    // destructor
    ~Variable_Selector();

    // member functions
    void BuildCuts(Cut_Data *cuts);

    void SetRange(std::vector<Flavour> fl,
		  std::vector<std::pair<double,double> > &bounds);

    bool Trigger(const Vec4D *p);

  };// end of class Variable_Selector

}// end of namespace ATOOLS

#endif

#ifndef Splitting_Function_Group_h
#define Splitting_Function_Group_h

#include "Splitting_Function_Base.H"
#include <list>

namespace CS_SHOWER {
  class Splitting_Function_Group : public Splitting_Function_Base {
  protected:
    std::list<Splitting_Function_Base *>           m_splittings;
    std::list<Splitting_Function_Base *>::iterator m_splitter;
    Splitting_Function_Base *            p_selected;
  public:
    Splitting_Function_Group() : Splitting_Function_Base(), p_selected(NULL) {}

    virtual ~Splitting_Function_Group();
    virtual void Add(Splitting_Function_Base *);
    virtual void SelectOne();

    virtual double operator() (const double,const double,const double,const double);
    virtual double OverIntegrated(const double,const double,const double,const double);
    virtual double Overestimated(const double,const double);
    virtual double RejectionWeight(const double,const double,
				   const double eta=1.,const double scale=0.);
    virtual double Z();

    virtual ATOOLS::Flavour & GetFlavourA() { return p_selected->GetFlavourA(); }
    virtual ATOOLS::Flavour & GetFlavourB() { return p_selected->GetFlavourB(); }
    virtual ATOOLS::Flavour & GetFlavourC() { return p_selected->GetFlavourC(); }

    friend std::ostream& operator<<(std::ostream &,Splitting_Function_Group &);
  };
}

#endif

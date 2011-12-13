#ifndef CSSHOWER_Showers_Splitting_Function_Group_H
#define CSSHOWER_Showers_Splitting_Function_Group_H

#include "CSSHOWER++/Showers/Splitting_Function_Base.H"
#include <list>

namespace CSSHOWER {
  class Splitting_Function_Group : public Splitting_Function_Base {
  protected:
    std::vector<Splitting_Function_Base *>           m_splittings;
    std::vector<Splitting_Function_Base *>::iterator m_splitter;
    Splitting_Function_Base *            p_selected;
    std::vector<double> m_partint;
  public:
    Splitting_Function_Group(): p_selected(NULL) {}

    virtual ~Splitting_Function_Group();
    virtual void Add(Splitting_Function_Base *);
    virtual void SelectOne();

    virtual double operator() (const double,const double,const double,
			       const double,const double);
    virtual double OverIntegrated(const double,const double,const double,const double);
    virtual double Overestimated(const double,const double);
    virtual double RejectionWeight(const double,const double,
				   const double eta=1.,const double scale=0.,
				   const double Q2=0.);
    virtual double Z();
    virtual void ResetLastInt();
    virtual void ClearSpecs();

    virtual const ATOOLS::Flavour & GetFlavourA() const    { return p_selected->GetFlavourA(); }
    virtual const ATOOLS::Flavour & GetFlavourB() const    { return p_selected->GetFlavourB(); }
    virtual const ATOOLS::Flavour & GetFlavourC() const    { return p_selected->GetFlavourC(); }
    virtual const ATOOLS::Flavour & GetFlavourSpec() const { return p_selected->GetFlavourSpec(); }

    virtual int GetCol() const { return p_selected->GetCol(); }

    inline Splitting_Function_Base *Selected() const { return p_selected; }

    friend std::ostream& operator<<(std::ostream &,Splitting_Function_Group &);
  };
}

#endif

#ifndef Splitting_Group_H
#define Splitting_Group_H

#include "Splitting_Function.H"
#include "FastList.H"
#include "Random.H"


namespace APACIC {
  typedef AORGTOOLS::FastList<Splitting_Function> SplFunList; 
  typedef AORGTOOLS::FastList<Splitting_Function>::Iterator SplFunIter;

  /*!
    A group of splitting functions - if for instance more
    than one splitting channel is accessible for a flavour.
  */
  class Splitting_Group : public Splitting_Function {
  protected:
    SplFunList group;
    Splitting_Function * selected;
    // partial sums
    double * partsums;
  public:
    /*! 
      Standard Constructor, if a Splitting_Function is specified it will
      be appended to the list of this particular group.
    */
    Splitting_Group(Splitting_Function * spl=0):partsums(0) {
      if (spl) group.Append(spl);
      selected = spl;
    };
    //! Standard Destructor
    ~Splitting_Group() {
      if (partsums) delete [] partsums;
    };
    //! Summing over all splittings, filling individual values and total sum
    virtual double CrudeInt(double _zmin, double _zmax) {
      if (!partsums) partsums = new double[group.GetLength()];
      lastint = 0;
      int i   = 0;
      for (SplFunIter iter(group);iter();++iter,++i)
	partsums[i] = lastint += iter()->CrudeInt(_zmin,_zmax);
      return lastint;
    };        
    /*!
      A load of varios functions to return values of the particular
      selected function out of the group, see Splitting_Function.H
      for more details.
    */
    double operator()(double z)                        { return (*selected)(z);};
    double GetZ()                                      { return selected->GetZ();};             
    double GetCoupling()                               { return selected->GetCoupling();};
    double GetCoupling(double t)                       { return selected->GetCoupling(t);};
    double GetWeight(double z,double pt2,bool masses)  { return selected->GetWeight(z,pt2,masses);};
    APHYTOOLS::Flavour & GetFlA()                      { return selected->GetFlA();};
    APHYTOOLS::Flavour & GetFlB()                      { return selected->GetFlB();};
    APHYTOOLS::Flavour & GetFlC()                      { return selected->GetFlC();};  

    //! Adds a splitting function to some exisiting group.
    virtual void Add(Splitting_Function * spl) {
      group.Append(spl);
      selected = spl;
    };
    //! Selects one splitting function out of the group.
    virtual void SelectOne() {
      double rr = lastint*AMATOOLS::Ran.get();
      int i;
      for (i=0;partsums[i]<rr;++i) {};
      selected=group[i];
    };


    void PrintStat(int mode=0) {
      if (mode>0) for(int i=0;i<mode;++i) std::cout<<' ';
      std::cout<<"Splitting Group:"<<GetFlA()<<" -> "<<GetFlB()<<" + "<<GetFlC()<<std::endl<<std::endl;
      for (SplFunIter iter(group);iter();++iter) {
	iter()->PrintStat(mode+4);
      };
    
    };
  };
}
#endif




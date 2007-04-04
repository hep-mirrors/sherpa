#ifndef Selector_Bias_h
#define Selector_Bias_h

#include "Selector.H"



namespace ATOOLS {
  class ET_Bias : public Selector_Base {
    std::vector<std::pair<double,double> > m_bounds;
    std::vector<int>   m_sels;
    std::vector<Vec4D> m_moms;
  public:
    ET_Bias(int,int,Flavour *);
    ~ET_Bias();
    void SetRange(std::vector<Flavour>,
		  std::vector<std::pair<double,double> > &);
    bool Trigger(const Vec4D * );
  };


  class PT_Bias : public Selector_Base {
    std::vector<std::pair<double,double> > m_bounds;
    std::vector<int>   m_sels;
    std::vector<Vec4D> m_moms;
  public:
    PT_Bias(int,int,Flavour *);
    ~PT_Bias();
    void SetRange(std::vector<Flavour>,
		  std::vector<std::pair<double,double> > &);
    bool Trigger(const Vec4D * );
  };


  class Eta_Bias : public Selector_Base {
    std::vector<std::pair<double,double> > m_bounds;
    std::vector<int>   m_sels;
    std::vector<Vec4D> m_moms;
  public:
    Eta_Bias(int,int,Flavour *);
    ~Eta_Bias();
    void SetRange(std::vector<Flavour>,
		  std::vector<std::pair<double,double> > &);
    bool Trigger(const Vec4D * );
  };

  struct ordering {
    enum code {
      E  = 0,
      ET = 1,
      PT = 2
    };
  };

  class Leading_Delta_Eta_Bias : public Selector_Base {
    ordering::code m_mode;
    std::vector<std::pair<double,double> > m_bounds;
    std::vector<std::pair<int,int> >       m_sels;
    std::vector<std::pair<Vec4D,Vec4D> >   m_moms;
  public:
    Leading_Delta_Eta_Bias(int,int,Flavour *,ordering::code=ordering::ET);
    ~Leading_Delta_Eta_Bias();
    void SetRange(std::vector<Flavour>,
		  std::vector<std::pair<double,double> > &);
    bool Trigger(const Vec4D * );
  };

  /*
  class Leading_Mass_Bias : public Selector_Base {
    ordering::code m_mode;
    std::vector<std::pair<double,double> > m_bounds;
    std::vector<int>   m_sels;
    std::vector<Vec4D> m_moms;
  public:
    Leading_Mass_Bias(int,int,Flavour *,ordering::code=ordering::ET);
    ~Leading_Mass_Bias();
    void SetRange(std::vector<Flavour>,
		  std::vector<std::pair<double,double> > &);
    bool Trigger(const Vec4D * );
  };

  class Leading_Delta_Phi_Bias : public Selector_Base {
    ordering::code m_mode;
    std::vector<std::pair<double,double> > m_bounds;
    std::vector<int>   m_sels;
    std::vector<Vec4D> m_moms;
  public:
    Leading_Delta_Phi_Bias(int,int,Flavour *,ordering::code=ordering::ET);
    ~Leading_Delta_Phi_Bias();
    void SetRange(std::vector<Flavour>,
		  std::vector<std::pair<double,double> > &);
    bool Trigger(const Vec4D * );
  };

  class Leading_Delta_R_Bias : public Selector_Base {
    ordering::code m_mode;
    std::vector<std::pair<double,double> > m_bounds;
    std::vector<int>   m_sels;
    std::vector<Vec4D> m_moms;
  public:
    Leading_Delta_R_Bias(int,int,Flavour *,ordering::code=ordering::ET);
    ~Leading_Delta_R_Bias();
    void SetRange(std::vector<Flavour>,
		  std::vector<std::pair<double,double> > &);
    bool Trigger(const Vec4D * );
  };

  class Mass_Bias : public Selector_Base {
    std::vector<std::pair<double,double> > m_bounds;
    std::vector<int>   m_sels;
    std::vector<Vec4D> m_moms;
  public:
    Mass_Bias(int,int,Flavour *);
    ~Mass_Bias();
    void SetRange(std::vector<Flavour>,
		  std::vector<std::pair<double,double> > &);
    bool Trigger(const Vec4D * );
  };

  class Delta_Eta_Bias : public Selector_Base {
    std::vector<std::pair<double,double> > m_bounds;
    std::vector<int>   m_sels;
    std::vector<Vec4D> m_moms;
  public:
    Delta_Eta_Bias(int,int,Flavour *);
    ~Delta_Eta_Bias();
    void SetRange(std::vector<Flavour>,
		  std::vector<std::pair<double,double> > &);
    bool Trigger(const Vec4D * );
  };

  class Delta_Phi_Bias : public Selector_Base {
    std::vector<std::pair<double,double> > m_bounds;
    std::vector<int>   m_sels;
    std::vector<Vec4D> m_moms;
  public:
    Eta_Bias(int,int,Flavour *);
    ~Eta_Bias();
    void SetRange(std::vector<Flavour>,
		  std::vector<std::pair<double,double> > &);
    bool Trigger(const Vec4D * );
  };

  class Delta_R_Bias : public Selector_Base {
    std::vector<std::pair<double,double> > m_bounds;
    std::vector<int>   m_sels;
    std::vector<Vec4D> m_moms;
  public:
    Delta_R_Bias(int,int,Flavour *);
    ~Delta_R_Bias();
    void SetRange(std::vector<Flavour>,
		  std::vector<std::pair<double,double> > &);
    bool Trigger(const Vec4D * );
  };
  */
}




/*!
  \class No_Bias
  \brief Dummy Selector, does not select anything.

  This is a dummy selector that does not apply any cuts on the
  phase space. Consequently, the methods are all empty.
*/

/*!
  \class Energy_Bias
  \brief Selector used for energy cuts on final state particles

  Selects according to energy ranges for the particles.
  These ranges are set with the method SetRange for each flavour, the 
  default choice is from 0 to the c.m. energy. However, all energy selection 
  takes place within one selector, this results in the following trigger 
  algorithm: The trigger checks via a loop for all particles whether they 
  are inside their energy range. If not, it returns a 0.
*/

/*!
  \class Angle_Bias
  \brief Selector used to apply cuts on angle between two final state particles

  Selects according to angular ranges (cosines) for pairs of particles.
  These ranges are set with the method SetRange for each flavour pair
  or for each flavour w.r.t. a beam specified by its number, the 
  default choice is from -1 to 1. Again, all angular selection 
  takes place within one selector, this results in the following trigger 
  algorithm: The trigger checks via a loop for all pairs of particles whether 
  they are inside their angular range. If not, it returns a 0.
*/


/*!
  \class Mass_Bias
  \brief Selector used to apply cuts on invariant mass of two final state particles

  Selects according to invariant-mass ranges for pairs of particles.
  These ranges are set with the method SetRange for each flavour pair, the default 
  choice is from 0 to the square of the c.m. energy. Again, all mass selection 
  takes place within one selector, this results in the following trigger 
  algorithm: The trigger checks via a loop for all pairs of particles whether 
  they are inside their mass range. If not, it returns a 0.
*/
#endif

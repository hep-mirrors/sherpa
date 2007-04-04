#ifndef Selector_Bias_h
#define Selector_Bias_h

#include "Selector.H"

namespace ATOOLS {

  struct ordering {
    enum code {
      E  = 0,
      ET = 1,
      PT = 2
    };
  };

  class ET_Bias : public Selector_Base {
    ordering::code m_mode;
    std::vector<std::pair<double,double> > m_bounds;
    std::vector<int>   m_sels;
    std::vector<Vec4D> m_moms;
  public:
    ET_Bias(int,int,Flavour *,ordering::code);
    ~ET_Bias();
    void BuildCuts(Cut_Data *);
    void SetRange(std::vector<Flavour>,
		  std::vector<std::pair<double,double> > &);
    bool Trigger(const Vec4D * );
  };

  class PT_Bias : public Selector_Base {
    ordering::code m_mode;
    std::vector<std::pair<double,double> > m_bounds;
    std::vector<int>   m_sels;
    std::vector<Vec4D> m_moms;
  public:
    PT_Bias(int,int,Flavour *,ordering::code);
    ~PT_Bias();
    void BuildCuts(Cut_Data *);
    void SetRange(std::vector<Flavour>,
		  std::vector<std::pair<double,double> > &);
    bool Trigger(const Vec4D * );
  };

  class Eta_Bias : public Selector_Base {
    ordering::code m_mode;
    std::vector<std::pair<double,double> > m_bounds;
    std::vector<int>   m_sels;
    std::vector<Vec4D> m_moms;
  public:
    Eta_Bias(int,int,Flavour *,ordering::code);
    ~Eta_Bias();
    void BuildCuts(Cut_Data *);
    void SetRange(std::vector<Flavour>,
		  std::vector<std::pair<double,double> > &);
    bool Trigger(const Vec4D * );
  };

  class Mass_Bias : public Selector_Base {
    ordering::code m_mode;
    std::vector<std::pair<double,double> > m_bounds;
    std::vector<int>   m_sels[2];
    std::vector<Vec4D> m_moms[2];
    bool m_idf;
  public:
    Mass_Bias(int,int,Flavour *,ordering::code=ordering::ET);
    ~Mass_Bias();
    void BuildCuts(Cut_Data *);
    void SetRange(std::vector<Flavour>,
		  std::vector<std::pair<double,double> > &);
    bool Trigger(const Vec4D * );
  };

  class Delta_Eta_Bias : public Selector_Base {
    ordering::code m_mode;
    std::vector<std::pair<double,double> > m_bounds;
    std::vector<int>   m_sels[2];
    std::vector<Vec4D> m_moms[2];
    bool m_idf;
  public:
    Delta_Eta_Bias(int,int,Flavour *,ordering::code=ordering::ET);
    ~Delta_Eta_Bias();
    void BuildCuts(Cut_Data *);
    void SetRange(std::vector<Flavour>,
		  std::vector<std::pair<double,double> > &);
    bool Trigger(const Vec4D * );
  };

  class Delta_Phi_Bias : public Selector_Base {
    ordering::code m_mode;
    std::vector<std::pair<double,double> > m_bounds;
    std::vector<int>   m_sels[2];
    std::vector<Vec4D> m_moms[2];
    bool m_idf;
  public:
    Delta_Phi_Bias(int,int,Flavour *,ordering::code=ordering::ET);
    ~Delta_Phi_Bias();
    void BuildCuts(Cut_Data *);
    void SetRange(std::vector<Flavour>,
		  std::vector<std::pair<double,double> > &);
    bool Trigger(const Vec4D * );
  };

  class Delta_R_Bias : public Selector_Base {
    ordering::code m_mode;
    std::vector<std::pair<double,double> > m_bounds;
    std::vector<int>   m_sels[2];
    std::vector<Vec4D> m_moms[2];
    bool m_idf;
  public:
    Delta_R_Bias(int,int,Flavour *,ordering::code=ordering::ET);
    ~Delta_R_Bias();
    void BuildCuts(Cut_Data *);
    void SetRange(std::vector<Flavour>,
		  std::vector<std::pair<double,double> > &);
    bool Trigger(const Vec4D * );
  };

}// end of namespace ATOOLS

#endif

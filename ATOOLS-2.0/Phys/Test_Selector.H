#ifndef Selector_h
#define Selector_h

#include "Vector.H"
#include "Flavour.H"
#include <iostream>
#include <vector>
#include "Selector.H"

namespace ATOOLS {

  // old version to be replaced soon !

  class Test_Selector : public Selector_Base {
    //    Flavour * flavs;
    //    int Nin,Nout;
    const Vec4D* cur_vecs;
    double  max_cos_l_beam, max_cos_l_l, max_cos_l_p, max_cos_l_g;
    double  max_cos_pg_l ,max_cos_g_beam, max_cos_p_beam;
    double  min_m2_q_q;
    double  min_E_l, min_E_q, min_E_p, min_E_g;

    const double m2_ij(const int i,const int j)  const 
    { return (cur_vecs[i] + cur_vecs[j]).Abs2();}
    const double cos_ij(const int i, const int j)  const 
    { Vec3D a= cur_vecs[i];
      Vec3D b= cur_vecs[j];
     return (a*b/(a.Abs()*b.Abs())); }

    std::vector<Selector_Log*> sel_logs;
  public:
    Test_Selector(int _Nin,int _Nout, Flavour * _flavs);
    bool Trigger(const Vec4D*);
    virtual void Output() { 
      int count = 0;
      for (std::vector<Selector_Log*>::iterator it = sel_logs.begin(); it!=sel_logs.end(); ++it) {
	(*it)->Output();
	count += (*it)->Rejections();
      }
    }
    ~Test_Selector() { 
      delete [] m_fl; 
      m_fl = 0;
      for (std::vector<Selector_Log*>::iterator it = sel_logs.begin(); it!=sel_logs.end(); ++it) {
	delete (*it);
      }
    }
    void BuildCuts(Cut_Data *);
  };

}
#endif



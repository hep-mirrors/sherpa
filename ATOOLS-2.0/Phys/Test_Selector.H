#ifndef Selector_h
#define Selector_h

#include "Vector.H"
#include "Flavour.H"
#include <iostream>
#include <vector>
#include "Selector.H"

namespace APHYTOOLS {

  // old version to be replaced !

  class Test_Selector : public Selector_Base {
    Flavour * flavs;
    int Nin,Nout;
    const AMATOOLS::vec4d* cur_vecs;
    double  max_cos_l_beam, max_cos_l_l, max_cos_l_p, max_cos_l_g;
    double  max_cos_pg_l ,max_cos_g_beam, max_cos_p_beam;
    double  min_m2_q_q;
    double  min_E_l, min_E_q, min_E_p, min_E_g;

    const double m2_ij(const int i,const int j)  const 
    { return (cur_vecs[i] + cur_vecs[j]).abs2();}
    const double cos_ij(const int i, const int j)  const 
    { AMATOOLS::vec3d a= cur_vecs[i];
      AMATOOLS::vec3d b= cur_vecs[j];
     return (a*b/(a.abs()*b.abs())); }

    std::vector<Selector_Log*> sel_logs;
  public:
    Test_Selector(int _Nin,int _Nout, Flavour * _flavs);
    bool Trigger(const AMATOOLS::vec4d*);
    virtual void Output() { 
      int count = 0;
      for (std::vector<Selector_Log*>::iterator it = sel_logs.begin(); it!=sel_logs.end(); ++it) {
	(*it)->Output();
	count += (*it)->Rejections();
      }
      std::cout<<"Total Number of Rejections: "<<count<<std::endl;
    }
    ~Test_Selector() { 
      std::cout<<" deleting Test_Selector "<<std::endl;
      delete [] flavs; 
      for (std::vector<Selector_Log*>::iterator it = sel_logs.begin(); it!=sel_logs.end(); ++it) {
	delete (*it);
      }
      std::cout<<" ... done "<<std::endl;
    }
    void BuildCuts(Cut_Data *);
  };

}
#endif



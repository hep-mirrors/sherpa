#ifndef Vegas_H
#define Vegas_H

#include "Message.H"
#include <list>

namespace PHASIC {

  class Vegas {
    std::string m_name;
    int m_nd;
    unsigned long m_nevt,m_nopt,m_snevt;
    double m_alpha;
    double m_weight,m_sum,m_sum2;
    double m_nc;
    double **p_xi;
    double *p_x,*p_xin,*p_r;
    double **p_d,**p_di,*p_dt;
    double *p_chi,*p_bestchi,**p_bestxi;
    int *p_ia,*p_opt;
    int **p_hit;
    int m_dim,m_mode,m_on,m_autooptimize;
    static int s_onext;
    void Rebin(double rc, double * xi);
  public:
    Vegas(int dim,int ndx,const std::string &,int=1);
    ~Vegas();
    double* GeneratePoint(const double * _ran);
    double GenerateWeight(double *);
    void AddPoint(double,double *);

    void Reset();
    void Optimize();
    void EndOptimize();
    void ConstChannel(int i) { p_opt[i]=0; }

    void WriteOut(const std::string &);
    void ReadIn(const std::string &);
    void WriteHistos(const std::string &);

    void SetAutoOptimize(int ao)    { m_autooptimize=ao*m_nd; }

    static inline void SetOnExternal(const bool on) { s_onext=on;     }
    static inline bool OnExternal()                 { return s_onext; }
  };

}

#endif

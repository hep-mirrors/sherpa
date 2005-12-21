#ifndef Statistics_Observable_H
#define Statistics_Observable_H

#include "Primitive_Observable_Base.H"
#include <iostream>

namespace ANALYSIS {
  struct Statistics_Data {
    unsigned long nevt, xsnevt;
    double  weightsum, xssum, sudsum;
    unsigned long nblobssum;
    unsigned long nplsum;

    unsigned long nevt_cut;
    double  weightsum_cut;    
    unsigned long nblobssum_cut;
    unsigned long nplsum_cut;
    Statistics_Data() :
     nevt(0), xsnevt(0), weightsum(0.), xssum(0.0), sudsum(0.0), nblobssum(0), nplsum(0),
     nevt_cut(0), weightsum_cut(0.), nblobssum_cut(0), nplsum_cut(0) {} 

    void EndEvaluation(double scale=1.);
    void Output(const std::string & pname);
  };

  typedef std::map<std::string,Statistics_Data>  Statistics_Map;


  class Statistics_Observable : public Primitive_Observable_Base {  
  protected:
    unsigned long m_nevt;
    Statistics_Map  m_signal_process_statistics;
  public:
    Statistics_Observable(const std::string & inlistname, int mode=0);
    void Evaluate(const ATOOLS::Blob_List & blobs,double weight=1.,int ncount=1);

    void EndEvaluation(double scale=1.);
    void Output(const std::string & pname);

    Primitive_Observable_Base * Copy() const;

  };

  std::ostream & operator<<(std::ostream & str, const Statistics_Data & sd);
}

#endif

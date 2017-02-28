#ifndef ATOOLS__Org__My_MPI_H
#define ATOOLS__Org__My_MPI_H

#include "ATOOLS/Org/CXXFLAGS.H"
#ifdef USING__MPI
#include "mpi.h"
#endif
#ifdef USING__Threading
inline int pthread_cond_signal
(pthread_cond_t *c,pthread_mutex_t *m)
{
  pthread_mutex_lock(m);
  int r(pthread_cond_signal(c));
  pthread_mutex_unlock(m);
  return r;
}
#endif

#include <vector>

namespace ATOOLS {

  class Data_Reader;

  class My_MPI {
  private:

#ifdef USING__MPI
    MPI::Intracomm *p_comm;
#endif

    void SetMPIRecv(std::vector<int> r);

  public:

    My_MPI();

    ~My_MPI();

    void SetUpSendRecv(Data_Reader *const read);

#ifdef USING__MPI
    inline MPI::Intracomm *MPIComm() { return p_comm; }
#endif

  };// end of class My_MPI

  extern My_MPI *mpi;

  void Abort(const int mode=0);

}// end of namespace ATOOLS

#endif

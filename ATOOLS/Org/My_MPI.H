#ifndef ATOOLS__Org__My_MPI_H
#define ATOOLS__Org__My_MPI_H

#include "ATOOLS/Org/CXXFLAGS.H"

#ifdef USING__MPI
#include "mpi.h"
#endif

#ifdef USING__Threading
#include <pthread.h>
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

  class My_MPI {

  public:

    My_MPI();

    void PrintRankInfo();

#ifdef USING__MPI

    void Barrier() {
      MPI_Barrier(m_comm);
    }

    int Rank() {
      int rank;
      MPI_Comm_rank(m_comm, &rank);
      return rank;
    }

    int Size() {
      int size;
      MPI_Comm_size(m_comm, &size);
      return size;
    }

    void Bcast(void* buffer, int count, MPI_Datatype type, int root=0) {
      MPI_Bcast(buffer, count, type, root, m_comm);
    }

    void Allreduce(void* buffer, int count, MPI_Datatype type, MPI_Op op) {
      MPI_Allreduce(MPI_IN_PLACE, buffer, count, type, op, m_comm);
    }

    void Allgather(const void *sendbuf, int sendcount, MPI_Datatype sendtype,
                         void *recvbuf, int recvcount, MPI_Datatype recvtype) {
      MPI_Allgather(sendbuf, sendcount, sendtype, recvbuf, recvcount, recvtype,
                    m_comm);
    }

#endif

  private:

#ifdef USING__MPI
    MPI_Comm m_comm;
#endif

  };// end of class My_MPI

  extern My_MPI* mpi;

  void Abort(const int mode=0);

}// end of namespace ATOOLS

#endif

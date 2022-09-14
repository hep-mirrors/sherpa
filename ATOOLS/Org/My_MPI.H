#ifndef ATOOLS__Org__My_MPI_H
#define ATOOLS__Org__My_MPI_H

#include "ATOOLS/Org/CXXFLAGS.H"

#ifdef USING__MPI
#include "mpi.h"
#include <string>
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
  private:

#ifdef USING__MPI
    MPI_Comm m_comm;
#endif
    int m_rank, m_size, m_myrank, m_mysize;

  public:

    My_MPI();

    void PrintRankInfo();
    void PrintRank();

    inline void SetMyRank(const int rank) { m_myrank=rank; }
    inline void SetMySize(const int size) { m_mysize=size; }

#ifdef USING__MPI

    void Barrier() {
      MPI_Barrier(m_comm);
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

    void Recv(void *buf, int count, MPI_Datatype datatype, int source, int tag)
    {
      MPI_Recv(buf, count, datatype, source, tag, m_comm, MPI_STATUS_IGNORE);
    }

    int Send(const void *buf, int count, MPI_Datatype datatype, int dest,
             int tag)
    {
      return MPI_Send(buf, count, datatype, dest, tag, m_comm);
    }

    void Gather(void *sendbuf, int sendcount, MPI_Datatype sendtype,
		void *recvbuf, int recvcount, MPI_Datatype recvtype,
		int root)
    {
      MPI_Gather(sendbuf, sendcount, sendtype,
		 recvbuf, recvcount, recvtype, root, m_comm);
    }

    // Convenience method to get the maximum of some value from all ranks.
    int Allmax(int);

    // Convenience method to gather strings from all ranks into an array.
    // Allgather is used, such that the returned vector is filled for all
    // ranks. Note that this involves a non-trivial amount of communication
    // and should not be used too often.
    std::vector<std::string> AllgatherStrings(const std::string&);

#endif

    inline int Rank() { return m_rank; }
    inline int Size() { return m_size; }

    inline int MyRank() { return m_myrank; }
    inline int MySize() { return m_mysize; }

  };// end of class My_MPI

  extern My_MPI* mpi;

  void Abort(const int mode=0);

  class MPI_Object {
  protected:

    static std::vector<MPI_Object*> s_objects;

  public:

    virtual void MPISync() = 0;

    int Communicate(const int mode=0);

  };// end of class MPI_Object

}// end of namespace ATOOLS

#endif

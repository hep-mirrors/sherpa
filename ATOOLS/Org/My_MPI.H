#ifndef ATOOLS__Org__My_MPI_H
#define ATOOLS__Org__My_MPI_H

#include "ATOOLS/Org/CXXFLAGS.H"

#ifdef USING__MPI
#include "mpi.h"
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

#endif

    inline int Rank() { return m_rank; }
    inline int Size() { return m_size; }

    inline int MyRank() { return m_myrank; }
    inline int MySize() { return m_mysize; }

  };// end of class My_MPI

  extern My_MPI* mpi;

  class MPI_Object {
  protected:

    static std::vector<MPI_Object*> s_objects;

  public:

    virtual void MPISync() = 0;

    int Communicate(const int mode=0);

  };// end of class MPI_Object

}// end of namespace ATOOLS

#endif

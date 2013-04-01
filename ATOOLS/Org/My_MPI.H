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

    bool m_hassend, m_hasrecv;

#ifdef USING__MPI
    MPI::Intracomm m_send, m_recv;
#endif

  public:

    My_MPI();

    ~My_MPI();

    void SetMPIRecv(std::vector<int> r);

    bool HasMPISend() const;
    bool HasMPIRecv() const;

#ifdef USING__MPI
    inline MPI::Intracomm &MPISend() { return m_send; }
    inline MPI::Intracomm &MPIRecv() { return m_recv; }
#endif

  };// end of class My_MPI

  extern My_MPI *mpi;

}// end of namespace ATOOLS

#endif

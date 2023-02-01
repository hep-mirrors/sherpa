#include "ATOOLS/Org/My_MPI.H"

#include "ATOOLS/Org/Message.H"

#include <csignal>
#include <unistd.h>

using namespace ATOOLS;

My_MPI *ATOOLS::mpi(NULL);

My_MPI::My_MPI(): m_rank(0), m_size(1)
{
#ifdef USING__MPI
  m_comm = MPI_COMM_WORLD;
  MPI_Comm_rank(m_comm,&m_rank);
  MPI_Comm_size(m_comm,&m_size);
#endif
  m_myrank=m_rank;
  m_mysize=m_size;
}

void My_MPI::PrintRankInfo()
{
#ifdef USING__MPI
  const int size = Size();
  if (size > 1)
    msg_Info() << METHOD << "(): Running on " << size << " ranks." << std::endl;
#endif
}

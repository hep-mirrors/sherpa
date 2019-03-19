#include "ATOOLS/Org/My_MPI.H"

#include "ATOOLS/Org/Message.H"

#include <csignal>
#include <unistd.h>

using namespace ATOOLS;

My_MPI *ATOOLS::mpi(NULL);

My_MPI::My_MPI()
{
#ifdef USING__MPI
  m_comm = MPI_COMM_WORLD;
#endif
}

void My_MPI::PrintRankInfo()
{
#ifdef USING__MPI
  const int size = Size();
  if (size > 1)
    msg_Info() << METHOD << "(): Running on " << size << " ranks." << std::endl;
#endif
}

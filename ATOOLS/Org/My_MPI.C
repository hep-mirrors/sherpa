#include "ATOOLS/Org/My_MPI.H"

#include "ATOOLS/Org/Message.H"

#include <csignal>
#include <unistd.h>

using namespace ATOOLS;

My_MPI* ATOOLS::mpi {nullptr};

My_MPI::My_MPI()
{
#ifdef USING__MPI
  m_comm = MPI_COMM_WORLD;
#endif
}

void My_MPI::PrintRankInfo()
{
#ifdef USING__MPI
  const auto size = Size();
  if (size > 1)
    msg_Info() << METHOD << "(): Running on " << size << " ranks." << std::endl;
#endif
}

void ATOOLS::Abort(const int mode)
{
#ifdef USING__MPI
  MPI_Abort(MPI_COMM_WORLD, 1 + mode);
#else
  if (mode)
    kill(getpid(), 9);
  abort();
#endif
}

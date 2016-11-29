#include "ATOOLS/Org/My_MPI.H"

#include "ATOOLS/Org/Run_Parameter.H"
#include "ATOOLS/Org/Data_Reader.H"
#include "ATOOLS/Org/Shell_Tools.H"
#include "ATOOLS/Org/Message.H"

#include <csignal>
#include <cstddef>
#include <cstring>
#include <unistd.h>

using namespace ATOOLS;

My_MPI *ATOOLS::mpi(NULL);

My_MPI::My_MPI()
{
#ifdef USING__MPI
  p_comm=&MPI::COMM_WORLD;
#endif
}

My_MPI::~My_MPI()
{
}

void My_MPI::SetUpSendRecv(Data_Reader *const read)
{
#ifdef USING__MPI
  int size=MPI::COMM_WORLD.Get_size();
  if (size>1) {
    msg_Info()<<METHOD<<"(): Running on "<<size<<" ranks."<<std::endl;
  }
#endif
}

void ATOOLS::Abort(const int mode)
{
#ifdef USING__MPI
  MPI::COMM_WORLD.Abort(1+mode);
#else
  if (mode) kill(getpid(),9);
  abort();
#endif  
}

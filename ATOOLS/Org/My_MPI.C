#include "ATOOLS/Org/My_MPI.H"

#include "ATOOLS/Org/Run_Parameter.H"
#include "ATOOLS/Org/Data_Reader.H"
#include "ATOOLS/Org/Shell_Tools.H"
#include "ATOOLS/Org/Message.H"

#include <stddef.h>
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
    double starttime=rpa->gen.Timer().RealTime();
    msg_Info()<<METHOD<<"(): Analyzing MPI environment {\n";
    int rank=MPI::COMM_WORLD.Get_rank(), pid(getpid()), hlen;
    char host[MPI_MAX_PROCESSOR_NAME];
    MPI_Get_processor_name(host,&hlen);
    if (rank==0) {
      msg_Info()<<"  Rank "<<rank<<", pid "<<pid
		<<" running on "<<host<<".\n";
      char mhost[MPI_MAX_PROCESSOR_NAME];
      MPI_Get_processor_name(mhost,&hlen);
      for (int tag=1;tag<size;++tag) {
	MPI::COMM_WORLD.Recv(&pid,1,MPI::INT,MPI::ANY_SOURCE,tag);
	MPI::COMM_WORLD.Recv(host,MPI_MAX_PROCESSOR_NAME,
			     MPI::CHAR,MPI::ANY_SOURCE,tag);
	msg_Info()<<"  Rank "<<tag<<", pid "<<pid
		  <<" running on "<<host<<"."<<std::endl;
      }
      double diff=rpa->gen.Timer().RealTime()-starttime;
      msg_Info()<<"} -> "<<FormatTime(size_t(diff))<<" elapsed"<<std::endl;
    }
    else {
      MPI::COMM_WORLD.Send(&pid,1,MPI::INT,0,rank);
      MPI::COMM_WORLD.Send(host,MPI_MAX_PROCESSOR_NAME,MPI::CHAR,0,rank);
    }
  }
#endif
}

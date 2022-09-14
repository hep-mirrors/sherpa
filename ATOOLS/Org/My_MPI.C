#include "ATOOLS/Org/My_MPI.H"

#include "ATOOLS/Org/Message.H"
#include "ATOOLS/Org/Run_Parameter.H"
#include "ATOOLS/Org/Exception.H"

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

std::vector<MPI_Object*> MPI_Object::s_objects;

int MPI_Object::Communicate(const int mode)
{
#ifndef USING__MPI
  MPISync();
#else
#ifdef DEBUG__MPI_Sync
  PRINT_FUNC("");
#endif
  size_t nobj(s_objects.size());
  std::vector<int> ids(2*nobj,0);
  size_t idx(0);
  for (;idx<nobj;++idx)
    if (s_objects[idx]==this) break;
  if (idx==nobj) THROW(fatal_error,"Internal error");
  ids[idx]=1;
  ids[nobj+idx]=mode;
#ifdef DEBUG__MPI_Sync
  std::cout<<"Rank "<<mpi->Rank()<<" targets "<<typeid(*this).name()
	   <<", index = "<<idx<<" at "<<rpa->gen.NumberOfGeneratedEvents()<<" events "<<std::endl;
#endif
  mpi->Allreduce(&ids[0],ids.size(),MPI_INT,MPI_MAX);
#ifdef DEBUG__MPI_Sync
  std::cout<<"Rank "<<mpi->Rank()<<" has "<<ids<<std::endl;
#endif
  for (size_t i(0);i<nobj;++i)
    if (ids[i]) {
#ifdef DEBUG__MPI_Sync
      std::cout<<"Rank "<<mpi->Rank()<<" processes "<<typeid(*s_objects[i]).name()
	       <<", index = "<<i<<" at "<<rpa->gen.NumberOfGeneratedEvents()<<" events "<<std::endl;
#endif
      s_objects[i]->MPISync();
    }
  int nbreak(0);
  for (size_t i(nobj);i<2*nobj;++i) nbreak+=ids[i];
  if (nbreak) {
    for (size_t i(0);i<nobj;++i)
      if (ids[i] && !ids[nobj+i]) {
#ifdef DEBUG__MPI_Sync
	std::cout<<"Rank "<<mpi->Rank()<<" continues"<<std::endl;
#endif
	return -1;
      }
  }
#endif
  return 0;
}

#include "ATOOLS/Org/My_MPI.H"

#include <stddef.h>

using namespace ATOOLS;

My_MPI *ATOOLS::mpi(NULL);

My_MPI::My_MPI():
  m_hassend(false), m_hasrecv(false)
{
}

My_MPI::~My_MPI()
{
#ifdef USING__MPI
  if (m_hassend) m_send.Free();
  if (m_hasrecv) m_recv.Free();
#endif  
}

void My_MPI::SetMPIRecv(std::vector<int> r)
{
#ifdef USING__MPI
  int rank=MPI::COMM_WORLD.Get_rank();
  if (rank==0) {
    m_hasrecv=true;
    m_recv=MPI::COMM_WORLD.Split(rank,rank);
    m_send=MPI::COMM_WORLD.Split(MPI_UNDEFINED,rank);
  }
  else {
    if (r[0]==0) {
      m_hassend=m_hasrecv=true;
      m_send=MPI::COMM_WORLD.Split(r[0],rank);
      m_recv=MPI::COMM_WORLD.Split(rank,rank);
    }
    else {
      m_hassend=true;
      m_recv=MPI::COMM_WORLD.Split(MPI_UNDEFINED,rank);
      m_send=MPI::COMM_WORLD.Split(r[0],rank);
    }
  }
#endif
}

bool My_MPI::HasMPISend() const
{
  return m_hassend;
}

bool My_MPI::HasMPIRecv() const
{
#ifdef USING__MPI
  if (m_hasrecv) return m_recv.Get_size()>1;
#endif
  return false;
}


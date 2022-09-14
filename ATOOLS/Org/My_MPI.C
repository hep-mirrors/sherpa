#include "ATOOLS/Org/My_MPI.H"

#include "ATOOLS/Org/Message.H"
#include "ATOOLS/Org/Run_Parameter.H"
#include "ATOOLS/Org/Exception.H"

#include <csignal>
#include <unistd.h>
#include <algorithm>

using namespace ATOOLS;

My_MPI* ATOOLS::mpi {nullptr};

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
  const auto size = Size();
  if (size > 1)
    msg_Info() << METHOD << "(): Running on " << size << " ranks." << std::endl;
#endif
}

void My_MPI::PrintRank() {
#ifdef USING__MPI
  // We can not use msg_Out and its friends, because they ignore calls from
  // all but the 0th rank.
  std::cout << "MPI Rank: " << Rank() << "\n";
#endif
}

#ifdef USING__MPI

int My_MPI::Allmax(int i)
{
  const int n_ranks {mpi->Size()};
  std::vector<int> all_i;
  if (mpi->Rank() == 0) {
    all_i.resize(n_ranks, 0);
  }
  int max;
  mpi->Gather(&i, 1, MPI_INT, &(all_i[0]), 1, MPI_INT, 0);
  if (mpi->Rank() == 0) {
    max = *std::max_element(all_i.cbegin(), all_i.cend());
  }
  mpi->Bcast(&max, 1, MPI_INT);
  return max;
}

std::vector<std::string> My_MPI::AllgatherStrings(const std::string& s) {

  const int n_ranks {Size()};
  const int s_size {static_cast<int>(s.size())};

  /*
   * Now, we Allgather the string lengths, so we can create the buffer into
   * which we'll receive all the strings.
   */

  int* recvcounts {(int*)malloc(n_ranks * sizeof(int))};

  Allgather(&s_size, 1, MPI_INT, recvcounts, 1, MPI_INT);

  /*
   * Figure out the length of the resulting combined string, and the
   * displacements for each rank (i.e. where to put their individual strings).
   */

  int totlen = 0;
  int* displs = NULL;
  char* totalstring = NULL;

  displs = (int*)malloc(n_ranks * sizeof(int));

  displs[0] = 0;
  totlen += recvcounts[0] + 1;

  for (int i = 1; i < n_ranks; i++) {
    totlen += recvcounts[i] + 1; /* plus one for '\0' after each word */
    displs[i] = displs[i - 1] + recvcounts[i - 1] + 1;
  }

  /* allocate string, pre-fill with null terminators */
  totalstring = (char*)malloc(totlen * sizeof(char));
  for (int i = 0; i < totlen; i++)
    totalstring[i] = '\0';

  /*
   * Now we have the receive buffer, counts, and displacements, and can gather
   * all the strings into the combined totalstring buffer.
   */

  MPI_Allgatherv(s.c_str(), s_size, MPI_CHAR, totalstring, recvcounts, displs,
                 MPI_CHAR, MPI_COMM_WORLD);


  /* put substrings from totalstring into a vector */
  std::vector<std::string> allstrings(1, "");
  for (int i {0}; i < totlen - 1; i++) {
    if (totalstring[i] == '\0')
      allstrings.push_back("");
    else
      allstrings.back() += totalstring[i];
  }

  free(totalstring);
  free(displs);
  free(recvcounts);

  return allstrings;
}

#endif

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

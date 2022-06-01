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

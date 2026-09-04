#include "ATOOLS/Org/CXXFLAGS.H"

#include <iostream>
#include <string>
#include <vector>
#include <unistd.h>

#include <highfive/H5File.hpp>
#include <highfive/H5DataSet.hpp>

#include "AddOns/HDF5/LHEH5.H"
#include "AddOns/HDF5/LHEH5_Reader_Base.H"

#include "PHASIC++/Main/Event_Reader.H"
#include "ATOOLS/Org/Run_Parameter.H"
#include "ATOOLS/Org/My_MPI.H"
#include "ATOOLS/Org/Message.H"

using namespace HighFive;
using namespace PHASIC;
using namespace ATOOLS;

namespace LHEH5 {

  // Reads parton-level events from LHEH5 files on disk. Each rank works
  // on its own contiguous slice of the "events" dataset, pulling it in
  // in chunks of m_ncache events; once a file is exhausted, the next one
  // from m_files takes over.
  class HDF5_Reader: public LHEH5_Reader_Base {
  private:

    LHEFile *p_file;
#if defined(USING__MPI) && defined(H5_HAVE_PARALLEL)
    MPI_Info m_info;
#endif
    size_t m_ifile, m_ilaststart, m_inextstart, m_ntotal;

    // Collective I/O is only meaningful for the parallel HDF5 driver, so
    // the transfer properties are built here rather than inside LHEFile.
    DataTransferProps TransferProps() const
    {
      auto xfer_props = DataTransferProps{};
#if defined(USING__MPI) && defined(H5_HAVE_PARALLEL)
      xfer_props.add(UseCollectiveIO{});
#endif
      return xfer_props;
    }

  protected:

    const LHEFile &CurrentFile() const override { return *p_file; }
    size_t NEventsInBuffer() const override { return p_file->NEvents(); }

    bool NextBuffer() override
    {
      if ((m_ifile + 1 >= m_files.size() && m_inextstart == 0) ||
          Communicate() < 0)
        return false;
      m_ievt = 0;
      return true;
    }

    std::string ExhaustionMessage() const override
    {
      return "There are no more events in the input file '" + m_files[m_ifile]
             + "' and there are no more input files for this jet "
               "multiplicity.";
    }

  public:

    HDF5_Reader(const Event_Reader_Key &key):
      LHEH5_Reader_Base(key, "HDF5_CACHE_SIZE", "an HDF5 file"),
      m_ifile(0), m_ilaststart(0), m_inextstart(0), m_ntotal(0)
    {
      Settings& s {Settings::GetMainSettings()};

#if defined(USING__MPI) && defined(H5_HAVE_PARALLEL)
      MPI_Info_create(&m_info);
      for (const auto& key : s["HDF5_MPIIO_PARAMS"].GetKeys()) {
        const auto val {
            s["HDF5_MPIIO_PARAMS"][key].SetDefault("").Get<std::string>()};
	msg_Info()<<METHOD<<"(): Add MPIIO parameters '"
		  <<key<<"' -> '"<<val<<"'\n";
	MPI_Info_set(m_info,key.c_str(),val.c_str());
      }
#endif

      p_file = OpenFile(m_files[m_ifile]);
    }

    ~HDF5_Reader()
    {
      delete p_file;
    }

    void PrintStatistics(std::ostream& o) override
    {
      for (int i{0}; i < m_files.size(); ++i) {
        o << "    " << m_files[i] << ": ";
        if (i < m_ifile)
          o << "all events read";
        else if (i == m_ifile)
	  o << m_ilaststart + m_ievt << " of " << m_ntotal / mpi->MySize()
	    << " events read ("
	    << (m_ilaststart + m_ievt) * 1000 / (m_ntotal / mpi->MySize()) /
		   10.0
	    << " %)";
        else
          o << "not yet opened";
        if (mpi->MySize() > 1) o << " on rank 0";
        o << '\n';
      }
    }

    void MPISync() override
    {
      delete p_file;
      if (m_inextstart == 0) ++m_ifile;
      p_file = OpenFile(m_files[m_ifile]);
    }

    LHEFile *OpenFile(const std::string &fname)
    {
      m_ievt=0;
      int size(mpi->MySize()), rank(mpi->MyRank());
#if defined(USING__MPI) && defined(H5_HAVE_PARALLEL)
      FileAccessProps fapl;
      fapl.add(MPIOFileAccess{MPI_COMM_WORLD,m_info});
      fapl.add(MPIOCollectiveMetadata{});
      File file(fname,File::ReadOnly,fapl);
#else
      File file(fname,File::ReadOnly);
#endif
      LHEFile *e(new LHEFile());
      e->ReadHeader(file,TransferProps());
      m_totalxs=e->TotalXS();
      m_unitwgt=e->UnitWeight();
      if (e->Version()[0]==2 &&
	  e->Version()[1]==0 &&
	  e->Version()[2]==0) m_unitwgt*=rpa->Picobarn();
      m_ntotal = file.getDataSet("events").getSpace().getDimensions().front();
      msg_Info() << "Reading events from file: '" << m_files[m_ifile]
                 << "' (total per rank: " << m_ntotal / mpi->MySize();
      if (m_inextstart > 0)
        msg_Info() << ", read: " << m_inextstart;
      msg_Info() << ").\n";
#if defined(USING__MPI) && defined(H5_HAVE_PARALLEL)
      mpi->Bcast(&m_ntotal,1,MPI_LONG_INT);
#endif
      size_t iStart(rank*m_ntotal/size);
      size_t iStop((rank+1)*m_ntotal/size-1);
      if (rank==size-1) iStop=m_ntotal-1;
      size_t nread(std::min(m_ncache,iStop-iStart-m_inextstart+1));
      e->ReadEvents(file,iStart+m_inextstart,nread,TransferProps());
      m_ilaststart = m_inextstart;
      if (iStart+(m_inextstart+=nread)>iStop) m_inextstart=0;
      return e;
    }

  };// end of class HDF5_Reader

}// end of namespace LHEH5

using namespace LHEH5;

DECLARE_GETTER(HDF5_Reader,"HDF5",Event_Reader,Event_Reader_Key);

Event_Reader *ATOOLS::Getter<Event_Reader,Event_Reader_Key,HDF5_Reader>::
operator()(const Event_Reader_Key &args) const
{
  return new HDF5_Reader(args);
}

void ATOOLS::Getter<Event_Reader,Event_Reader_Key,HDF5_Reader>::
PrintInfo(std::ostream &str,const size_t width) const
{
  str<<"HDF5 reader (version 2)";
}

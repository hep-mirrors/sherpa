#include "ATOOLS/Org/CXXFLAGS_PACKAGES.H"
#include "SHERPA/Tools/Output_Base.H"
#ifdef USING__GZIP
#include "ATOOLS/Org/Gzip_Stream.H"
#else
#include <fstream>
#endif
#include "ATOOLS/Phys/Variations.H"
#include "ATOOLS/Phys/NLO_Subevt.H"
#include "ATOOLS/Org/Shell_Tools.H"
#include "ATOOLS/Org/MyStrStream.H"
#include "ATOOLS/Org/Exception.H"
#include "ATOOLS/Org/My_MPI.H"
#include "ATOOLS/Org/Run_Parameter.H"
#include "ATOOLS/Org/Scoped_Settings.H"

using namespace SHERPA;
using namespace ATOOLS;

namespace SHERPA {

  class Output_Weights: public Output_Base {
  private:

    std::string m_basename, m_ext;

#ifdef USING__GZIP
    ATOOLS::ogzstream m_outstream;
#else
    std::ofstream     m_outstream;
#endif

  public:

    Output_Weights(const Output_Arguments &args):
      Output_Base{ "Weights" }
    {
      MyStrStream basename;
      basename << args.m_outpath << "/" << args.m_outfile;
      m_ext=".wts";
#ifdef USING__GZIP
      m_ext += ".gz";
#endif
#ifdef USING__MPI
      if (mpi->Size()>1) {
        basename << "_" << mpi->Rank();
      }
#endif
      m_basename = basename.str();
      m_outstream.open((m_basename+m_ext).c_str());
      if (!m_outstream.good())
	THROW(fatal_error, "Could not open event file "+m_basename+m_ext+".");
      const int precision =
        Settings::GetMainSettings()["OUTPUT_PRECISION"].Get<int>();
      m_outstream.precision(precision);
    }

    void Header()
    {
      m_outstream << "# EventNumber Nominal";
      static std::array<Variations_Type, 2> types = {Variations_Type::qcd, Variations_Type::qcut};
      for (const auto type : types) {
        size_t numqcdvars = s_variations->Size(type);
        for (size_t i(0);i<numqcdvars;++i) {
          m_outstream << ' ' << s_variations->GetVariationNameAt(i, type);
        }
      }
      m_outstream << '\n';
    }

    ~Output_Weights()
    {
      m_outstream.close();
    }

    void Output(Blob_List* blobs)
    {
      Blob_Data_Base *sd((*blobs->FindFirst(btp::Signal_Process))["NLO_subeventlist"]);
      NLO_subevtlist *subs=sd?sd->Get<NLO_subevtlist*>():NULL;
      if (subs==NULL) {
        Output(blobs->WeightsMap());
        m_outstream<<"\n";
      } else {
	for (size_t j(0);j<subs->size();++j) {
          Output((*subs)[j]->m_results);
          m_outstream<<"\n";
	}
      }
    }

    void Output(const Weights_Map& wgtmap)
    {
      static std::array<Variations_Type, 2> types = {
          Variations_Type::qcd, Variations_Type::qcut};
      const auto nom = wgtmap.Nominal();
      m_outstream<<rpa->gen.NumberOfGeneratedEvents();
      if (IsZero(nom))
        m_outstream<<' '<<0.0;
      else
        m_outstream<<' '<<nom;
      for (const auto type : types) {
        size_t numvars = s_variations->Size(type);
        auto wgts = wgtmap.Combine(type);
        const auto relfac = wgtmap.NominalIgnoringVariationType(type);
        for (size_t i(0);i<numvars;++i) {
          if (IsZero(wgts.var(i) * relfac))
            m_outstream << " " << 0.0;
          else
            m_outstream << " " << wgts.var(i) * relfac;
        }
      }
    }

    void ChangeFile()
    {
      m_outstream.close();
      std::string newname(m_basename+m_ext);
      for (size_t i(0);FileExists(newname);
	   newname=m_basename+"."+ToString(++i)+m_ext);
      m_outstream.open(newname.c_str());
      if (!m_outstream.good())
	THROW(fatal_error, "Could not open event file "+newname+".");
    }

  };

}

DECLARE_GETTER(Output_Weights,"Weights",Output_Base,Output_Arguments);
Output_Base *ATOOLS::Getter<Output_Base,Output_Arguments,Output_Weights>::
operator()(const Output_Arguments &args) const
{ return new Output_Weights(args); }

void ATOOLS::Getter<Output_Base,Output_Arguments,Output_Weights>::
PrintInfo(std::ostream &str,const size_t width) const
{ str<<"weight output"; }


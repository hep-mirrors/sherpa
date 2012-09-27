#include "AddOns/Analysis/Observables/Primitive_Observable_Base.H"
#include "ATOOLS/Org/MyStrStream.H"
#include "ATOOLS/Org/Data_Reader.H"
#include "ATOOLS/Org/Exception.H"

using namespace ANALYSIS;
using namespace ATOOLS;

namespace ANALYSIS {

  class Blobdata : public Primitive_Observable_Base {
    std::string m_datakey;

  public:

    Blobdata(const std::string& datakey,
             int type, double xmin, double xmax, int nbins) :
      Primitive_Observable_Base(type, xmin, xmax, nbins), m_datakey(datakey)
    {
      m_name="Blobdata_"+m_datakey+".dat";
    }

    ~Blobdata()
    {
    }


    void Evaluate(const ATOOLS::Blob_List & blobs, double weight, double ncount)
    {
      if (m_datakey=="") return;
      Blob *signal(blobs.FindFirst(btp::Signal_Process));
      if (signal) {
        Blob_Data_Base *facscale((*signal)[m_datakey]);
        if (facscale) {
          p_histo->Insert(sqrt(facscale->Get<double>()), weight, ncount); 
        }
        else PRINT_INFO("Key "<<m_datakey<<" not found in event.");
      }
    }


    Primitive_Observable_Base * Copy() const
    {
      // don't duplicate
      return new Blobdata("", m_type, m_xmin, m_xmax, m_nbins);
    }

  };// end of class Blobdata

}



DECLARE_GETTER(Blobdata_Getter,"Blobdata",
	       Primitive_Observable_Base,Argument_Matrix);

Primitive_Observable_Base *
Blobdata_Getter::operator()(const Argument_Matrix &parameters) const
{
  if (parameters.size()<1) return NULL;
  if (parameters.size()==1) {
    if (parameters[0].size()<5) return NULL;
    return new Blobdata(parameters[0][4],
                        HistogramType(parameters[0][3]),
                        ATOOLS::ToType<double>(parameters[0][0]),
                        ATOOLS::ToType<double>(parameters[0][1]),
                        ATOOLS::ToType<int>(parameters[0][2]));
  }
  else {
    return NULL;
  }
}

void Blobdata_Getter::PrintInfo(std::ostream &str,const size_t width) const
{
  str<<"e.g. Blobdata  0.0  500.0  50  LinErr  Factorisation_Scale";
}


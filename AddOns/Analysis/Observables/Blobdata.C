#include "AddOns/Analysis/Observables/Primitive_Observable_Base.H"
#include "ATOOLS/Org/MyStrStream.H"
#include "ATOOLS/Org/Exception.H"
#include "ATOOLS/Org/Scoped_Settings.H"

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



DECLARE_GETTER(Blobdata,"Blobdata",
	       Primitive_Observable_Base,Analysis_Key);

Primitive_Observable_Base *
ATOOLS::Getter<Primitive_Observable_Base,Analysis_Key,Blobdata>::operator()(const Analysis_Key &key) const
{
  Scoped_Settings s{ key.m_settings };
  const auto parameters = s.SetDefault<std::string>({}).GetVector<std::string>();
  if (parameters.size()<5) return NULL;
  return new Blobdata(parameters[4],
                      HistogramType(parameters[3]),
                      s.Interprete<double>(parameters[0]),
                      s.Interprete<double>(parameters[1]),
                      s.Interprete<int>(parameters[2]));
}

void ATOOLS::Getter<Primitive_Observable_Base,Analysis_Key,Blobdata>::PrintInfo(std::ostream &str,const size_t width) const
{
  str<<"e.g. [Blobdata, 0.0, 500.0, 50, LinErr, Factorisation_Scale]";
}


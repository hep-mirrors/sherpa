#include "ATOOLS/Org/Exception.H"
#include "ATOOLS/Org/Message.H"
#include "PDF/Main/PDF_Base.H"

namespace PDF {
  class LHAPDF_Dummy : public PDF_Base {
  private:
  public:
    LHAPDF_Dummy()  {}
    ~LHAPDF_Dummy() {}
    PDF_Base * GetCopy() { return new LHAPDF_Dummy; }

    void   CalculateSpec(const double&,const double&) {}

    double GetXPDF(const ATOOLS::Flavour&) { return 0.; }
    double GetXPDF(const kf_code&, bool)   { return 0.; }

  };
}

using namespace PDF;
using namespace ATOOLS;

DECLARE_PDF_GETTER(LHAPDF_Getter);

PDF_Base *LHAPDF_Getter::operator()
  (const Parameter_Type &args) const
{
  return new LHAPDF_Dummy();
}

void LHAPDF_Getter::PrintInfo
(std::ostream &str,const size_t width) const
{
  str<<"LHAPDF Dummy";
}

extern "C" void InitPDFLib()
{
  THROW(fatal_error,"Sherpa not compiled with LHAPDF support.");
}

extern "C" void ExitPDFLib()
{
}

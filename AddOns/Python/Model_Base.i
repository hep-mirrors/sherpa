%{
#include "ATOOLS/Math/MyComplex.H"
#include "MODEL/Main/Model_Base.H"
#include "PDF/Main/PDF_Base.H"

using namespace MODEL;
%}
 
%include "std_complex.i"

namespace MODEL {

  class Model_Base {

  protected:

    virtual void InitVertices() = 0;
    virtual void ParticleInit() = 0;

  public:

    static void ShowSyntax(const size_t mode);
    void InitializeInteractionModel();
    virtual bool ModelInit(PDF::PDF_Base* pdf) = 0;

    std::string Name();

    virtual int                     ScalarNumber(const std::string);
    virtual double                  ScalarConstant(const std::string);
    virtual std::complex<double>    ComplexConstant(const std::string);
    virtual std::string             MappedLorentzName(const std::string& label);

  };

}

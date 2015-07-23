#ifndef SHERPA_Tools_Scale_Variations_H
#define SHERPA_Tools_Scale_Variations_H

#include <map>
#include <string>
#include <vector>
#include "ATOOLS/Org/MyStrStream.H"
#include "ATOOLS/Phys/Flavour.H"
#include "ATOOLS/Phys/Weight_Info.H"

namespace MODEL {
  class One_Running_AlphaS;
}

namespace PDF {
  class PDF_Base;
}

namespace PHASIC {
  class Process_Base;
}

namespace SHERPA {

  //! This class is a book-keeping entity for each scale variation
  class Scale_Variation {
  private:
    //! Switches whether the PDFs were self-created and have to be deleted
    bool m_deletepdfs,m_deleteas;
    //! Factors by which the nominal scales are to be varied by
    double m_muR2fac,m_muF2fac;
    //! Value
    double m_val;
    //! Value for RS events
    std::vector<double> m_RSvals;
    //! Pointers to the PDFs to be used
    PDF::PDF_Base *p_pdf1,*p_pdf2;
    //! Ids of the PDFs to be used
    int m_pdf1id,m_pdf2id;
    //! Names of the PDFs to be used
    std::string m_pdf1set,m_pdf2set;
    //! PDF set members to be used
    int m_pdf1setmember,m_pdf2setmember;
    //! Pointer to AlphaS to be used
    MODEL::One_Running_AlphaS *p_as;
    //! Name of the scale variation
    std::string m_name;

    /*!
     * @return Generated name of the scale variation
     *
     * This method generates the name of this scale variation
     */
    std::string GenerateName();

  public:
    /*!
     * @param muR2fac double Factor to vary muR2 by
     * @param muF2fac double Factor to vary muF2 by
     * @param pdf1 PDF_Base PDF to vary to for beam 1
     * @param pdf2 PDF_Base PDF to vary to for beam 2
     * @param as One_Running_AlphaS Corresponding AlphaS
     * @param deletepdfs bool Whether to delete the pdfs at the end or not
     * @param deleteas bool Whether to delete AlphaS at the end or not
     *
     * Constructor
     */
    Scale_Variation(const double &muR2fac, const double &muF2fac,
                    PDF::PDF_Base * pdf1, PDF::PDF_Base * pdf2,
                    MODEL::One_Running_AlphaS * as,
                    bool deletepdfs, bool deleteas);
    /*!
     * Destructor
     */
    ~Scale_Variation();

    /*!
     * @return Returns the variation factor for muR2
     *
     * Returns the variation factor for muR2
     */
    inline double MuR2Fac() const { return m_muR2fac; }
    /*!
     * @return Returns the variation factor for muF2
     *
     * Returns the variation factor for muF2
     */
    inline double MuF2Fac() const { return m_muF2fac; }

    /*!
     * @return Returns the value of the varied ME
     *
     * Returns the value of the varied ME
     */
    inline double Value() const                 { return m_val; }
    /*!
     * @param i size_t Index of the RS weight
     * @return Returns the value of the varied ME
     *
     * Returns the value of the varied RS ME
     */
    inline double Value(const size_t& i) const  { return m_RSvals[i]; }

    /*!
     * @return Returns the ID of the PDF of beam 1
     *
     * Returns the ID of the PDF of beam 1
     */
    inline int    Pdf1Id() const  { return m_pdf1id; }
    /*!
     * @return Returns the ID of the PDF of beam 2
     *
     * Returns the ID of the PDF of beam 2
     */
    inline int    Pdf2Id() const  { return m_pdf2id; }
    /*!
     * @return Returns the ID of the PDFs
     *
     * Returns the ID of the PDFs, if both are identical, -1 otherwise
     */
    inline int    PdfId() const
    { if (m_pdf1id==m_pdf2id) return m_pdf1id; else return -1; }

    /*!
     * @return Returns the PDF of beam 1
     *
     * Returns the PDF of beam 1
     */
    inline PDF::PDF_Base * PDF1() { return p_pdf1; }
    /*!
     * @return Returns the PDF of beam 2
     *
     * Returns the PDF of beam 2
     */
    inline PDF::PDF_Base * PDF2() { return p_pdf2; }

    /*!
     * @return Returns the AlphaS
     *
     * Returns the PDF of AlphaS
     */
    inline MODEL::One_Running_AlphaS * AlphaS() { return p_as; }

    /*!
     * @return Returns the name
     *
     * Returns the name of the scale variation
     */
    inline std::string Name() const { return m_name; }

    /*!
     * @param val double Value of the varied ME
     *
     * Sets the value of the varied ME
     */
    inline void SetValue(const double &val)                { m_val=val; }
    /*!
     * @param i size_t Index of RS ME
     * @param val double Value of the varied RS ME
     *
     * Sets the value of the varied RS ME
     */
    inline void SetValue(const size_t& i,const double &vl) { m_RSvals[i]=vl; }

    /*!
     * Deletes all stored RS ME values
     */
    inline void DeleteRSValues()                   { m_RSvals.clear(); }
    /*!
     * @param n size_t Number of RS ME elements
     *
     * Initialises the RS ME value storage vector
     */
    inline void InitialisRSValues(const size_t& n) { m_RSvals.resize(n,0.); }
  };

  typedef std::map<std::string, Scale_Variation *> NamedScaleVariationMap;

  std::ostream& operator<<(std::ostream &s,const Scale_Variation &sv);
  std::ostream& operator<<(std::ostream &s,const NamedScaleVariationMap &nsvm);
  std::ostream& operator<<(std::ostream &s,const NamedScaleVariationMap *nsvm);

  //! This class performs all a posteriori scale variations
  class Scale_Variations {
  private:
    //! Switches stearing the behaviour of the scale variations
    bool m_on,m_loadlhapdf,m_ckkw,m_kpnegativepdf;
    //! Generic quark and gluon flavours
    ATOOLS::Flavour m_gluon,m_quark;
    //! Map holding all Scale_Variation book-keeping objects
    NamedScaleVariationMap * p_nsvmap;

    //! This struct holds all relevant objects for simple access during
    //! the recomputation of the varied weights
    struct Parameters {
      double B,VI,KP,renterm,facterm1,facterm2;
      double muR2,muF12,muF22;
      double x1,x2,x1p,x2p;
      std::vector<double> rswgts,rsmuR2s,rsmuF2s;
      std::vector<double> renwgts,kpwgts;
      size_t oqcd,oew;
      int fl1,fl2;
      ATOOLS::mewgttype::code type;
      std::vector<ATOOLS::DADS_Info> dads;
      std::vector<ATOOLS::RDA_Info> rda;
      ATOOLS::Cluster_Sequence_Info csi;
    } m_params;

    /*!
     * This method resets values stored in all Scale_Variation objects and
     * Parameters
     */
    void ResetValues();
    /*!
     * @param winfo Weight_Info Weight information of the event
     * @param proc Process_Base The process that generated the event
     *
     * This method extracts all parameters necessary for the scale variations
     * from the event
     */
    void ExtractParameters(const ATOOLS::Weight_Info &winfo,
                           PHASIC::Process_Base * proc);
    /*!
     * @param sv Scale_Variation Which variation to calculate
     * @param proc Process_Base The process that produced the event
     * @return Returns true
     *
     * This method loops over all variations to be calculated and calls
     * the actual Calculate method for each
     */
    bool Calculate(Scale_Variation * sv,
                   PHASIC::Process_Base * proc);
    /*!
     * @param type mewgttype::code Type of the event
     * @param B double Value of the Born matrix element
     * @param VI double Value of the Virtual correction and the I operator
     * @param renwgts vector<double> Scale dependence coefficients of VI
     * @param kpwgts vector<double> Coefficients of the KP terms
     * @param dads vector<DADS_Info> DADS term info objects
     * @param rda vector<RDA_Info> RDA info objects
     * @param csi Cluster_Sequence_Info Cluster sequence info object
     * @param x1 double Partonic x of beam 1
     * @param x2 double Partonic x of beam 2
     * @param x1p double Partonic x' of beam 1 for KP terms
     * @param x2p double Partonic x' of beam 2 for KP terms
     * @param fl1 int Flavour of incoming parton of beam 1
     * @param fl2 int Flavour of incoming parton of beam 2
     * @param muR2 double Reference renormalisation scale
     * @param muF12 double Reference renormalisation scale of beam 1
     * @param muF22 double Reference renormalisation scale of beam 2
     * @param muR2fac double Variation factor for renormalisation scale
     * @param muF2fac double Variation factor for factorisation scale
     * @param oqcd size_t QCD order of the process computed
     * @param oew size_t EW order of the process computed
     * @param pdf1 PDF_Base PDF of beam 1
     * @param pdf2 PDF_Base PDF of beam 2
     * @param as One_Running_AlphaS AlphaS parametrisation
     * @return Returns the varied weight of the event
     *
     * This method does the actual calculation of the variation
     */
    double Calculate(const ATOOLS::mewgttype::code& type,
                     const double& B, const double& VI,
                     const std::vector<double>& renwgts,
                     const std::vector<double>& kpwgts,
                     const std::vector<ATOOLS::DADS_Info>& dads,
                     const std::vector<ATOOLS::RDA_Info>& rda,
                     const ATOOLS::Cluster_Sequence_Info& csi,
                     const double& x1, const double& x2,
                     const double& x1p, const double& x2p,
                     const int& fl1, const int& fl2,
                     const double& muR2,
                     const double& muF12, const double& muF22,
                     const double& muR2fac, const double& muF2fac,
                     const size_t& oqcd, const size_t& oew,
                     PDF::PDF_Base * pdf1, PDF::PDF_Base * pdf2,
                     MODEL::One_Running_AlphaS * as);
    //! Not implemented yet, will be used to simplify Calculate
    double ReweightLO();
    double ReweightBVI();
    double ReweightKP();
    double ReweightDADS();
    double ReweightNLOCT();
    double ReweightH();
    /*!
     * @param fa double Value of PDF of beam 1 at factorisation scale
     * @param fb double Value of PDF of beam 2 at factorisation scale
     * @param csi Cluster_Sequence_Info Cluster sequence info object
     * @param muF2fac doubel Variation factor of factorisation scale
     * @param pdf1 PDF_Base PDF of beam 1
     * @param pdf2 PDF_Base PDF of beam 2
     * @param skip bool If true, simply return 1.
     * @return Returns PDF ratio factor
     *
     * This method calculates the PDF ration factor for MEPS setups
     */
    double PDFRatioFactor(const double& fa, const double& fb,
                          const ATOOLS::Cluster_Sequence_Info& csi,
                          const double& muF2fac,
                          PDF::PDF_Base * pdf1, PDF::PDF_Base * pdf2,
                          bool skip=false);
  public:
    /*!
     * Constructor
     */
    Scale_Variations();
    /*!
     * Destructor
     */
    ~Scale_Variations();

    /*!
     * @param winfo Weight_Info Weight info of the generated event
     * @param proc Process_Base The process that produced the event
     *
     * This method computes the scale variations
     */
    bool ComputeVariations(const ATOOLS::Weight_Info &winfo,
                           PHASIC::Process_Base * proc);

    /*!
     * @return Named scale variations map
     *
     * This method returns the scale variations map
     */
    inline NamedScaleVariationMap * GetNamedScalesMap() const
    { return p_nsvmap; }

  };

  std::ostream& operator<<(std::ostream &s,const Scale_Variations &svs);
}
#endif

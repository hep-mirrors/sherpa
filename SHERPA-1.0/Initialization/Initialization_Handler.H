#ifndef Initialization_Handler_h
#define Initialization_Handler_h

#include "Model_Base.H"
#include "Beam_Spectra_Handler.H"
#include "ISR_Handler.H"
#include "Matrix_Element_Handler.H"
#include "Shower_Handler.H"
#include "Beam_Remnant_Handler.H"
#include "Fragmentation_Handler.H"
#include "Hadron_Decay_Handler.H"

#include "Data_Read.H"

#include <string>
#include <iostream> 

namespace SHERPA {
  class Initialization_Handler {
  private :
    std::string m_path,m_file;
    std::string m_modeldat, m_beamdat, m_isrdat, m_medat, m_showerdat;
    std::string m_beamremnantdat, m_fragmentationdat,m_hadrondecaysdat;

    double                 m_beam_splimits[2];
    double                 m_bunch_splimits[2];
    APHYTOOLS::Flavour     m_beam_particles[2];
    APHYTOOLS::Flavour     m_bunch_particles[2];

    AORGTOOLS::Data_Read * p_dataread;

    MODEL::Model_Base                  * p_model;
    BEAM::Beam_Spectra_Handler         * p_beamspectra;
    ISR::ISR_Handler                   * p_isrhandler;
    Matrix_Element_Handler             * p_mehandler;
    Shower_Handler                     * p_showerhandler;
    Beam_Remnant_Handler               * p_beamremnants;
    Fragmentation_Handler              * p_fragmentation;
    Hadron_Decay_Handler               * p_hadrondecays;
  private :
    bool InitializeTheModel();
    bool InitializeTheBeams();
    bool InitializeThePDFs();
    bool InitializeTheMatrixElements();
    bool InitializeTheShowers();
    bool InitializeTheBeamRemnants();
    bool InitializeTheFragmentation();
    bool InitializeTheHadronDecays();
    //bool InitializeTheHardDecays();
    //bool InitializeTheAnalysis();
  public :
    Initialization_Handler(std::string=std::string("./"),
			   std::string=std::string("Run.dat"));
    ~Initialization_Handler();

    bool InitializeTheFramework();
    bool CalculateTheHardProcesses();

    MODEL::Model_Base           * GetModel()                { return p_model;           }
    ISR::ISR_Handler            * GetISRHandler()           { return p_isrhandler;      }
    BEAM::Beam_Spectra_Handler  * GetBeamSpectraHandler()   { return p_beamspectra;     }
    Matrix_Element_Handler      * GetMatrixElementHandler() { return p_mehandler;       }
    Shower_Handler              * GetShowerHandler()        { return p_showerhandler;   }
    Beam_Remnant_Handler        * GetBeamRemnantHandler()   { return p_beamremnants;    }
    Fragmentation_Handler       * GetFragmentationHandler() { return p_fragmentation;   }
    Hadron_Decay_Handler        * GetHadronDecayHandler()   { return p_hadrondecays;    }
    APHYTOOLS::Flavour          * GetBeams()                { return m_beam_particles;  }
    APHYTOOLS::Flavour          * GetBunches()              { return m_bunch_particles; }
  };
}


#endif

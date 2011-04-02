#ifndef AMISIC_Main_MI_Base_H
#define AMISIC_Main_MI_Base_H
/*!
  \file MI_Base.H
  \brief Declares the class MI_Base
*/

#include "ATOOLS/Org/CXXFLAGS.H"
#include "ATOOLS/Org/File_IO_Base.H"
#include "ATOOLS/Phys/Particle_List.H"
#include "ATOOLS/Phys/Blob_List.H"
#include "ATOOLS/Org/Exception.H"
#include <map>

namespace PHASIC { class Process_Base; }

namespace AMISIC {

  class MI_Base: public ATOOLS::File_IO_Base, 
		 public ATOOLS::Terminator_Object {
  public:

    enum TypeID { HardEvent =  1,
		  SoftEvent =  2,
		  Unknown   = 99 };

    typedef std::pair<const std::string,MI_Base*> String_MI_Base_Pair;
    typedef std::map<const std::string,MI_Base*> String_MI_Base_Map;

  private:

    std::string m_name;
    TypeID      m_type;
    
    bool m_jetveto;

  protected:

    static String_MI_Base_Map *s_bases;

    double       *m_start, *m_stop, *m_last;
    unsigned int  m_nparameter;

    bool m_generatedparameter, m_generatedprocess, m_filledblob;

    ATOOLS::Particle_List m_inparticles, m_outparticles;

    PHASIC::Process_Base *p_xs;
    static MI_Base       *s_hard, *s_soft;

    static bool s_stophard, s_stopsoft, s_cleaned;

  public:

    // constructor
    MI_Base(std::string name,TypeID type,unsigned int nparameter=0,
	    unsigned int infiles=1,unsigned int outfiles=1);
    /*!
      \fn MI_Base(std::string name,TypeID type,unsigned int nparameter=0,
      unsigned int infiles=1,unsigned int outfiles=1)
      \brief The standard constructor. Each MI_Base informs the system 
      about its name, its type, the number of its parameters
      and the number of input and output files.
     */

    // destructor
    virtual ~MI_Base();

    // member functions
    static void  CleanUp();
    /*!
      \fn static void  CleanUp()
      \brief Cleans up all MI_Bases.
    */
    virtual bool Initialize();
    /*!
      \fn virtual bool Initialize()
      \brief Initializes the MI_Base.
    */

    virtual void Update(const MI_Base *mibase);
    /*!
      \fn virtual void Update(const MI_Base *mibase)
      \brief Updates the MI_Base with information 
      from another MI_Base.
     */
    virtual void Reset();
    /*!
      \fn virtual void Reset()
      \brief Resets the MI_Base before event generation.
    */

    virtual bool VetoProcess(ATOOLS::Blob *blob);
    /*!
      \fn virtual bool VetoProcess(ATOOLS::Blob *blob)
      \brief Vetoes a given process according to the 
      probability density used for event generation.
    */
    virtual bool GenerateProcess();
    /*!
      \fn virtual bool GenerateProcess()
      \brief Generates one process according to the 
      specific probability density.
    */

    bool FillBlob(ATOOLS::Blob *blob);
    /*!
      \fn bool FillBlob(ATOOLS::Blob *blob)
      \brief Converts the information from the event 
      generation step into the Sherpa standard format.
    */

    static bool StopGeneration(TypeID type=Unknown);
    /*!
      \fn static bool StopGeneration(TypeID type=Unknown)
      \brief System flag. Indicates that no further interactions
      of the given type shall be generated.
    */

    static void UpdateAll(const MI_Base *mibase);
    /*!
      \fn static void UpdateAll(const MI_Base *mibase)
      \brief Updates all active MI_Bases with the information from
      the given MI_Base.
    */
    static void ResetAll();  
    /*!
      \fn static void ResetAll()
      \brief Resets all MI_Bases.
    */

    static std::string TypeToString(TypeID type);
    static TypeID StringToType(std::string type);

    static void SetStopGeneration(TypeID type=Unknown,const bool stop=false);
    /*!
      \fn static void SetStopGeneration(TypeID type=Unknown,
      const bool stop=false)
      \brief Setter method to force stopping the event generation.
    */

    // inline functions
    inline void SetStart(const double start,unsigned int i=0)
    { if (i<m_nparameter) m_start[i]=start; }
    inline void SetStop(const double stop,unsigned int i=0)
    { if (i<m_nparameter) m_stop[i]=stop; }

    inline void SetJetVeto(const bool jv) { m_jetveto=jv; }

    inline double Start(unsigned int i=0) const { return m_start[i]; }
    inline double Stop(unsigned int i=0) const  { return m_stop[i];  }

    inline TypeID Type() const { return m_type; }

    inline bool JetVeto() const { return m_jetveto; }

    inline PHASIC::Process_Base *XS() const { return p_xs; }

  }; // end of class MI_Base
  /*!
    \class MI_Base
    \brief The base class for all underlying event simulation classes

    This class is the base class for all underlying event 
    simulation classes. It is an interface to the overall steering of 
    the package, which is provided by the class Amisic.<br>
    It keeps track of the initialized models and checks the consistency 
    of the underlying event generation procedure. During the 
    event generation step it handles the filling of underlying event 
    blobs and updates the generator.
  */

  class MI_None: public MI_Base {
  protected:

    void Update(const MI_Base *mibase);
    /*!
      \fn void Update(const MI_Base *mibase)
      \brief Dummy method
    */

  public:

    // constructor
    MI_None(TypeID type);
    /*!
      \fn MI_None(TypeID type)
      \brief The standard constructor.
    */

    // destructor
    ~MI_None();

    // member functions
    bool Initialize();
    /*!
      \fn bool Initialize()
      \brief Dummy method
    */
    void Reset();
    /*!
      \fn bool Reset()
      \brief Dummy method
    */

    bool GenerateOrderingParameter();
    /*!
      \fn bool GenerateOrderingParameter()
      \brief Dummy method
    */
    bool VetoProcess(ATOOLS::Blob *blob);
    /*!
      \fn bool VetoProcess(ATOOLS::Blob *blob)
      \brief Dummy method
    */
    bool GenerateProcess();
    /*!
      \fn bool GenerateProcess()
      \brief Dummy method
    */

  }; // end of class MI_None
  /*!
    \class MI_None
    \brief Dummy class to represent the absence of underlying events.

    This class represents the absence of underlying events. 
    It is employed if either no underlying events are 
    to be generated at all or if e.g. only hard underlying events 
    are desired and no soft underlying event model is to be used.
  */

} // end of namespace AMISIC

#endif

#ifndef MI_Base_H
#define MI_Base_H

#include "File_IO_Base.H"
#include "Blob_List.H"
#include "Exception.H"
#include <map>

#ifdef USING__Sherpa
namespace SHERPA { class Matrix_Element_Handler; }
#endif

namespace EXTRAXS { class XS_Base; }

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

  protected:

    static String_MI_Base_Map s_bases;

    double       *m_start, *m_stop, *m_last;
    unsigned int  m_nparameter;

    bool m_dicedparameter, m_dicedprocess, m_filledblob, m_weighted;

    ATOOLS::Blob *p_blob;

#ifdef USING__Sherpa
    SHERPA::Matrix_Element_Handler *p_mehandler;
#endif

    EXTRAXS::XS_Base *p_xs;

    static bool s_stophard, s_stopsoft, s_cleaned;

  public:

    // constructor
    MI_Base(std::string _m_name,TypeID _m_type,
	    unsigned int _m_nparameter=0,
	    unsigned int infiles=1,unsigned int outfiles=1);

    // destructor
    virtual ~MI_Base();

    // member functions
    static void  CleanUp();
    virtual bool Initialize();

    virtual void Update(const MI_Base *mibase);
    virtual void Reset();

    virtual bool DiceOrderingParameter();
    virtual bool DiceProcess();

    bool CreateBlob(ATOOLS::Blob *blob);

    static bool StopGeneration(TypeID type=Unknown);

    static void UpdateAll(const MI_Base *mibase);
    static void ResetAll();  

    static std::string TypeToString(TypeID type);
    static TypeID StringToType(std::string type);

    // inline functions
    inline void SetStart(const double _m_start,unsigned int i=0)
    { if (i<m_nparameter) m_start[i]=_m_start; }
    inline void SetStop(const double _m_stop,unsigned int i=0)
    { if (i<m_nparameter) m_stop[i]=_m_stop; }

    inline double Start(unsigned int i=0) const { return m_start[i]; }
    inline double Stop(unsigned int i=0) const  { return m_stop[i];  }

#ifdef USING__Sherpa
    inline SHERPA::Matrix_Element_Handler *MEHandler() 
    { return p_mehandler; }
#endif

    inline EXTRAXS::XS_Base *XS() const { return p_xs; }

  }; // end of class MI_Base

  class MI_None: public MI_Base {
  protected:

    void Update(const MI_Base *mibase);

  public:

    // constructor
    MI_None(TypeID _m_type);

    // destructor
    ~MI_None();

    // member functions
    bool Initialize();
    void Reset();

    bool DiceOrderingParameter();
    bool DiceProcess();

  }; // end of class MI_None

} // end of namespace AMISIC

#endif

#ifndef Grid_Handler_H
#define Grid_Handler_H

#include "Data_To_Function.H"
#include "Type.H"

namespace AMISIC {

  template <class Argument_Type,class Result_Type>
  class Grid_Handler {
  public:

    typedef Argument_Type GridArgumentType;
    typedef Result_Type   GridResultType;

    typedef ATOOLS::Data_To_Function<GridArgumentType,GridResultType> 
    GridFunctionType;

  private:

    GridFunctionType *p_grid;

    ATOOLS::Type::ID m_streamtype;

    std::string m_streamname, m_datatag;

    void Init();

  protected:

    bool ReadFromFile(std::string tempname);
    bool WriteToFile(std::string tempname,
		     std::vector<std::string> &comments);

  public:

    // constructors
    Grid_Handler();
    Grid_Handler(ATOOLS::Type::ID streamtype,std::string streamname);
    Grid_Handler(GridFunctionType *_p_grid);

    // destructor
    ~Grid_Handler();

    // member functions
    bool ReadIn(ATOOLS::Type::ID temptype=ATOOLS::Type::TUnknown,
		std::string tempname=std::string());
    bool WriteOut(ATOOLS::Type::ID temptype=ATOOLS::Type::TUnknown,
		  std::string tempname=std::string(),std::vector<std::string> 
		  comments=std::vector<std::string>(0));

    // inline functions
    inline ATOOLS::Type::ID   StreamType() const { return m_streamtype; }
    inline const std::string &StreamName() const { return m_streamname; }
    inline const std::string &DataTag() const    { return m_datatag;    }

    inline const std::vector<std::string> &Comments() const 
    { return m_comments; }

    inline void SetStreamType(const ATOOLS::Type::ID streamtype) 
    { m_streamtype=streamtype; }
    inline void SetStreamName(const std::string streamname)
    { m_streamname=streamname; }
    inline void SetDataTag(const std::string datatag)
    { m_datatag=datatag; }

    inline GridFunctionType *Grid() { return p_grid; }

    inline GridResultType operator()(GridArgumentType x) 
    { return (*p_grid)(x); }
    
    inline GridArgumentType operator[](GridResultType y)
    { return (*p_grid)[y]; }

  }; // end of Class Grid_Handler

} // end of namespace AMISIC

#include "Grid_Handler.C"

#endif

#ifndef Grid_Handler_H
#define Grid_Handler_H

#include "Data_To_Function.H"
#include "Grid_Creator_Base.H"
#include "Type.H"

namespace AMISIC { template <class Argument_Type,class Result_Type>
		   class Grid_Creator_Base;                              } 

namespace AMISIC {

  template <class Argument_Type,class Result_Type>
  class Grid_Handler {
  public:

    typedef Argument_Type GridArgumentType;
    typedef Result_Type   GridResultType;

    typedef ATOOLS::Data_To_Function<GridArgumentType,GridResultType> GridFunctionType;
    typedef Grid_Creator_Base<GridArgumentType,GridResultType> GridCreatorType;

  private:

    GridFunctionType *p_grid;
    GridCreatorType *p_gridcreator;

    ATOOLS::Type::ID m_streamtype;
    std::string m_streamname;

    std::string m_datatag;

    void Init();

  protected:

    bool ReadFromFile(std::string tempname);
    bool WriteToFile(std::string tempname,std::vector<std::string>& comments);

  public:

    // constructors
    Grid_Handler();
    Grid_Handler(ATOOLS::Type::ID _m_streamtype,std::string _m_streamname);
    Grid_Handler(GridFunctionType *_p_grid);

    // destructor
    ~Grid_Handler();

    // setters
    void SetStreamType(const ATOOLS::Type::ID _m_streamtype);
    void SetStreamName(const std::string _m_streamname);
    void SetDataTag(const std::string _m_datatag);

    // getters
    const ATOOLS::Type::ID         StreamType() const;
    const std::string              StreamName() const;
    const std::string              DataTag() const;
    const std::vector<std::string> Comment() const;

    // member functions
    bool ReadIn(ATOOLS::Type::ID temptype=ATOOLS::Type::TUnknown,std::string tempname=std::string());
    bool WriteOut(ATOOLS::Type::ID temptype=ATOOLS::Type::TUnknown,std::string tempname=std::string(),
		  std::vector<std::string> comments=std::vector<std::string>(0));

    GridFunctionType *Grid();

    GridResultType   operator()(GridArgumentType x);
    GridArgumentType operator[](GridResultType y);

  }; // end of Class Grid_Handler

} // end of namespace AMISIC

#include "Grid_Handler.inl.H"
#include "Grid_Handler.C"

#endif

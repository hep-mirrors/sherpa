#ifndef Data_To_Function_H
#define Data_To_Function_H

#include "Axis.H"

#include <map>
#include <vector>

namespace ATOOLS {

  template <class Argument_Type,class Result_Type>
  class Data_To_Function {
  public:

    enum MonotonyID { 
      Increasing =  1,
      Decreasing =  2,
      None       =  3,
      MUnknown    = 99 
    };

    enum AcquisitionModeID { 
      LowerData     =  2,
      Data          =  3,
      UpperData     =  4,
      Interpolation =  5,
      AUnknown       = 99 
    };

    enum InterpolationModeID { 
      Linear        =  1,
      IUnknown       = 99 
    };

  public:

    typedef Argument_Type ArgumentType;
    typedef Result_Type   ResultType;

    typedef std::vector<ArgumentType> ArgumentVector;
    typedef std::vector<ResultType>   ResultVector;

    typedef std::pair<ArgumentType,ResultType> XYPair;
    typedef std::pair<ResultType,ArgumentType> YXPair;

    typedef std::vector<XYPair> XYVector;
    typedef std::vector<YXPair> YXVector;

    typedef Data_To_Function<Argument_Type,ResultType> DataToFunctionType;

  private:

    ATOOLS::Axis<ArgumentType> *p_xaxis;
    ATOOLS::Axis<ResultType>   *p_yaxis;

    XYVector *p_xydata;
    YXVector *p_yxdata;

    MonotonyID m_monotony;

    AcquisitionModeID   m_acquisitionmode;
    InterpolationModeID m_interpolationmode;

    void Init();
    void Resize(unsigned int newsize);

    inline bool Sort()
    { bool sortx=SortX(), sorty=SortY(); return sortx&&sorty; }
    bool SortX();
    bool SortY();

    unsigned int ClosestX(Argument_Type x,
			  unsigned int &left,unsigned int &right);
    unsigned int ClosestY(Result_Type y,
			  unsigned int &left,unsigned int &right);

    ResultType   LinearY(ArgumentType x,
			 unsigned int &left,unsigned int &right);
    ArgumentType LinearX(ResultType y,
			 unsigned int &left,unsigned int &right);

  protected:

    unsigned int ClosestX(Argument_Type x);
    unsigned int ClosestY(Result_Type y);

    ResultType   InterpolateY(ArgumentType x);
    ArgumentType InterpolateX(ResultType y);

    ResultType   DataY(ArgumentType x,AcquisitionModeID tempmode);
    ArgumentType DataX(ResultType y,AcquisitionModeID tempmode);

  public:

    // constructors
    Data_To_Function();
    Data_To_Function(const DataToFunctionType &reference,
		     bool adoptscaling=true);
    Data_To_Function(const ArgumentVector *_m_xdata,
		     const ResultVector *_m_ydata);

    // destructor
    ~Data_To_Function();

    // setters
    inline void SetMonotony(const MonotonyID _m_monotony)
    { if (m_monotony==MUnknown) m_monotony = _m_monotony; }
    inline void SetAcquisitionMode(const AcquisitionModeID _m_acquisitionmode)
    { m_acquisitionmode = _m_acquisitionmode; }
    inline void SetInterpolationMode(const InterpolationModeID _m_interpolationmode)
    { m_interpolationmode = _m_interpolationmode; }

    // getters
    inline const MonotonyID          Monotony() const
    { return m_monotony; }
    inline const AcquisitionModeID   AcquisitionMode() const
    { return m_acquisitionmode; }
    inline const InterpolationModeID InterpolationMode() const
    { return m_interpolationmode; }

    inline ArgumentType XMin()
    { return p_xydata->size()>0?(*p_xaxis)[(*p_xydata).front().first]:0.0; }
    inline ArgumentType XMax()
    { return p_xydata->size()>0?(*p_xaxis)[(*p_xydata).back().first]:0.0; }
    inline ResultType   YMin()
    { return p_yxdata->size()>0?(*p_yaxis)[(*p_yxdata).front().first]:0.0; }
    inline ResultType   YMax()
    { return p_yxdata->size()>0?(*p_yaxis)[(*p_yxdata).back().first]:0.0; }

    inline ArgumentType DeltaXMin()
    { ResultType dummyleft, dummyright;
      return DeltaXMin(dummyleft,dummyright); }
    inline ArgumentType DeltaXMax()
    { ResultType dummyleft, dummyright; 
      return DeltaXMax(dummyleft,dummyright); }
    inline ResultType   DeltaYMin()
    { ArgumentType dummyleft, dummyright; 
      return DeltaYMin(dummyleft,dummyright); }
    inline ResultType   DeltaYMax()
    { ArgumentType dummyleft, dummyright; 
      return DeltaYMax(dummyleft,dummyright); }
    
    ArgumentType DeltaXMin(ResultType& left,ResultType& right);
    ArgumentType DeltaXMax(ResultType& left,ResultType& right);
    ResultType   DeltaYMin(ArgumentType& left,ArgumentType& right);
    ResultType   DeltaYMax(ArgumentType& left,ArgumentType& right);
    
    // access functions
    inline ATOOLS::Axis<ArgumentType> *XAxis() { return p_xaxis; }
    inline ATOOLS::Axis<ResultType>   *YAxis() { return p_yaxis; }

    // member functions
    bool Import(const DataToFunctionType &reference,bool adoptscaling=true);
    bool Import(const ArgumentVector *_p_xdata,
		const ResultVector *_p_ydata,bool normal=true);

    void Export(ArgumentVector *_p_xdata,ResultVector *_p_ydata,
		bool normal=true);

    bool AddPoint(ArgumentType x,ResultType y);
    bool DeleteXPoint(ArgumentType x);
    bool DeleteYPoint(ResultType y);
    bool ReplaceXPoint(ArgumentType x,ResultType y);
    bool ReplaceYPoint(ArgumentType x,ResultType y);

    void RescaleY(const std::string scalename);
    void RescaleX(const std::string scalename);

    inline void Clear()
    { p_xydata->clear(); p_yxdata->clear(); }

    inline ResultType   YData(unsigned int i)
    { return (*p_yaxis)[(*p_yxdata)[i%p_yxdata->size()].first]; }
    inline ArgumentType XData(unsigned int i)
    { return (*p_xaxis)[(*p_xydata)[i%p_xydata->size()].first]; }

    inline XYPair XYData(unsigned int i)
    { return XYPair((*p_xaxis)[(*p_xydata)[i%p_xydata->size()].first],
		    (*p_yaxis)[(*p_xydata)[i%p_xydata->size()].second]); }
    inline YXPair YXData(unsigned int i)
    { return YXPair((*p_yaxis)[(*p_yxdata)[i%p_yxdata->size()].first],
		    (*p_xaxis)[(*p_yxdata)[i%p_yxdata->size()].second]); }

    inline unsigned int YDataSize() { return p_yxdata->size(); }
    inline unsigned int XDataSize() { return p_xydata->size(); }

    int YPosition(ResultType y);
    int XPosition(ArgumentType x);

    void ScaleY(ResultType scalefactor);
    void ScaleX(ArgumentType scalefactor);

    void MoveY(ResultType distance);
    void MoveX(ArgumentType distance);

    void NormalizeY(ArgumentType xmin=0.0,ArgumentType xmax=0.0,
		    ResultType norm=1.0);
    void NormalizeX(ResultType ymin=0.0,ResultType ymax=0.0,
		    ArgumentType norm=1.0);

    DataToFunctionType *IntegralY(ArgumentType xmin=0.0,ArgumentType xmax=0.0,
				  std::string xscaling=nullstring,
				  std::string yscaling=nullstring,
				  bool forward=true,const MonotonyID=None);
    DataToFunctionType *IntegralX(ResultType ymin=0.0,ResultType ymax=0.0,
  				  std::string yscaling=nullstring,
				  std::string xscaling=nullstring,
				  bool forward=true,const MonotonyID=None);

    ResultType   IntegrateY(ArgumentType xmin=0.0,ArgumentType xmax=0.0);
    ArgumentType IntegrateX(ResultType ymin=0.0,ResultType ymax=0.0);

    DataToFunctionType *DerivativeY(ArgumentType xmin=0.0,
				    ArgumentType xmax=0.0,
				    std::string xscaling="",
				    std::string yscaling="");
    DataToFunctionType *DerivativeX(ResultType ymin=0.0,
				    ResultType ymax=0.0,
				    std::string xscaling="",
				    std::string yscaling="");

    ResultType   DifferentialY(ArgumentType x);
    ArgumentType DifferentialX(ResultType y);
    
    ResultType   Y(ArgumentType x,AcquisitionModeID tempmode=AUnknown);
    ArgumentType X(ResultType y,AcquisitionModeID tempmode=AUnknown);

    inline ResultType   operator()(ArgumentType x)
      { return Y(x,m_acquisitionmode); }
    inline ArgumentType operator[](ResultType y)
      { return X(y,m_acquisitionmode); }

  }; // end of class Data_To_Function

} // end of namespace ATOOLS

#endif

#ifndef Data_To_Function_H
#define Data_To_Function_H

#include "Axis.H"

#include "My_Functional.H"
#include "My_Map.H"
#include "My_Vector.H"

namespace ATOOLS {

  template <class Argument_Type,class Result_Type>
  class Data_To_Function: public std::unary_function<Argument_Type,Result_Type> {
  public:

    enum AcquisitionModeID { // Analytic      =  1,
                             Data          =  2,
			     Interpolation =  3,
			     // Fit           =  4,
			     AUnknown       = 99 };

    enum InterpolationModeID { Linear        =  1,
			       IUnknown       = 99 };

  public:

    typedef Argument_Type ArgumentType;
    typedef Result_Type   ResultType;

    typedef std::vector<ArgumentType> ArgumentVector;
    typedef std::vector<ResultType>   ResultVector;

#ifdef __GNUC__
#if __GNUC__>2
    typedef __gnu_cxx::__normal_iterator<ArgumentType*, std::vector<ArgumentType> > ArgumentIterator;
    typedef __gnu_cxx::__normal_iterator<ResultType*, std::vector<ResultType> > ResultIterator;
#else
    typedef ArgumentType *ArgumentIterator;
    typedef ResultType *ResultIterator;
#endif
#else
    typedef ArgumentType *ArgumentIterator;
    typedef ResultType *ResultIterator;
#endif

    typedef std::pair<ArgumentType,ResultType> XYPair;
    typedef std::pair<ResultType,ArgumentType> YXPair;

    typedef std::vector<XYPair> XYVector;
    typedef std::vector<YXPair> YXVector;

    typedef std::map<ArgumentType,Result_Type> XYMap;
    typedef std::map<ResultType,Argument_Type> YXMap;

#ifdef __GNUC__
#if __GNUC__>2
    typedef __gnu_cxx::__normal_iterator<XYPair*, std::vector<XYPair> > XYVectorIterator;
    typedef __gnu_cxx::__normal_iterator<YXPair*, std::vector<YXPair> > YXVectorIterator;
#else
    typedef XYPair *XYVectorIterator;
    typedef YXPair *YXVectorIterator;
#endif
#else
    typedef XYPair *XYVectorIterator;
    typedef YXPair *YXVectorIterator;
#endif
    typedef std::_Rb_tree_iterator<XYPair,XYPair&,XYPair*> XYMapIterator;
    typedef std::_Rb_tree_iterator<YXPair,YXPair&,YXPair*> YXMapIterator;

    typedef Data_To_Function<Argument_Type,ResultType> DataToFunctionType;

  private:

    ATOOLS::Axis<ArgumentType> *p_xaxis;
    ATOOLS::Axis<ResultType>   *p_yaxis;

    XYVector *p_xydata;
    YXVector *p_yxdata;

    AcquisitionModeID   m_acquisitionmode;
    InterpolationModeID m_interpolationmode;

    void Init();

    void Resize(unsigned int newsize);

    void Sort();
    void SortX();
    void SortY();

    unsigned int ClosestX(Argument_Type x,unsigned int &left,unsigned int &right);
    unsigned int ClosestY(Result_Type y,unsigned int &left,unsigned int &right);

    ResultType   LinearY(ArgumentType x,unsigned int &left,unsigned int &right);
    ArgumentType LinearX(ResultType y,unsigned int &left,unsigned int &right);

  protected:

    unsigned int ClosestX(Argument_Type x);
    unsigned int ClosestY(Result_Type y);

    ResultType   InterpolateY(ArgumentType x);
    ArgumentType InterpolateX(ResultType y);

    ResultType   DataY(ArgumentType x);
    ArgumentType DataX(ResultType y);

  public:

    // constructors
    Data_To_Function();
    Data_To_Function(const Data_To_Function<ArgumentType,ResultType> &reference);
    Data_To_Function(XYMap *_m_xydata);
    Data_To_Function(ArgumentVector *_m_xdata,ResultVector *_m_ydata);

    // destructor
    ~Data_To_Function();

    // setters
    void SetAcquisitionMode(AcquisitionModeID _m_acquisitionmode);
    void SetInterpolationMode(InterpolationModeID _m_interpolationmode);

    // getters
    AcquisitionModeID   AcquisitionMode();
    InterpolationModeID InterpolationMode();

    ArgumentType XMin();
    ArgumentType XMax();
    ResultType   YMin();
    ResultType   YMax();

    ArgumentType DeltaXMin();
    ArgumentType DeltaXMax();
    ResultType   DeltaYMin();
    ResultType   DeltaYMax();

    ArgumentType DeltaXMin(ResultType& left,ResultType& right);
    ArgumentType DeltaXMax(ResultType& left,ResultType& right);
    ResultType   DeltaYMin(ArgumentType& left,ArgumentType& right);
    ResultType   DeltaYMax(ArgumentType& left,ArgumentType& right);

    // access functions
    ATOOLS::Axis<ArgumentType> *XAxis();
    ATOOLS::Axis<ResultType>   *YAxis();

    // member functions
    void Import(Data_To_Function<ArgumentType,ResultType> &reference);
    void Import(XYMap *_p_xydata,bool normal=true);
    void Import(ArgumentVector *_p_xdata,ResultVector *_p_ydata,bool normal=true);

    void Export(XYMap *_p_xydata,bool normal=true);
    void Export(ArgumentVector *_p_xdata,ResultVector *_p_ydata,bool normal=true);

    bool AddPoint(ArgumentType x,ResultType y);
    bool DeleteXPoint(ArgumentType x);
    bool DeleteYPoint(ResultType y);
    bool ReplaceXPoint(ArgumentType x,ResultType y);
    bool ReplaceYPoint(ArgumentType x,ResultType y);

    ResultType   YData(unsigned int i);
    ArgumentType XData(unsigned int i);

    XYPair XYData(unsigned int i);
    YXPair YXData(unsigned int i);

    unsigned int YDataSize();
    unsigned int XDataSize();

    int YPosition(ResultType y);
    int XPosition(ArgumentType x);

    void ScaleY(ResultType scalefactor);
    void ScaleX(ArgumentType scalefactor);

    void MoveY(ResultType distance);
    void MoveX(ArgumentType distance);

    void NormalizeY(ArgumentType xmin=0.0,ArgumentType xmax=0.0,ResultType norm=1.0);
    void NormalizeX(ResultType ymin=0.0,ResultType ymax=0.0,ArgumentType norm=1.0);

    DataToFunctionType *IntegralY(ArgumentType xmin=0.0,ArgumentType xmax=0.0);
    DataToFunctionType *IntegralX(ResultType ymin=0.0,ResultType ymax=0.0);

    ResultType   IntegrateY(ArgumentType xmin,ArgumentType xmax);
    ArgumentType IntegrateX(ResultType ymin,ResultType ymax);

    ResultType   Y(ArgumentType x,AcquisitionModeID tempmode = AUnknown);
    ArgumentType X(ResultType y,AcquisitionModeID tempmode = AUnknown);

    ResultType   operator()(ArgumentType x);
    ArgumentType operator[](ResultType y);

  }; // end of class Data_To_Function

} // end of namespace ATOOLS

#include "Data_To_Function.IC"
#include "Data_To_Function.C"

#endif

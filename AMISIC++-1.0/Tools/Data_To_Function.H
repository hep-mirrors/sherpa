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

    enum AcquisitionMode { // Analytic      =  1,
	                   Data          =  2,
			   Interpolation =  3,
			   // Fit           =  4,
			   AUnknown       = 99 };

    enum InterpolationMode { Linear        =  1,
			     IUnknown       = 99 };

  public:

    typedef Argument_Type ArgumentType;
    typedef Result_Type   ResultType;

    typedef std::vector<ArgumentType> ArgumentVector;
    typedef std::vector<ResultType>   ResultVector;

    typedef std::pair<const ArgumentType,ResultType> XYPair;
    typedef std::pair<const ResultType,ArgumentType> YXPair;

    typedef std::map<ArgumentType,Result_Type> XYDataMap;
    typedef std::map<ResultType,Argument_Type> YXDataMap;

#ifdef __GNUC__
#if __GNUC__>2
    typedef __gnu_cxx::__normal_iterator<ArgumentType*, std::vector<ArgumentType> > ArgumentIterator;
    typedef __gnu_cxx::__normal_iterator<ResultType*, std::vector<ResultType> >     ResultIterator;
#else
    typedef ArgumentType *ArgumentIterator;
    typedef ResultType *ResultIterator;
#endif
#else
    typedef ArgumentType *ArgumentIterator;
    typedef ResultType *ResultIterator;
#endif
    typedef std::_Rb_tree_iterator<XYPair,XYPair&,XYPair*> XYDataIterator;
    typedef std::_Rb_tree_iterator<YXPair,YXPair&,YXPair*> YXDataIterator;

  private:

    ATOOLS::Axis<ArgumentType> *p_xaxis;
    ATOOLS::Axis<ResultType>   *p_yaxis;

    ArgumentVector *p_xdata;
    ResultVector   *p_ydata;

    XYDataMap      *p_xydata;
    YXDataMap      *p_yxdata;

    AcquisitionMode   m_acquisitionmode;
    InterpolationMode m_interpolationmode;

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
    Data_To_Function(XYDataMap *_m_xydata);
    Data_To_Function(ArgumentVector *_m_xdata,ResultVector *_m_ydata);

    // destructor
    ~Data_To_Function();

    // setters
    void SetAcquisitionMode(AcquisitionMode _m_acquisitionmode);
    void SetInterpolationMode(InterpolationMode _m_interpolationmode);

    // getters
    AcquisitionMode   GetAcquisitionMode();
    InterpolationMode GetInterpolationMode();

    ArgumentType GetXMin();
    ArgumentType GetXMax();
    ResultType   GetYMin();
    ResultType   GetYMax();

    ArgumentType GetDeltaXMin();
    ArgumentType GetDeltaXMax();
    ResultType   GetDeltaYMin();
    ResultType   GetDeltaYMax();

    ArgumentType GetDeltaXMin(ResultType& left,ResultType& right);
    ArgumentType GetDeltaXMax(ResultType& left,ResultType& right);
    ResultType   GetDeltaYMin(ArgumentType& left,ArgumentType& right);
    ResultType   GetDeltaYMax(ArgumentType& left,ArgumentType& right);

    // access functions
    ATOOLS::Axis<ArgumentType> *XAxis();
    ATOOLS::Axis<ResultType>   *YAxis();

    // member functions
    void Import(XYDataMap *_p_xydata);
    void Import(ArgumentVector *_p_xdata,ResultVector *_p_ydata);
    void Import(Data_To_Function<ArgumentType,ResultType> &reference);

    void Export(XYDataMap *_p_xydata);
    void Export(ArgumentVector *_p_xdata,ResultVector *_p_ydata);

    bool AddPoint(ArgumentType x,ResultType y);
    bool DeleteXPoint(ArgumentType x);
    bool DeleteYPoint(ResultType y);
    bool ReplaceXPoint(ArgumentType x,ResultType y);
    bool ReplaceYPoint(ArgumentType x,ResultType y);

    ResultType   YData(unsigned int i);
    ArgumentType XData(unsigned int i);

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

    ResultType   IntegrateY(ArgumentType xmin,ArgumentType xmax);
    ArgumentType IntegrateX(ResultType ymin,ResultType ymax);

    ResultType   Y(ArgumentType x,AcquisitionMode tempmode = AUnknown);
    ArgumentType X(ResultType y,AcquisitionMode tempmode = AUnknown);

    ResultType   operator()(ArgumentType x);
    ArgumentType operator[](ResultType y);

  }; // end of class Data_To_Function

} // end of namespace ATOOLS

#include "Data_To_Function.IC"
#include "Data_To_Function.C"

#endif

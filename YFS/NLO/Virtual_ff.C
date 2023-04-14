#include "YFS/NLO/Virtual_ff.H"
#include "YFS/Main/YFS_Handler.H"

#include "ATOOLS/Org/Exception.H"  
#include "ATOOLS/Org/Message.H"  

using namespace YFS;


Virtual_ff::Virtual_ff(int _order):
m_order(_order)
{   
    m_orderMin = 0;
    m_orderMax = 3;
}

Virtual_ff::~Virtual_ff(){

}

bool Virtual_ff::CheckOrder(){
    if(m_order < m_orderMin){
        msg_Error()<<"YFS pragmatic order = "<<m_order
        <<std::endl;
        THROW(fatal_error,"YFS Pragmatic order can't be negative");
    }
    if(m_order > m_orderMax){
        msg_Error()<<"YFS pragmatic order = "<<m_order
        <<std::endl;
        THROW(fatal_error,"YFS not implemented at this order yet");
    }
    return true;
}

double Virtual_ff::Calculate(double beta,int order){
    m_beta  = beta;
    m_order = order;
    if(CheckOrder()){
        switch(m_order){
            case 0: m_virtual =  0.;
            case 1: m_virtual = Del1();
            case 2: m_virtual = Del1() + Del2();
            case 3: m_virtual = Del1() + Del2() + Del3();
        }
    }
    if(IsBad(m_virtual)){
        msg_Error()<<"YFS m_virtual = "<<m_virtual<<std::endl;
    }
    msg_Debugging()<<METHOD<<"Correction Order "<<m_order<<std::endl
                    <<"Virtual weight = "<<m_virtual<<std::endl;
    return m_virtual;
}

double Virtual_ff::Del1(){
    m_del1 = 0.5*m_beta;
    if(IsBad(m_del1)){
        THROW(fatal_error,"Error in YFS Virtual delta1");
    }
    else return m_del1;
}


double Virtual_ff::Del2(){
    m_del2 = 0.125*m_beta*m_beta;
    if(IsBad(m_del2)){
        THROW(fatal_error,"Error in YFS Virtual delta2");
    }
    else return m_del2;
}


double Virtual_ff::Del3(){
    m_del3 = 1./48.*pow(m_beta,3);
    if(IsBad(m_del3)){
        THROW(fatal_error,"Error in YFS Virtual delta3");
    }
    else return m_del3;
}
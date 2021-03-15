#include "fvCFD.H"


//test primitive type: label, scalar
using namespace Foam;

template<class T>
void printTraits()
{
    Info<< pTraits<T>::typeName
        << ": zero=" << pTraits<T>::zero
        << " one=" << pTraits<T>::one << endl;
}


int main(int argc,char* argv[]) {

    //label: A label is an int32_t or int64_t as specified by the pre-processor macro WM_LABEL_SIZE.
    Info<<"labelMin: "<<labelMin<<endl;
    Info<<"labelMax: "<<labelMax<<endl;

    IOobject::writeDivider(Info);
    //scalar: Single floating point number identical to float or double depending on whether WM_SP, WM_DP or WM_LP is defined.
    Info<<"GREAT: "<<GREAT<<endl; //浮点误差的倒数
    Info<<"ROOTGREAT: "<<ROOTGREAT<<endl; //浮点误差的倒数 的根号
    Info<<"VGREAT: "<<VGREAT<<endl; // 能表示的最大的数字/10
    Info<<"SMALL: "<<SMALL<<endl; //返回目标数据类型能表示的最逼近1的正数和1的差的绝对值,浮点误差
    Info<<"ROOTSMALL: "<<ROOTSMALL<<endl;
    Info<<"VSMALL: "<<VSMALL<<endl; // 能表示的最小的数字
    Info<<"ROOTVSMALL: "<<ROOTVSMALL<<endl; // 能表示的最小的数字 // 能表示的最小的数字

    IOobject::writeDivider(Info);
    //zero: A class representing the concept of 0 used to avoid unnecessary manipulations for objects that are known to be zero at compile-time.
    //Zero: Global zero

    Info<<"Zero: "<<static_cast<label>(Zero)<<endl;
    label index=5-Zero;
    Info<<"index: "<<index<<endl;

    //PTraits： Traits class for primitives， All primitives need a specialised version of this class. The specialised version will normally also require a conversion method
    
    printTraits<bool>();
    printTraits<label>();
    printTraits<scalar>();
    printTraits<vector>();
    printTraits<tensor>();

    


    return 0;
}
#include "fvCFD.H"


using namespace Foam;
int main(int argc,char* argv[]) {

    //tmp 就像c++ stl中的shared_ptr,但是其指针，最多只能两个tmp共享
    //tmp的设计主要用于函数对象返回值，当函数返回时，返回的对象常常需要两次拷贝构造，使用tmp，可以减少大对象的构造开销，而仅是构造tmp对象
    //相比autoPtr,tmp不仅仅可以有指针构造，也可以是引用构造，并且autoPtr的指针，不能由2个autoPtr共享，autoPtr像c++ stl中的unique_ptr

    scalarField f0(1.0,2.0);
    tmp<scalarField> tf0(f0);
    Info<<"typeName: "<<tf0.typeName()<<endl;
    Info<<"is Tmp: "<<tf0.isTmp()<<endl; // 如果tmp是通过ref构造，则不是IsTmp，而是CONST_REF
    Info<<"is empty: "<<tf0.empty()<<endl; //如果tmp是通过ref构造，是empty
    Info<<"is valid: "<<tf0.valid()<<endl;


    IOobject::writeDivider(Info);
    scalarField* f1=new scalarField(1.0, 2.0);
    tmp<scalarField> tf1(f1);

    // Info<<"is Tmp: "<<tf1.isTmp()<<endl;
    Info<<"is Empty before calling ref(): "<<tf1.empty()<<endl; //是.empty()来判断tmp是否为空，如果是->empty则是call scalarField的empty
    scalarField f1Ref=tf1.ref(); //发生了赋值构造
    Info<<"is Empty after calling ref(): "<<tf1.empty()<<endl;
    Info<<"is valid after calling ref(): "<<tf1.valid()<<endl;


    IOobject::writeDivider(Info);
    Info<<"is Empty before calling operator(): "<<tf1.empty()<<endl; //是.empty()来判断tmp是否为空，如果是->empty则是call scalarField的empty
    const scalarField& f1Ref2=tf1();
    Info<<"is Empty after calling operator(): "<<tf1.empty()<<endl;
    Info<<"is valid after calling operator(): "<<tf1.valid()<<endl;


    IOobject::writeDivider(Info);

    scalarField* f2=new scalarField(1.0, 2.0);
    tmp<scalarField> tf2(f2);
     Info<<"is Empty before calling ptr(): "<<tf2->empty()<<endl;
    tmp<scalarField> f2Transfer(tf2.ptr());
    Info<<"is Tmp after calling ptr(): "<<tf2.isTmp()<<endl;
    Info<<"is Empty after calling ptr(): "<<tf2.empty()<<endl; //call ptr, tmp就被置空了
    Info<<"is valid after calling ptr(): "<<tf2.valid()<<endl;

    IOobject::writeDivider(Info);

    tmp<scalarField> tf3(new scalarField(1.0, 2.0));
    tmp<scalarField> tf4(tf3);
    Info<<"is Empty tf3: "<<tf3.empty()<<endl; 
    Info<<"is Empty tf4: "<<tf4.empty()<<endl; 
    Info<<"is unique tf4: "<<tf4().unique()<<endl; //not unique
    Info<<"ref count: "<<tf4().count()<<endl; //1

    //tmp<scalarField> tf5(tf3);  //FOAM FATAL ERROR: Attempt to create more than 2 tmp's referring to the same object of type tmp<N4Foam5FieldIdEE>
    //最多两个tmp指向某个对象





    // delete f1;

    return 0;
}
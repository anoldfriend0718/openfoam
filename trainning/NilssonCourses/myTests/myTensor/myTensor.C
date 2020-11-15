#include "fvCFD.H"

using namespace Foam;
int main() {
    tensor t1 {11,12,13,21,22,23,31,32,33};
    Info<<"T1: "<<t1<<endl;
    Info<<"T1 xx: "<<t1.xx()<<endl;
    Info<<"T1 xy: "<<t1.xy()<<endl;
    Info<<"T1 xz: "<<t1.xz()<<endl;
    Info<<"T1 x:  "<<t1.x()<<endl;
    Info<<"T1 y:  "<<t1.y()<<endl;
    Info<<"T1 z:  "<<t1.z()<<endl;
    Info<<"T1.T():  "<<t1.T()<<endl;
    Info<<"det(T1):  "<<det(t1)<<endl;

    tensor t2 {1,2,3,4,5,6,7,8,9};
    Info<<"T2: "<<t2<<endl;

    Info<<"T1+T2: "<<t1+t2<<endl;
    scalar r1=t1&&t2;
    Info<<"T1:T2: "<<r1<<endl;

    dimensionedTensor sigma("sigma",dimensionSet(1,-1,-2,0,0,0,0),t2);
    Info<<"sigma: "<<sigma<<endl;
    
    Info<< "Sigma name: " << sigma.name() << endl;
    Info<< "Sigma dimensions: " << sigma.dimensions() << endl;
    Info<< "Sigma value: " << sigma.value() << endl;
    Info<< "Sigma yy: "<<sigma.value().yy()<<endl;

    tensorField tf1(2,tensor::one);
    Info<<"tf1: "<<tf1<<endl;

    tensorField tf2(2,t2);
    Info<<"tf2: "<<tf2<<endl;
    Info<<"2*tf2: "<<2*tf2<<endl;

    tensorField tf3{tf2};
    tf3[0]=tensor{1,1,1,1,1,1,1,1,1};
    Info<<"tf3: "<<tf3<<endl;






    return 0;
}
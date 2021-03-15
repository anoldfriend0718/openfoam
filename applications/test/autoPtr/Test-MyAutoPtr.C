#include "fvCFD.H"


using namespace Foam;
int main(int argc,char* argv[]) {    
    autoPtr<scalarField> f1(new scalarField(1.0, 2.0));
    Info<<"is empty of f1: "<<f1.empty()<<endl;
    autoPtr<scalarField> f2(f1);
    Info<<"is empty of f1: "<<f1.empty()<<endl;
    Info<<"is empty of f2: "<<f2.empty()<<endl;

    IOobject::writeDivider(Info);
   
    return 0;
}
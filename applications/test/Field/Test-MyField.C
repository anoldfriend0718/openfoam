#include "fvCFD.H"
#include "scalarField.H"

using namespace Foam;
int main() {
    scalarField field1({1,2,0});
    Info<<"field1: "<<field1<<endl;

    scalar minValue=1e-3;
    tmp<scalarField> tField2=max(field1,minValue);
    const scalarField& field2=tField2.ref();
    Info<<"field2: "<<field2<<endl;

    return 0;

}
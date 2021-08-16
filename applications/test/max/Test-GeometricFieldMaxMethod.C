#include "dimensionedScalarFwd.H"
#include "fvCFD.H"
#include "dimensionSets.H"
#include "linear.H"
#include "surfaceInterpolation.H"
#include "volFieldsFwd.H"
#include "zero.H"

using namespace Foam;
int main(int argc, char *argv[]) 
{

    #include "setRootCaseLists.H"

    #include "createTime.H"
    #include "createMesh.H"



    //Test faceInterpolated

    volScalarField T
    (
        IOobject(
            "T",
            runTime.timeName(),
            mesh,
            IOobject::NO_READ,
            IOobject::NO_WRITE
        ),
        mesh,
        dimensionedScalar(dimTemperature,Zero)
    );

    forAll(T, i)
    {
        T[i]=i;
    }
    T[0]=-0.1;
    T[1]=-0.2;


    Info<<"T: "<<T.internalField().field()<<endl;

    Info<<"apply max method"<<endl;
    // T.max(0.0);
    T=max(T,dimensionedScalar(dimTemperature,Zero));
    Info<<"T: "<<T.internalField().field()<<endl;

    IOobject::writeDivider(Info);

    return 0;

}
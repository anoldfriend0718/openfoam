#include "dimensionedScalarFwd.H"
#include "fvCFD.H"
#include "scalar.H"
#include "volFieldsFwd.H"
#include "dimensionSets.H"
#include "specialSurfaceArea.H"

using namespace Foam;
int main(int argc, char *argv[]) 
{
    #include "setRootCaseLists.H"
    #include "createTime.H"
    #include "createMesh.H"


    volScalarField eps
    (
        IOobject
        (
            "eps",
            runTime.timeName(),
            mesh,
            IOobject::MUST_READ,
            IOobject::AUTO_WRITE
        ),
        mesh
    );
    volScalarField coke
    (
        IOobject
        (
            "coke",
            runTime.timeName(),
            mesh,
            IOobject::MUST_READ,
            IOobject::AUTO_WRITE
        ),
        mesh
    );

    specialSurfaceArea surfaceArea(mesh);
    runTime++;

    surfaceArea.correct();

    volScalarField SS
    (
        IOobject
        (
            "SS",
            runTime.timeName(),
            mesh,
            IOobject::NO_READ,
            IOobject::AUTO_WRITE
        ),
        mesh,
        dimensionedScalar("SSi",dimless/dimLength,Zero)
    );
    SS.primitiveFieldRef()=surfaceArea.SS().field();

    runTime.write();

    
    return 0;


    









}

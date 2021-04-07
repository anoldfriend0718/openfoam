#include "fvCFD.H"
#include "volFieldsFwd.H"

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

    volScalarField epsGrad //without filter depth
    (
        "epsGrad",
        mag(fvc::grad(eps))
    );

    epsGrad.write();

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

    volScalarField cokeGrad //without filter depth
    (
        "cokeGrad",
        mag(fvc::grad(coke))
    );

    cokeGrad.write();

    volScalarField speciesSurfaceArea
    {
        "speciesSurfaceArea",
        2.0*mag(fvc::grad(eps))*(1-(1-coke)*(1-coke))
    };
    speciesSurfaceArea.write();

    volScalarField surfaceArea
    {
        "surfaceArea",
        speciesSurfaceArea*Foam::pow(1e-6,1)
    }; 

    surfaceArea.write();

    volScalarField filterDepth
    {
        "filterDepth",
        1-(1-coke)*(1-coke)
    };
    filterDepth.write();






    









}

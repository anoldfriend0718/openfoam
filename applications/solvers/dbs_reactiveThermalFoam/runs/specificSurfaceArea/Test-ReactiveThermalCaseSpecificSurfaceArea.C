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

    volScalarField speciesSurfaceArea
    {
        "speciesSurfaceArea",
        2*mag(fvc::grad(1-eps))*(1-eps*eps)
    };

    volScalarField surfaceArea
    {
        "surfaceArea",
        speciesSurfaceArea*Foam::pow(1e-6,1)
    }; 

    surfaceArea.write();





    









}

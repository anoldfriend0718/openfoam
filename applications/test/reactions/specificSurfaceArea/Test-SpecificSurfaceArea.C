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


    volScalarField epsGrad //without filter depth
    (
        "epsGrad",
        mag(fvc::grad(eps))
    );

    epsGrad.write();

    volScalarField speciesSurfaceArea1
    {
        "speciesSurfaceArea1",
        mag(fvc::grad(eps))*4*eps*(1-eps)
    };

    speciesSurfaceArea1.write();

    volScalarField speciesSurfaceArea2
    {
        "speciesSurfaceArea2",
        2*mag(fvc::grad(eps))*(1-eps*eps)
    };

    speciesSurfaceArea2.write();








    // epsSpeciesSurfaceArea.write();


    









}

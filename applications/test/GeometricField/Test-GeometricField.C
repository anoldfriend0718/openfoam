#include "fvCFD.H"
#include "dimensionSets.H"
#include "volFieldsFwd.H"

using namespace Foam;
int main(int argc, char *argv[]) 
{

    #include "setRootCaseLists.H"

    #include "createTime.H"
    #include "createMesh.H"

    volScalarField eps
    (
        IOobject(
            "eps",
            runTime.timeName(),
            mesh,
            IOobject::NO_READ,
            IOobject::NO_WRITE
        ),
        mesh,
        dimensionedScalar(dimless,0.0)
    );

    forAll(eps, i)
    {
        eps[i]=min(i*0.1,1.0);
    }

    volScalarField tEps(max(eps,0.001));
    // eps=tEps;

    // Info<<"eps: "<<eps<<endl;

    volScalarField rho
    (
        IOobject(
            "rho",
            runTime.timeName(),
            mesh,
            IOobject::NO_READ,
            IOobject::NO_WRITE
        ),
        mesh,
        dimensionedScalar(dimDensity,1.0)
    );

    volScalarField rhoEps("rhoEpsEps",rho*tEps*tEps);
    Info<<"rho*eps*eps: "<<rhoEps<<endl;


    return 0;

}
#include "fvCFD.H"
#include "dimensionSets.H"
#include "linear.H"
#include "surfaceInterpolation.H"
#include "volFieldsFwd.H"

using namespace Foam;
int main(int argc, char *argv[]) 
{

    #include "setRootCaseLists.H"

    #include "createTime.H"
    #include "createMesh.H"


    // //Test operator *
    // volScalarField eps
    // (
    //     IOobject(
    //         "eps",
    //         runTime.timeName(),
    //         mesh,
    //         IOobject::NO_READ,
    //         IOobject::NO_WRITE
    //     ),
    //     mesh,
    //     dimensionedScalar(dimless,0.0)
    // );

    // forAll(eps, i)
    // {
    //     eps[i]=min(i*0.1,1.0);
    // }

    // volScalarField tEps(max(eps,0.001));
    // // eps=tEps;

    // // Info<<"eps: "<<eps<<endl;

    // volScalarField rho
    // (
    //     IOobject(
    //         "rho",
    //         runTime.timeName(),
    //         mesh,
    //         IOobject::NO_READ,
    //         IOobject::NO_WRITE
    //     ),
    //     mesh,
    //     dimensionedScalar(dimDensity,1.0)
    // );

    // volScalarField rhoEps("rhoEpsEps",rho*tEps*tEps);
    // Info<<"rho*eps*eps: "<<rhoEps<<endl;

    // //Test operator gSumProd
    // volScalarField Yi
    // (
    //     IOobject(
    //         "Yi",
    //         runTime.timeName(),
    //         mesh,
    //         IOobject::NO_READ,
    //         IOobject::NO_WRITE
    //     ),
    //     mesh,
    //     dimensionedScalar(dimDensity,2.0)
    // );

    // Info<<"gSumProd Yi: "<<gSumProd(Yi,Yi)<<endl;

    // volVectorField U
    // (
    //     IOobject(
    //         "U",
    //         runTime.timeName(),
    //         mesh,
    //         IOobject::NO_READ,
    //         IOobject::NO_WRITE
    //     ),
    //     mesh,
    //     dimensionedVector(dimVelocity,vector(2.0,1.0,0))
    // );

    // Info<<"gSumCompProd U: "<<gSumCmptProd(U,U)<<endl;



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

    surfaceScalarField surfT
    (
        IOobject
        (
            "surfT",
            runTime.timeName(),
            mesh,
            IOobject::NO_READ,
            IOobject::NO_WRITE
        ),
        linearInterpolate(T)
    );

    // T.write();
    // surfT.write();
    Info<<"T: "<<endl;
    Info<<T<<endl;

    IOobject::writeDivider(Info);
    Info<<"surfT: "<<endl;
    Info<<surfT<<endl;






    return 0;

}
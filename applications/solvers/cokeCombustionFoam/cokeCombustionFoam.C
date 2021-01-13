#include "IOdictionary.H"
#include "IFstream.H"
#include "Time.H"
#include "argList.H"
#include "fvMesh.H"
#include "rhoReactionThermo.H"
#include "thermo.H"
#include "volFieldsFwd.H"
#include "constSolidThermo.H"
#include "cokeCombustion.H"
#include "fvc.H"

using namespace Foam;

int main(int argc, char *argv[]) {

    //#include "setRootCase.H"
    argList args(argc, argv);
    if (!args.checkRootCase()) {
        Foam::FatalError.exit();
    }
    //#include "createTime.H"
    Foam::Time runTime(Foam::Time::controlDictName, args);
    //#include "createMesh.H"
    fvMesh mesh(Foam::IOobject(fvMesh::defaultRegion, runTime.timeName(), runTime, IOobject::MUST_READ));

    // Info << "Reading gas thermophysical properties\n" << endl;
    // autoPtr<rhoReactionThermo> pThermo(rhoReactionThermo::New(mesh));
    // rhoReactionThermo& thermo = pThermo();
    // thermo.validate(args.executable(), "h", "e");
    

    // Info << "Reading solid thermophysical properties\n" << endl;
    // IOdictionary thermoSolid
    // (
    //     IOobject
    //     (
    //         "thermo.Solid",
    //         runTime.constant(),
    //         mesh,
    //         IOobject::MUST_READ_IF_MODIFIED,
    //         IOobject::NO_WRITE
    //     )
    // );
    // const constSolidThermo cokeThermo(thermoSolid.subDict("coke"));
    // const constSolidThermo rockThermo(thermoSolid.subDict("rock"));
    // Info<<"coke thermo: "<<cokeThermo<<endl;
    // Info<<"rock thermo: "<<rockThermo<<endl;
    
    // Info<<"read U and create phi..."<<endl;
    // volVectorField U
    // (
    //     IOobject
    //     (
    //         "U",
    //         runTime.timeName(),
    //         mesh,
    //         IOobject::MUST_READ,
    //         IOobject::AUTO_WRITE
    //     ),
    //     mesh
    // );

    // volScalarField rhof
    // (
    //     IOobject
    //     (
    //         "rho",
    //         runTime.timeName(),
    //         mesh,
    //         IOobject::NO_READ,
    //         IOobject::AUTO_WRITE
    //     ),
    //     thermo.rho()
    // );

    // surfaceScalarField phi
    // (
    //     IOobject
    //     (
    //         "phi",
    //         runTime.timeName(),
    //         mesh,
    //         IOobject::READ_IF_PRESENT,
    //         IOobject::AUTO_WRITE
    //     ),
    //     linearInterpolate(rhof*U) & mesh.Sf()
    // );


    //   //read eps, coke, rock fractions
    // volScalarField eps
    // (
    //     IOobject
    //     (
    //         "eps",
    //         runTime.timeName(),
    //         mesh,
    //         IOobject::MUST_READ,
    //         IOobject::AUTO_WRITE
    //     ),
    //     mesh
    // );
    // volScalarField coke
    // (
    //     IOobject
    //     (
    //         "coke",
    //         runTime.timeName(),
    //         mesh,
    //         IOobject::MUST_READ,
    //         IOobject::AUTO_WRITE
    //     ),
    //     mesh
    // );
    // volScalarField rock
    // (
    //     IOobject
    //     (
    //         "rock",
    //         runTime.timeName(),
    //         mesh,
    //         IOobject::NO_READ,
    //         IOobject::NO_WRITE
    //     ),
    //     1-eps-coke
    // );
    // Info<<"fluid porosity: "<<coke.field()<<", solid fraction: "<<(1-eps)->field()<<endl;
    // Info<<"coke fraction: "<<coke.field()<<", rock fraction: "<<rock.field()<<endl;

    // cokeCombustion combustion(mesh,thermo,cokeThermo,rockThermo);

    // for(label i=0;i<2;i++)
    // {
    //     Info<<endl<<"new step..."<<endl;
    //     combustion.correct();

    //     Info<<"Test the RR..."<<endl;
    //     Info<<"Time-splitting RR O2: "<<combustion.RRO2().field()<<endl;
    //     Info<<"Time-splitting RR CO2: "<<combustion.RRCO2().field()<<endl;
    //     Info<<"Time-splitting RR coke: "<<combustion.RRCoke().field()<<endl;
    //     Info<<"Normal RR O2: "<<combustion.calculateTransientRRO2()().field()<<endl;

    //     Info<<"Test the su matrix..."<<endl;
    //     volScalarField& YO2=thermo.composition().Y(0);
    //     volScalarField& YCO2=thermo.composition().Y(2);
    //     Info<<"R O2: "<<(combustion.R(YO2)&YO2)->field()<<endl;
    //     Info<<"R CO2: "<<(combustion.R(YCO2)&YCO2)->field()<<endl;
    //     Info<<"R coke: "<<(combustion.Rs(coke)&coke)->field()<<endl;

    //     Info<<"Test Qdot..."<<endl;
    //     Info<<"Qdot: "<<combustion.Qdot()->field()<<endl;

    //     runTime++;

    //     Info<<"simulate solving solid equation by update the eps and coke..."<<endl;
    //     eps+=0.1;
    //     coke-=0.1;

    //     Info<<"simulate solving energy equation by update the enthalpy..."<<endl;
    //     // Info<<"Cp: "<<thermo.Cp()<<endl;
    //     volScalarField& he=thermo.he();
    //     scalar he1=1200000;
    //     forAll(he, celli)
    //     {
    //         he[celli]=he1;
    //     }

    //     volScalarField::Boundary& heBf = he.boundaryFieldRef();
    //     forAll(heBf, patchi)
    //     {
    //         scalarField& hep=heBf[patchi];
    //         forAll(hep,facei)
    //         {
    //             hep[facei]=he1;
    //         }
    //     }

    
    //     thermo.correct();

        
    //    Info<<"simulate solving pressure equation by update the pressure..."<<endl;
    //    volScalarField& p=thermo.p();
    //    p*=1.05;

    //    Info<<"simulate solving the species equations by updating Y..."<<endl;
    //    YO2-=0.03;
    //    YCO2 +=0.03;



        

        
    // }



    return 0;
}
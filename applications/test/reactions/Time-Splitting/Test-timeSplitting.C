#include "IOdictionary.H"
#include "IFstream.H"
#include "Time.H"
#include "argList.H"
#include "fvMesh.H"
#include "fvmSup.H"
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

    Info << "Reading gas thermophysical properties\n" << endl;
    autoPtr<rhoReactionThermo> pThermo(rhoReactionThermo::New(mesh));
    rhoReactionThermo& thermo = pThermo();
    thermo.validate(args.executable(), "h", "e");
    

    Info << "Reading solid thermophysical properties\n" << endl;
    IOdictionary thermoSolid
    (
        IOobject
        (
            "thermo.Solid",
            runTime.constant(),
            mesh,
            IOobject::MUST_READ_IF_MODIFIED,
            IOobject::NO_WRITE
        )
    );
    const constSolidThermo cokeThermo(thermoSolid.subDict("coke"));
    const constSolidThermo rockThermo(thermoSolid.subDict("rock"));
    Info<<"coke thermo: "<<cokeThermo<<endl;
    Info<<"rock thermo: "<<rockThermo<<endl;
    
    Info<<"read U and create phi..."<<endl;
    volVectorField U
    (
        IOobject
        (
            "U",
            runTime.timeName(),
            mesh,
            IOobject::MUST_READ,
            IOobject::AUTO_WRITE
        ),
        mesh
    );

    volScalarField rhof
    (
        IOobject
        (
            "rho",
            runTime.timeName(),
            mesh,
            IOobject::NO_READ,
            IOobject::AUTO_WRITE
        ),
        thermo.rho()
    );

    surfaceScalarField phi
    (
        IOobject
        (
            "phi",
            runTime.timeName(),
            mesh,
            IOobject::READ_IF_PRESENT,
            IOobject::AUTO_WRITE
        ),
        linearInterpolate(rhof*U) & mesh.Sf()
    );


      //read eps, coke, rock fractions
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
    volScalarField rock
    (
        IOobject
        (
            "rock",
            runTime.timeName(),
            mesh,
            IOobject::NO_READ,
            IOobject::NO_WRITE
        ),
        1-eps-coke
    );
    Info<<"fluid porosity: "<<eps.field()<<", solid fraction: "<<(1-eps)->field()<<endl;
    Info<<"coke fraction: "<<coke.field()<<", rock fraction: "<<rock.field()<<endl;

    cokeCombustion combustion(mesh,thermo,cokeThermo,rockThermo);
    combustion.correct();

    Info<<"Test the RR..."<<endl;
    Info<<"Time-splitting RR O2: "<<combustion.RRO2().field()<<endl;
    Info<<"Time-splitting RR CO2: "<<combustion.RRCO2().field()<<endl;
    Info<<"Time-splitting RR coke: "<<combustion.RRCoke().field()<<endl;
    Info<<"Normal RR O2: "<<combustion.calculateTransientRRO2()().field()<<endl;

    
    
    
    Info<<"Test the su matrix..."<<endl;
    volScalarField YO2=thermo.composition().Y(0);
    volScalarField YCO2=thermo.composition().Y(2);
    Info<<"R O2: "<<(combustion.R(YO2)&YO2)->field()<<endl;
    Info<<"R CO2: "<<(combustion.R(YCO2)&YCO2)->field()<<endl;
    Info<<"R coke: "<<(combustion.Rs(coke)&coke)->field()<<endl;


    Info<<"R O2 source term: "<<combustion.R(YO2)->source()<<endl;
    Info<<"R O2 matrix: "<<combustion.R(YO2)()<<endl;

    volScalarField RO2_1=(combustion.R(YO2)&YO2)();
    auto RO2Matrix_1=fvm::Su(RO2_1,YO2);
    Info<<"R O2 matrix 1: "<<RO2Matrix_1()<<endl;



    Info<<"Test Qdot..."<<endl;
    Info<<"Qdot: "<<combustion.Qdot()->field()<<endl;

    
    volScalarField rhoCoke
    (
        IOobject
        (
            "solid",
            runTime.timeName(),
            mesh,
            IOobject::NO_READ,
            IOobject::AUTO_WRITE
        ),
        mesh,
        cokeThermo.density
    );

    const tmp<GeometricField<double, fvPatchField, volMesh>>& deltarRho
        =1./rhof-1.0/rhoCoke;
    Info<<"deltarRho<< "<<deltarRho->field()<<endl;
    Info<<"rhoCoke: "<<rhoCoke.field()<<endl;
    Info<<"rhof: "<<rhof.field()<<endl;

    Info<<"Rs source term: "<<combustion.Rs(coke)->source()<<endl;
    const tmp<fvMatrix<double>>& RRg=-deltarRho*rhof*combustion.Rs(coke);
    Info<<"RRg source term: "<<RRg->source()<<endl;

    Info<<"RRg: "<<RRg.ref()<<endl;

    Info<<"RRg2: "<<fvm::Su(RRg&coke,rhof)()<<endl;


    return 0;
}
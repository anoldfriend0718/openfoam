#include "fvCFD.H"
#include "PtrList.H"
#include "UList.H"
#include "dictionary.H"
#include "IFstream.H"
#include "Time.H"
#include "argList.H"
#include "fvMesh.H"
#include "rhoReactionThermo.H"
#include "specie.H"
#include "surfaceFieldsFwd.H"
#include "volFieldsFwd.H"
#include "reactingMixture.H"
#include "perfectGas.H"
#include "sensibleEnthalpy.H"
#include "thermo.H"
#include "thermodynamicConstants.H"
#include "constTransport.H"
#include "hConstThermo.H"
#include "turbulentFluidThermoModel.H"
#include "CombustionModel.H"

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

    Info << "Reading thermophysical properties\n" << endl;
    autoPtr<rhoReactionThermo> pThermo(rhoReactionThermo::New(mesh));
    rhoReactionThermo& thermo = pThermo();
    thermo.validate(args.executable(), "h", "e");

    Info << "Base class typeName: " << thermo.typeName << endl;
    Info << "type: " << thermo.type() << endl;
    Info << "thermal name: " << thermo.thermoName() << endl;

    Info<<"read the composition mass fraction..."<<endl;
    const basicSpecieMixture& compositions=thermo.composition();
    Info<<"species"<<compositions.species()<<endl;

    PtrList<volScalarField> Y = compositions.Y();
    forAll(Y,i)
    {
        Info<<Y[i].name()<<" is active: "<<compositions.active(i)
            <<" , mole weight: "<<compositions.Wi(i)
            <<" , field: "<<Y[i].primitiveField()<<endl;
    }

    typedef constTransport<species::thermo<hConstThermo<perfectGas<specie>>,sensibleEnthalpy>> gasThermoType;
    typedef reactingMixture<gasThermoType> reactingMixtureType;
    Info<< "check gas thermal type: "<<gasThermoType::typeName()<<endl;
    
    const reactingMixtureType& mixture=static_cast<const reactingMixtureType&>(compositions);
    Info<<"speciesCompositionTable: "<<mixture.specieComposition()<<endl;
    //get the thermo models of different composition in the mixture
    const PtrList<gasThermoType>& gasThermos=mixture.speciesData();

    Info<<"composition size: "<<gasThermos.size()<<endl;
    forAll(gasThermos, i)
    {
        Info<<"compostion ["<<i<<"]: "<<gasThermos[i].name()<<endl;
    }

    scalar p0=1e5;
    scalar T0=300;
    Info<<"p0: "<<p0<<", T0: "<<T0<<endl;
    Info<<"Pstd: "<<Pstd<<", Tstd: "<<Tstd<<endl;

    forAll(gasThermos, i)
    {
        Info<<gasThermos[i].name()
            <<", Cp: "<<compositions.Cp(i, p0, T0)
            <<", Cv: "<<compositions.Cv(i, p0, T0)
            <<", Ha: "<<compositions.Ha(i, p0, T0)
            <<", Hc: "<<compositions.Hc(i)
            <<", Hs: "<<compositions.Hs(i, p0, T0)
            <<", rho: "<<compositions.rho(i, p0,T0)
            <<", kappa: "<<compositions.kappa(i, p0, T0)
            <<", mu: "<<compositions.mu(i, p0, T0)
            <<", alphah: "<<compositions.alphah(i, p0, T0)
            <<endl;
    }

    Info<<"Mixture: "
        <<", Cp: "<<thermo.Cp()().field()
        <<", Cv: "<<thermo.Cv()().field()
        <<", Cpv: "<<thermo.Cpv()().field()
        <<", CpByCpv: "<<thermo.CpByCpv()().field()
        <<", He: "<<thermo.he()().field()
        <<", rho: "<<thermo.rho()().field()
        <<", psi: "<<thermo.psi()().field()
        <<", kappa: "<<thermo.kappa()().field()
        <<", mu: "<<thermo.mu()().field()
        <<", alphahe: "<<thermo.alphahe()().field()
        <<endl;

      const word inertSpecie(thermo.lookup("inertSpecie"));
    if (!compositions.species().found(inertSpecie))
    {
        FatalIOErrorIn(args.executable().c_str(), thermo)
            << "Inert specie " << inertSpecie << " not found in available species "
            << compositions.species()
            << exit(FatalIOError);
    }
    Info<<"inertSpecies: "<<inertSpecie<<endl;

    volScalarField rho
    (
        IOobject
        (
            "rho",
            runTime.timeName(),
            mesh
        ),
        thermo.rho()
    );

    Info<< "Reading field U\n" << endl;
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


    // volScalarField& p = thermo.p();

    Info<< "Reading/calculating face flux field phi\n" << endl;

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
        linearInterpolate(rho*U) & mesh.Sf()
    );

    Info << "Creating turbulence model.\n" << nl;
    autoPtr<compressible::turbulenceModel> turbulence
    (
        compressible::turbulenceModel::New
        (
            rho,
            U,
            phi,
            thermo
        )
    );


    Info<< "Creating reaction model\n" << endl;
    autoPtr<CombustionModel<rhoReactionThermo>> reaction
    (
        CombustionModel<rhoReactionThermo>::New(thermo, turbulence())
    );
    reaction->correct();


    PtrList<volScalarField> reactionRates;
    reactionRates.setSize(Y.size());
    forAll(Y,i)
    {
        reactionRates.set
        (
            i,
            new volScalarField
            (
                    IOobject
                    (
                        "R_" + Y[i].name(),
                        runTime.timeName(),
                        mesh,
                        IOobject::NO_READ,
                        IOobject::AUTO_WRITE
                    ),
                    mesh,
			        dimensionedScalar("Ri", dimensionSet(1, -3, -1, 0, 0, 0, 0), 0)
            )
        );
        reactionRates[i] = reaction->R(Y[i]) & Y[i]; 
        Info<<"reaction rate of Y, "<<Y[i].name()<<" :"<<reactionRates[i].field()<<endl;
    }








    

    // const auto& thermoTypes=compositions.

    return 0;
}
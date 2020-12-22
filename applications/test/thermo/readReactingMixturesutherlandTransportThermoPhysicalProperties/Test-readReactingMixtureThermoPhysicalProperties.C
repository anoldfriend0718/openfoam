#include "PtrList.H"
#include "UList.H"
#include "basicSpecieMixture.H"
#include "dictionary.H"
#include "IFstream.H"
#include "Time.H"
#include "argList.H"

#include "sutherlandTransport.H"
#include "janafThermo.H"
#include "fvMesh.H"
#include "pureMixture.H"
#include "rhoReactionThermo.H"
#include "rhoThermo.H"
#include "specie.H"
#include "thermo.H"
#include "thermodynamicConstants.H"
#include "volFieldsFwd.H"
#include "Boussinesq.H"
#include "reactingMixture.H"
#include "perfectGas.H"
#include "sensibleEnthalpy.H"

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

    typedef sutherlandTransport<species::thermo<janafThermo<perfectGas<specie>>,sensibleEnthalpy>> gasThermoType;
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

    



    

    // const auto& thermoTypes=compositions.

    return 0;
}
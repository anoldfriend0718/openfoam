#include "PtrList.H"
#include "UList.H"
#include "basicSpecieMixture.H"
#include "dictionary.H"
#include "IFstream.H"
#include "Time.H"
#include "argList.H"

#include "constTransport.H"

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
#include "hConstThermo.H"

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

    Info<<"Temperature: "<<thermo.T().field()<<endl;
    Info<<"Pressure: "<<thermo.p().field()<<endl;
    Info<<"Density: "<<thermo.rho().field()<<endl;

    PtrList<volScalarField> Y = compositions.Y();
    forAll(Y,i)
    {
        Info<<"mass fraction of species "<<Y[i].name()<<": "<<Y[i].field()<<endl;
    }
    
    PtrList<volScalarField::Internal> c(compositions.species().size());
    forAll(c, i)
    {
        c.set(
            i,
            new volScalarField::Internal
            (
                IOobject
                (
                    "c." + Y[i].name(),
                    runTime.timeName(),
                    mesh,
                    IOobject::NO_READ,
                    IOobject::NO_WRITE
                ),
                mesh,
                dimensionedScalar(dimMoles/dimVolume,0)
            ));
    }

    const volScalarField& rho=thermo.rho();
    forAll(Y,si)
    {
        forAll(rho, celli)
        {
            c[si][celli]=rho[celli]*Y[si][celli]/compositions.Wi(si);
        }
        Info<<"mole concentration of species "<<Y[si].name()<<": "<<c[si].field()<<endl;
    }

    Info<<"Cp: "<<thermo.Cp()().field()<<endl;
    Info<<"ha: "<<thermo.he().field()<<endl;
    Info<<"hc: "<<thermo.hc()().field()<<endl;

    const scalarField h(1,303008);
    const scalarField p(1,1e5);
    const scalarField T0(1,200);
    const labelList cells(1,0);

    Info<<"THE: "<<thermo.THE(h,p,T0,cells)<<endl;

    // PtrList<volScalarField> Y = compositions.Y();
    // forAll(Y,i)
    // {
    //     Info<<Y[i].name()<<" is active: "<<compositions.active(i)
    //         <<" , mole weight: "<<compositions.Wi(i)
    //         <<" , field: "<<Y[i].primitiveField()<<endl;
    // }

    // typedef constTransport<species::thermo<hConstThermo<perfectGas<specie>>,sensibleEnthalpy>> gasThermoType;
    // typedef reactingMixture<gasThermoType> reactingMixtureType;
    // Info<< "check gas thermal type: "<<gasThermoType::typeName()<<endl;
    
    // const reactingMixtureType& mixture=static_cast<const reactingMixtureType&>(compositions);
    // Info<<"speciesCompositionTable: "<<mixture.specieComposition()<<endl;
    // //get the thermo models of different composition in the mixture
    // const PtrList<gasThermoType>& gasThermos=mixture.speciesData();

    // Info<<"composition size: "<<gasThermos.size()<<endl;
    // forAll(gasThermos, i)
    // {
    //     Info<<"compostion ["<<i<<"]: "<<gasThermos[i].name()<<endl;
    // }

    // scalar p0=1e5;
    // scalar T0=300;
    // Info<<"p0: "<<p0<<", T0: "<<T0<<endl;
    // Info<<"Pstd: "<<Pstd<<", Tstd: "<<Tstd<<endl;

    // forAll(gasThermos, i)
    // {
    //     Info<<gasThermos[i].name()
    //         <<", Cp: "<<compositions.Cp(i, p0, T0)
    //         <<", Cv: "<<compositions.Cv(i, p0, T0)
    //         <<", Ha: "<<compositions.Ha(i, p0, T0)
    //         <<", Hc: "<<compositions.Hc(i)
    //         <<", Hs: "<<compositions.Hs(i, p0, T0)
    //         <<", rho: "<<compositions.rho(i, p0,T0)
    //         <<", kappa: "<<compositions.kappa(i, p0, T0)
    //         <<", mu: "<<compositions.mu(i, p0, T0)
    //         <<", alphah: "<<compositions.alphah(i, p0, T0)
    //         <<endl;
    // }

    // Info<<"Mixture: "
    //     <<", Cp: "<<thermo.Cp()().field()
    //     <<", Cv: "<<thermo.Cv()().field()
    //     <<", Cpv: "<<thermo.Cpv()().field()
    //     <<", CpByCpv: "<<thermo.CpByCpv()().field()
    //     <<", He: "<<thermo.he()().field()
    //     <<", rho: "<<thermo.rho()().field()
    //     <<", psi: "<<thermo.psi()().field()
    //     <<", kappa: "<<thermo.kappa()().field()
    //     <<", mu: "<<thermo.mu()().field()
    //     <<", alphahe: "<<thermo.alphahe()().field()
    //     <<endl;

    



    

    // const auto& thermoTypes=compositions.

    return 0;
}
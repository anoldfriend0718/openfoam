#include "dictionary.H"
#include "IFstream.H"
#include "Time.H"
#include "argList.H"
#include "constTransport.H"
#include "eConstThermo.H"
#include "fvMesh.H"
#include "pureMixture.H"
#include "rhoReactionThermo.H"
#include "rhoThermo.H"
#include "specie.H"
#include "volFieldsFwd.H"
#include "Boussinesq.H"

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

    autoPtr<rhoThermo> pThermo(rhoThermo::New(mesh));
    rhoThermo &thermo = pThermo();
    thermo.validate(args.executable(), "h", "e");

    Info << "Base class typeName_: " << thermo.typeName_() << endl;
    Info << "Base class typeName: " << thermo.typeName << endl;
    Info << "type: " << thermo.type() << endl;
    Info << "thermal name: " << thermo.thermoName() << endl;

    const rhoReactionThermo &reactionThermo = static_cast<const rhoReactionThermo &>(thermo);
    const basicMixture &mixture = reactionThermo.composition();

    typedef Foam::pureMixture<constTransport<eConstThermo<Boussinesq<specie>>>> PureMixtureType;
    const PureMixtureType &pureMixture = static_cast<const PureMixtureType &>(mixture);
    Info << "mixture type: " << pureMixture.typeName() << endl;

    Info << "rho field: " << thermo.rho().primitiveField() << endl;
    Info << "psi field: " << thermo.psi().primitiveField() << endl;
    Info << "Cv field: " << thermo.Cv()->primitiveField() << endl;
    Info << "Cp field: " << thermo.Cp()->primitiveField() << endl;
    Info << "HE field: " << thermo.he().primitiveField() << endl;
    Info << "dynamic viscosity mu field: " << thermo.mu()().primitiveField() << endl;
    Info << "thermal conductivity kappa field: " << thermo.kappa()->primitiveField() << endl;
    Info << "thermal diffusivity alphahe field: " << thermo.alphahe()->primitiveField() << endl;

    return 0;
}
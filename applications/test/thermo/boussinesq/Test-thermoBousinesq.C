#include "dictionary.H"

#include "IFstream.H"
#include "constTransport.H"

#include "Boussinesq.H"
#include "eConstThermo.H"
#include "scalar.H"
#include "sensibleInternalEnergy.H"
#include "specie.H"
#include "thermo.H"

using namespace Foam;

int main(int argc, char *argv[]) {
    typedef constTransport<species::thermo<eConstThermo<Boussinesq<specie>>, sensibleInternalEnergy>> BoussinesqType;

    dictionary dict(IFstream("thermalDict")());

    // // Info << dict.subDict("water");
    BoussinesqType waterBoussinesqType(dict.subDict("water"));
    scalar refTemperature(readScalar(dict.subDict("water").subDict("equationOfState").lookup("T0")));
    scalar refDensity(readScalar(dict.subDict("water").subDict("equationOfState").lookup("rho0")));
    scalar beta(readScalar(dict.subDict("water").subDict("equationOfState").lookup("beta")));

    Info << "check water Boussinesq thermo physical properties..." << endl;

    Info << "type name: " << waterBoussinesqType.typeName() << endl;
    Info << "ref Temperature: " << refTemperature << endl;
    Info << "ref Density: " << refDensity << endl;
    Info << "beta: " << beta << endl;

    scalar T1 = 280;
    scalar p1 = 0;

    Info << "RR: " << RR << endl;
    Info << "R: " << waterBoussinesqType.R() << endl;

    Info << "molecular weight[kg/kmol]: " << waterBoussinesqType.W() << endl;
    Info << "psi(" << p1 << "," << T1 << "): " << waterBoussinesqType.psi(p1, T1) << endl;
    Info << "rho(" << p1 << "," << T1 << "): " << waterBoussinesqType.rho(p1, T1) << endl;

    Info << "Cv(" << p1 << "," << T1 << ") [J/kg/k]: " << waterBoussinesqType.Cv(p1, T1) << endl;
    Info << "cv(" << p1 << "," << T1 << ")= Cv*W [J/kmol/k]: " << waterBoussinesqType.cv(p1, T1) << endl;

    Info << "CpMCv(" << p1 << "," << T1 << "): " << waterBoussinesqType.CpMCv(p1, T1) << endl;

    Info << "Cp(" << p1 << "," << T1 << ") [J/kg/k]: " << waterBoussinesqType.Cp(p1, T1) << endl;
    Info << "cp(" << p1 << "," << T1 << ")= Cp*W [J/kmol/k]: " << waterBoussinesqType.cp(p1, T1) << endl;

    Info << "Cpv(" << p1 << "," << T1 << "): " << waterBoussinesqType.Cpv(p1, T1) << endl;
    Info << "CpByCpv(" << p1 << "," << T1 << "): " << waterBoussinesqType.CpByCpv(p1, T1) << endl;

    Info << "Energy type: " << waterBoussinesqType.heName() << endl;
    Info << "Enthalpy/Internal energy [J/kg], HE(" << p1 << "," << T1 << "): " << waterBoussinesqType.HE(p1, T1)
         << endl;
    Info << "Enthalpy/Internal energy [J/kmol], he(" << p1 << "," << T1 << "): " << waterBoussinesqType.he(p1, T1)
         << endl;

    Info << "dynamic viscosity"
         << "mu(" << p1 << "," << T1 << "): " << waterBoussinesqType.mu(p1, T1) << endl;
    Info << "thermal diffusivity of enthalpy"
         << "alphah(" << p1 << "," << T1 << ")=mu/Pr: " << waterBoussinesqType.alphah(p1, T1) << endl;

    Info << "thermal conductivity"
         << "kappa(" << p1 << "," << T1 << ")=rho*mu/Pr: " << waterBoussinesqType.kappa(p1, T1) << endl;

    return 0;
}
#include "IOdictionary.H"
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
#include "UniformField.H"
#include "fvc.H"

using namespace Foam;

class constSolidThermo;
Ostream& operator<<(Ostream&, const constSolidThermo&);

class constSolidThermo
{
    
    public:
        friend Ostream& operator<<(Ostream&, const constSolidThermo&);
        constSolidThermo(const dictionary& dict):
            molWeight(readScalar(dict.lookup("molWeight"))),
            density(readScalar(dict.lookup("density"))),
            kappa(readScalar(dict.lookup("kappa"))),
            Cp(readScalar(dict.lookup("cp")))
        {

        }
        scalar molWeight;
        scalar density;
        scalar kappa;
        scalar Cp;
};

Ostream& operator<<(Ostream& os, const constSolidThermo& thermo)
{
    os<<"W: "<<thermo.molWeight<< token::SPACE;
    os<<"rho: "<<thermo.density<< token::SPACE;
    os<<"kappa: "<<thermo.kappa<< token::SPACE;
    os<<"Cp: "<<thermo.Cp<< token::SPACE;
    return os;
}


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
    
    typedef constTransport<species::thermo<hConstThermo<perfectGas<specie>>,sensibleEnthalpy>> gasThermoType;
    typedef reactingMixture<gasThermoType> reactingMixtureType;
    const basicSpecieMixture& compositions=thermo.composition();
    const reactingMixtureType& mixture=static_cast<const reactingMixtureType&>(compositions);
    //get the thermo models of different composition in the mixture
    const PtrList<gasThermoType>& gasThermos=mixture.speciesData();


    const scalarField& rho=thermo.rho();
    const scalarField& T=thermo.T();
    const scalarField& p=thermo.p();
    forAll(rho,i)
    {
        Info<<"rho: "<<rho<<endl;
        Info<<"T: "<<T<<endl;
        Info<<"p: "<<p<<endl;
    }

    PtrList<volScalarField> Y = mixture.Y();
    label nGasSpecies=mixture.species().size();
    label nTotalSpecies=nGasSpecies+2; //gas species+ coke+rock
    Info<<"nGasSpecies: "<<nGasSpecies<<", nTotalSpecies: "<<nTotalSpecies<<endl;
    label cokeIndex=nTotalSpecies-2;
    label rockIndex=nTotalSpecies-1;

    label O2Index=-1;
    label CO2Index=-1;
    forAll(Y, i)
    {
        Info<<Y[i].name()<<" is active: "<<compositions.active(i)
            <<" , mole weight: "<<compositions.Wi(i)
            <<" , field: "<<Y[i].primitiveField()<<endl;
        if(Y[i].name()=="O2")
        {
            O2Index=i;
        }
        else if(Y[i].name()=="CO2")
        {
            CO2Index=i;
        }
    }
    if(O2Index==-1 || CO2Index==-1)
    {
        FatalErrorInFunction<<"No O2/CO2 in species"
        <<exit(FatalError);
    }
    List<label> reactantIndexs{O2Index,cokeIndex};


    Info<<"Oxygen Index: "<<O2Index<<", CO2 Index: "<<CO2Index<<", cokeIndex: "<<cokeIndex<<endl;
    Info<<"reactantIndex: "<<reactantIndexs<<endl;

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
    Info<<"fluid porosity: "<<coke.field()<<", solid fraction: "<<(1-eps)->field()<<endl;
    Info<<"coke fraction: "<<coke.field()<<", rock fraction: "<<rock.field()<<endl;

    //read the solid thermo
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



    //read the chemistry properties
    IOdictionary chemistryProperties
    (
        IOobject
        (
            thermo.phasePropertyName("chemistryProperties"),
            thermo.db().time().constant(),
            thermo.db(),
            IOobject::MUST_READ_IF_MODIFIED,
            IOobject::NO_WRITE
        )
    );

    const scalar deltaTChemIni(readScalar(chemistryProperties.lookup("initialChemicalTimeStep")));
    const scalar deltaTChemMax(chemistryProperties.lookupOrDefault("maxChemicalTimeStep", great));
    volScalarField::Internal deltaTChem
    (
        IOobject
        (
            thermo.phasePropertyName("deltaTChem"),
            runTime.constant(),
            mesh,
            IOobject::NO_READ,
            IOobject::NO_WRITE
        ),
        mesh,
        dimensionedScalar(dimTime, deltaTChemIni)

    );
    Info<<"deltaTChemIni: "<<deltaTChemIni<<", deltaTChemMax: "<<deltaTChemMax<<endl;

    //read reactions 
    IOdictionary chemicalDict
    (
        IOobject
        (
            "reactions",
            thermo.db().time().constant(),
            thermo.db(),
            IOobject::MUST_READ_IF_MODIFIED,
            IOobject::NO_WRITE
        )
    );
    const dictionary cokeCombustionDict=chemicalDict.subDict("customizedReactions").subDict("cokeCombustion");
    const scalar A(readScalar(cokeCombustionDict.lookup("A")));
    const scalar Ta(readScalar(cokeCombustionDict.lookup("Ta"))); //E/R
    const scalar hr(readScalar(cokeCombustionDict.lookup("hr")));
    Info<<"coke combustion kinetics: A: "<<A<<", Ta: "<<Ta<<", hr: "<<hr<<endl;

    volScalarField RRO2
    (
        IOobject
        (
            "R_O2",
            runTime.timeName(),
            mesh,
            IOobject::NO_READ,
            IOobject::AUTO_WRITE
        ),
        mesh,
        dimensionedScalar("Ri", dimensionSet(1, -3, -1, 0, 0, 0, 0), 0)
    );

    volScalarField RRCO2
    (
        IOobject
        (
            "R_CO2",
            runTime.timeName(),
            mesh,
            IOobject::NO_READ,
            IOobject::AUTO_WRITE
        ),
        mesh,
        dimensionedScalar("Ri", dimensionSet(1, -3, -1, 0, 0, 0, 0), 0)
    );

    volScalarField RRCoke
    (
        IOobject
        (
            "R_coke",
            runTime.timeName(),
            mesh,
            IOobject::NO_READ,
            IOobject::AUTO_WRITE
        ),
        mesh,
        dimensionedScalar("Ri", dimensionSet(1, -3, -1, 0, 0, 0, 0), 0)
    );



    scalarField c(nTotalSpecies,Zero);
    scalarField c0(nTotalSpecies,Zero);
    UniformField<scalar> deltaT(runTime.deltaTValue());
    Info<<"Flow Delta Time: "<<deltaT.field()<<endl;
    
    forAll(rho,celli)
    {
        const scalar rhoi = rho[celli];
        scalar Ti=T[celli]; // not reference here, donot update the temperature in T[celli], but update Ti
        const scalar pi=p[celli];
        const scalar epsi = eps[celli];
        for(label i=0;i<nGasSpecies;i++)
        {
            c[i]=epsi*rhoi*Y[i][celli]/gasThermos[i].W();
        }

        const scalar cokeFraci=coke[celli];
        const scalar rockFraci=rock[celli];
        c[cokeIndex]=cokeFraci*cokeThermo.density/cokeThermo.molWeight;
        c[rockIndex]=rockFraci*rockThermo.density/rockThermo.molWeight;

        for(label i=0;i<nTotalSpecies;i++)
        {
            c0[i]=c[i];
        }

        Info<<"c0: "<<c0<<endl;

        // Initialise time progress
        scalar timeLeft = deltaT[celli]; //flow time step
        scalar& subDeltaT=deltaTChem[celli]; //Latest estimation of integration step, it will be updated in the finer chemical step
        volScalarField tempCoke(coke);
        // Calculate the chemical source terms
        while (timeLeft > small)
        {
            scalar dt = timeLeft;

            scalar Cpf=c[0]*mixture.Wi(0)*mixture.Cp(0, pi, Ti);
            for(label i=1;i<nGasSpecies;i++)
            {
                Cpf+=c[i]*mixture.Wi(i)*mixture.Cp(i, pi, Ti);
            }
            scalar Cps=c[nTotalSpecies-2]*cokeThermo.molWeight*cokeThermo.Cp+
                       c[nTotalSpecies-1]*rockThermo.molWeight*rockThermo.Cp;
            scalar ha=(Cpf+Cps)*Ti;

            Info<<"timeLeft: "<<dt<<endl;
            Info<<">>>>>>>>>>>>>>>>>>>>>>>"<<endl;
            Info<<"Cpf: "<<Cpf<<", Cps: "<<Cps<<", Ti: "<<Ti<<", Ha: "<<ha<<endl;


            volScalarField ssArea
            (
                "ssArea",
                mag(fvc::grad(tempCoke))*(4.0*tempCoke*(1-tempCoke))
            );
            Info<<"mag(fvc::grad(tempCoke)): "<<mag(fvc::grad(tempCoke))->field()<<endl;
            Info<<"effective surface area: "<<ssArea.field()<<endl;

            const scalar ak=ssArea[celli]*A*std::exp(-Ta/Ti);
            const scalar cokeReactionRate=ak*c[O2Index];
            Info<<"coke reaction rate: "<<cokeReactionRate<<endl;
           
            // Calculate the stable/accurate time-step
            scalar tMin = great;
            forAll(reactantIndexs,i)
            {
                label si=reactantIndexs[i];
                if(cokeReactionRate>small) //if cokeReactionRate greater than zero
                {
                    tMin=min(tMin,c[si]/cokeReactionRate);
                }
            }
            subDeltaT=tMin; //Latest estimation of integration step was updated
            //compare the left flow time and the latest estimation of integration step to get the chemical time step
            dt=min(dt,subDeltaT); 

            Info<<"chemical time step: "<<dt<<", latest estimation of integration step: "<<subDeltaT<<endl;

            //solve for the new composition
            //Euler Implict
            // scalar c_O2_new=(c[O2Index]/dt)/(ak+1/dt);

            //2nd R-K method
            // scalar f1=-ak*c[O2Index];
            // scalar f2=-ak*(c[O2Index]+dt*f1);
            // scalar c_O2_new=c[O2Index]+0.5*dt*(f1+f2);

            //4th R-K method  
            scalar f1=-ak*c[O2Index];
            scalar f2=-ak*(c[O2Index]+dt/2.0*f1);
            scalar f3=-ak*(c[O2Index]+dt/2.0*f2);
            scalar f4=-ak*(c[O2Index]+dt*f3);
            scalar c_O2_new=c[O2Index]+dt/6.0*(f1+2.0*f2+2.0*f3+f4);


            scalar deltaC_O2=c_O2_new-c[O2Index];
            Info<<"c_O2_new: "<<c_O2_new<<", deltaC_O2: "<<deltaC_O2<<endl;
            c[O2Index]+=deltaC_O2;
            c[CO2Index]-=deltaC_O2;
            c[cokeIndex]+=deltaC_O2;
            // Limit the composition
            for (label i=0; i<nTotalSpecies; i++)
            {
                c[i] = max(0, c[i]);
            }
            Info<<"new c updated: "<<c<<endl;

            //solve for new coke fraction
            tempCoke[celli]=c[cokeIndex]*cokeThermo.molWeight/cokeThermo.density;
            Info<<"new coke fraction updated: "<<tempCoke[celli]<<endl;
           
            // solve for the new temperature
            Cpf=c[0]*mixture.Wi(0)*mixture.Cp(0, pi, Ti); //Species Cp is not related to p and T, but mixture depend on c
            for(label i=1;i<nGasSpecies;i++)
            {
                Cpf+=c[i]*mixture.Wi(i)*mixture.Cp(i, pi, Ti);
            }
            Cps=c[nTotalSpecies-2]*cokeThermo.molWeight*cokeThermo.Cp+
                       c[nTotalSpecies-1]*rockThermo.molWeight*rockThermo.Cp;

            Ti=(ha-deltaC_O2*hr)/(Cpf+Cps);
            Info<<"new T updated: "<<Ti<<", with cumulative increment T: "<<Ti-T[celli]<<endl;

            timeLeft -= dt;
            Info<<"======================"<<endl;
        }

        deltaTChem[celli] =min(deltaTChem[celli], deltaTChemMax);

        RRO2[celli]=(c[O2Index]-c0[O2Index])*mixture.Wi(O2Index)/deltaT[celli];
        RRCO2[celli]=(c[CO2Index]-c0[CO2Index])*mixture.Wi(CO2Index)/deltaT[celli];
        RRCoke[celli]=(c[cokeIndex]-c0[cokeIndex])*cokeThermo.molWeight/deltaT[celli];
        

    }

    Info<<"Time-splitting RR O2: "<<RRO2.field()<<endl;
    Info<<"Time-splitting RR CO2: "<<RRCO2.field()<<endl;
    Info<<"Time-splitting RR coke: "<<RRCoke.field()<<endl;


    volScalarField ssArea
    (
        "ssArea",
        mag(fvc::grad(coke))*(4.0*coke*(1-coke))
    );
    const scalar ak=ssArea[0]*A*std::exp(-Ta/T[0]);
    const scalar cokeReactionRate=ak*c0[O2Index];

    Info<<"Normal RR O2: "<<-cokeReactionRate*mixture.Wi(O2Index)<<endl;
    
    




  

    return 0;
}
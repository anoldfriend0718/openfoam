#include "IOdictionary.H"
#include "PtrList.H"
#include "UList.H"
#include "basicSpecieMixture.H"
#include "dictionary.H"
#include "IFstream.H"
#include "Time.H"
#include "argList.H"

#include "constTransport.H"

#include "fvMatricesFwd.H"
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

class cokeCombustion:
    public IOdictionary
{
    private:
        const Switch chemistry_;
        //- Initial chemical time step
        const scalar deltaTChemIni_;

        //- Maximum chemical time step
        const scalar deltaTChemMax_;

        //- Latest estimation of integration step
        volScalarField::Internal deltaTChem_;

        // volScalarField RRO2;
        const fvMesh& mesh_;
        const volScalarField& eps_;
        const volScalarField& coke_;
        const volScalarField& rock_;


        const rhoReactionThermo& thermo_;
        const basicSpecieMixture& compositions_;
        const constSolidThermo& cokeThermo_;
        const constSolidThermo& rockThermo_;
        const PtrList<volScalarField>& Y_;
        const label nGasSpecies_;
        const label nTotalSpecies_; 

         //- Temporary concentration field
        mutable scalarField c_;

         //- Temporary initial concentration field
        mutable scalarField c0_;


        label O2Index_;
        label CO2Index_;
        label cokeIndex_;
        label rockIndex_;
        List<label> reactantIndexs_;

        const IOdictionary chemicalDict_;
        const dictionary cokeCombustionDict_;
        const scalar A_;
        const scalar Ta_;
        const scalar hr_;
        const word odeSolver_;

        volScalarField::Internal RRO2_;
        volScalarField::Internal RRCO2_;
        volScalarField::Internal RRCoke_;

        void init()
        {
            sortIndex();
        }

        void sortIndex()
        {
            cokeIndex_=nTotalSpecies_-2;
            rockIndex_=nTotalSpecies_-1;

            O2Index_=-1;
            CO2Index_=-1;
            forAll(Y_, i)
            {
                // Info<<Y[i].name()<<" is active: "<<compositions.active(i)
                //     <<" , mole weight: "<<compositions.Wi(i)
                //     <<" , field: "<<Y[i].primitiveField()<<endl;
                if(Y_[i].name()=="O2")
                {
                    O2Index_=i;
                }
                else if(Y_[i].name()=="CO2")
                {
                    CO2Index_=i;
                }
            }
            if(O2Index_==-1 || CO2Index_==-1)
            {
                FatalErrorInFunction<<"No O2/CO2 in species"
                <<exit(FatalError);
            }
            reactantIndexs_.append(O2Index_);
            reactantIndexs_.append(cokeIndex_);
        }

        void solvei(scalarField& c,scalar& Ti,scalar& cokei,
                    const scalar pi, const scalar ssi,
                    scalar& dt,scalar& subDeltaT)
        {
            scalar Cpf=c[0]*compositions_.Wi(0)*compositions_.Cp(0, pi, Ti);
            for(label i=1;i<nGasSpecies_;i++)
            {
                Cpf+=c[i]*compositions_.Wi(i)*compositions_.Cp(i, pi, Ti);
            }
            scalar Cps=c[cokeIndex_]*cokeThermo_.molWeight*cokeThermo_.Cp+
                        c[rockIndex_]*rockThermo_.molWeight*rockThermo_.Cp;
            scalar ha=(Cpf+Cps)*Ti;

            Info<<"timeLeft: "<<dt<<endl;
            Info<<">>>>>>>>>>>>>>>>>>>>>>>"<<endl;
            Info<<"Cpf: "<<Cpf<<", Cps: "<<Cps<<", Ti: "<<Ti<<", Ha: "<<ha<<endl;

            const scalar filterDepth=4.0*cokei*(1-cokei);
            const scalar ak=ssi*filterDepth*A_*std::exp(-Ta_/Ti);
            const scalar cokeReactionRate=ak*c[O2Index_];
            Info<<"coke reaction rate: "<<cokeReactionRate<<endl;
        
            // Calculate the stable/accurate time-step
            scalar tMin = great;
            forAll(reactantIndexs_,i)
            {
                label si=reactantIndexs_[i];
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
            //4th R-K method 
            scalar c_O2_new=small;
            if(odeSolver_=="EulerImplicit") 
            {
                c_O2_new=solveODEByEulerImplicit(ak,c[O2Index_],dt);
            }
            else if(odeSolver_=="2ndRK")
            {
                c_O2_new=solveODEBy2ndRK(ak,c[O2Index_],dt);
            }
            else
            {
                c_O2_new=solveODEBy4thRK(ak,c[O2Index_],dt);
            }
            
            scalar deltaC_O2=c_O2_new-c[O2Index_];
            Info<<"c_O2_new: "<<c_O2_new<<", deltaC_O2: "<<deltaC_O2<<endl;
            c[O2Index_]+=deltaC_O2;
            c[CO2Index_]-=deltaC_O2;
            c[cokeIndex_]+=deltaC_O2;
            // Limit the composition
            for (label i=0; i<nTotalSpecies_; i++)
            {
                c[i] = max(0, c[i]);
            }
            Info<<"new c updated: "<<c<<endl;

            //solve for new coke fraction
            cokei=c[cokeIndex_]*cokeThermo_.molWeight/cokeThermo_.density;
            Info<<"new coke fraction updated: "<<cokei<<endl;
        
            // solve for the new temperature
            Cpf=c[0]*compositions_.Wi(0)*compositions_.Cp(0, pi, Ti); //Species Cp is not related to p and T, but mixture depend on c
            for(label i=1;i<nGasSpecies_;i++)
            {
                Cpf+=c[i]*compositions_.Wi(i)*compositions_.Cp(i, pi, Ti);
            }
            Cps=c[cokeIndex_]*cokeThermo_.molWeight*cokeThermo_.Cp+
                c[rockIndex_]*rockThermo_.molWeight*rockThermo_.Cp;

            Ti=(ha-deltaC_O2*hr_)/(Cpf+Cps);
            Info<<"new T updated: "<<Ti<<endl;
        }
 

        inline scalar solveODEByEulerImplicit(const scalar ak,const scalar ci,const scalar dt)
        {
            scalar c_O2_new=(ci/dt)/(ak+1/dt);
            return c_O2_new;
        }

        inline scalar solveODEBy2ndRK(const scalar ak,const scalar ci,const scalar dt)
        {
            scalar f1=-ak*ci;
            scalar f2=-ak*(ci+dt*f1);
            scalar c_O2_new=ci+0.5*dt*(f1+f2);
            return c_O2_new;
        }

        inline scalar solveODEBy4thRK(const scalar ak,const scalar ci,const scalar dt)
        {
            scalar f1=-ak*ci;
            scalar f2=-ak*(ci+dt/2.0*f1);
            scalar f3=-ak*(ci+dt/2.0*f2);
            scalar f4=-ak*(ci+dt*f3);
            scalar c_O2_new=ci+dt/6.0*(f1+2.0*f2+2.0*f3+f4);
            return c_O2_new;
        }

    public:
        cokeCombustion(const fvMesh& mesh, const rhoReactionThermo& thermo, const constSolidThermo& cokeThermo,
                       const constSolidThermo& rockThermo):
            IOdictionary
            (
                IOobject
                (
                    thermo.phasePropertyName("chemistryProperties"),
                    thermo.db().time().constant(),
                    thermo.db(),
                    IOobject::MUST_READ_IF_MODIFIED,
                    IOobject::NO_WRITE
                )
            ),
            chemistry_(lookup("chemistry")),
            deltaTChemIni_(readScalar(lookup("initialChemicalTimeStep"))),
            deltaTChemMax_(lookupOrDefault("maxChemicalTimeStep", great)),
            deltaTChem_
            (
                IOobject
                (
                    thermo.phasePropertyName("deltaTChem"),
                    mesh.time().constant(),
                    mesh,
                    IOobject::NO_READ,
                    IOobject::NO_WRITE
                ),
                mesh,
                dimensionedScalar(dimTime, deltaTChemIni_)
            ),
            mesh_(mesh),
            eps_(mesh.lookupObject<volScalarField>("eps")),
            coke_(mesh.lookupObject<volScalarField>("coke")),
            rock_(mesh.lookupObject<volScalarField>("rock")),
            thermo_(thermo),
            compositions_(thermo.composition()),
            cokeThermo_(cokeThermo),
            rockThermo_(rockThermo),
            Y_(thermo.composition().Y()),
            nGasSpecies_(Y_.size()),
            nTotalSpecies_(nGasSpecies_+2),
            c_(nTotalSpecies_,Zero),
            c0_(nTotalSpecies_,Zero),
            reactantIndexs_({}),
            chemicalDict_
            (
                IOobject
                (
                    "reactions",
                    thermo.db().time().constant(),
                    thermo.db(),
                    IOobject::MUST_READ_IF_MODIFIED,
                    IOobject::NO_WRITE
                )
            ),
            cokeCombustionDict_
            (   
                chemicalDict_
                    .subDict("customizedReactions")
                    .subDict("cokeCombustion")
            ),
            A_(readScalar(cokeCombustionDict_.lookup("A"))),
            Ta_(readScalar(cokeCombustionDict_.lookup("Ta"))),
            hr_(readScalar(cokeCombustionDict_.lookup("hr"))),
            odeSolver_(cokeCombustionDict_.lookup("odeSolver")),
            RRO2_
            (
                IOobject
                (
                    "R_O2",
                    mesh.time().timeName(),
                    mesh,
                    IOobject::NO_READ,
                    IOobject::AUTO_WRITE
                ),
                mesh,
                dimensionedScalar("Ri", dimensionSet(1, -3, -1, 0, 0, 0, 0), 0)
            ),
            RRCO2_
            (
                IOobject
                (
                    "R_CO2",
                    mesh.time().timeName(),
                    mesh,
                    IOobject::NO_READ,
                    IOobject::AUTO_WRITE
                ),
                mesh,
                dimensionedScalar("Ri", dimensionSet(1, -3, -1, 0, 0, 0, 0), 0)
            ),
            RRCoke_
            (
                IOobject
                (
                    "R_coke",
                    mesh.time().timeName(),
                    mesh,
                    IOobject::NO_READ,
                    IOobject::AUTO_WRITE
                ),
                mesh,
                dimensionedScalar("Ri", dimensionSet(1, -3, -1, 0, 0, 0, 0), 0)
            )
        {
            init();
            Info<<"chemistry on/off: "<<chemistry_<<endl;
            Info<<"coke combustion ode solver: "<<odeSolver_<<endl;
        }

        void solve(const scalar deltaTValue)
        {
            if (!this->chemistry_)
            {
                return;
            }

            UniformField<scalar> deltaT(deltaTValue);
            // Info<<"Flow Delta Time: "<<deltaTValue<<endl;

            const scalarField& rho=thermo_.rho();
            const scalarField& T=thermo_.T();
            const scalarField& p=thermo_.p();

            forAll(rho,i)
            {
                Info<<"rho: "<<rho<<endl;
                Info<<"T: "<<T<<endl;
                Info<<"p: "<<p<<endl;
            }

            volScalarField cokeSpeciesSurfaceArea //without filter depth
            (
                "cokeSpeciesSurfaceArea",
                mag(fvc::grad(coke_))
            );

            // Info<<"coke surface area without counting filter depth: "<<cokeSpeciesSurfaceArea.field()<<endl;

            forAll(rho,celli)
            {
                 // not reference here, donot update the temperature in T[celli], but update Ti
                scalar Ti=T[celli];
                //  donot update the coke fraction in coke[celli], but update cokei
                scalar cokei=coke_[celli];
                c_[cokeIndex_]=cokei*cokeThermo_.density/cokeThermo_.molWeight;
                const scalar rocki=rock_[celli];
                c_[rockIndex_]=rocki*rockThermo_.density/rockThermo_.molWeight;

                const scalar rhoi = rho[celli];  
                const scalar pi=p[celli];
                const scalar ssi=cokeSpeciesSurfaceArea[celli];
                const scalar epsi = eps_[celli];
                for(label i=0;i<nGasSpecies_;i++)
                {
                    c_[i]=epsi*rhoi*Y_[i][celli]/compositions_.Wi(i);
                }

                for(label i=0;i<nTotalSpecies_;i++)
                {
                    c0_[i]=c_[i];
                }

                Info<<"c0: "<<c0_<<endl;

                // Initialise time progress
                scalar timeLeft = deltaT[celli]; //flow time step
                scalar& subDeltaT=deltaTChem_[celli]; //Latest estimation of integration step, it will be updated in the finer chemical step
   
            // Calculate the chemical source terms
                while (timeLeft > small)
                {
                    scalar dt = timeLeft;
                    solvei(c_, Ti, cokei, pi, ssi, dt, subDeltaT);
                    timeLeft -= dt;
                    Info<<"======================"<<endl;   
                }

                deltaTChem_[celli] =min(deltaTChem_[celli], deltaTChemMax_);

                RRO2_[celli]=(c_[O2Index_]-c0_[O2Index_])*compositions_.Wi(O2Index_)/deltaT[celli];
                RRCO2_[celli]=(c_[CO2Index_]-c0_[CO2Index_])*compositions_.Wi(CO2Index_)/deltaT[celli];
                RRCoke_[celli]=(c_[cokeIndex_]-c0_[cokeIndex_])*cokeThermo_.molWeight/deltaT[celli];
            }

  
        }

        inline const volScalarField::Internal& RRO2() const
        {
            return RRO2_;
        };

        inline const volScalarField::Internal& RRCO2() const
        {
            return RRCO2_;
        };

        inline const volScalarField::Internal& RRCoke() const 
        {
            return RRCoke_;
        };

        tmp<fvScalarMatrix> R(volScalarField& Y) const
        {
            tmp<fvScalarMatrix> tSu(new fvScalarMatrix(Y, dimMass/dimTime));
            fvScalarMatrix& Su = tSu.ref();
            const word speciesName=Y.name();
            if(speciesName=="O2")
            {
                Su+=RRO2_;
            }
            else if(speciesName=="CO2")
            {
                Su+=RRCO2_;
            }
            else
            {
                FatalErrorInFunction<<"No valid species for cokeCombustion R method"
                    <<exit(FatalError);
            }
            return tSu;
        }

        tmp<fvScalarMatrix> Rs(volScalarField& frac) const
        {
            tmp<fvScalarMatrix> tSu(new fvScalarMatrix(frac, dimMass/dimTime));
            fvScalarMatrix& Su = tSu.ref();
            const word solidSpeciesName=frac.name();
            if(solidSpeciesName=="coke")
            {
                Su+=RRCoke_;
            }
            else
            {
                FatalErrorInFunction<<"No valid solid species for cokeCombustion Rs method"
                    <<exit(FatalError);
            }
            return tSu;
        }

        tmp<volScalarField> Qdot() const
        {
            tmp<volScalarField> tQdot
            (
                volScalarField::New
                (
                    "Qdot",
                    this->mesh_,
                    dimensionedScalar(dimEnergy/dimVolume/dimTime, 0)
                )
            );

            if (!this->chemistry_)
            {
                return tQdot;
            }
            scalarField& Qdot = tQdot.ref();

            forAll(Qdot, celli)
            {
                //RR is calculated by the standardChemistryModel solve method (integrate reaction rate is on) or calculate method (integrated reaction rate is off)
                Qdot[celli] -= hr_*RRCoke_[celli]/cokeThermo_.molWeight;
            }
            return tQdot;
        }

        tmp<volScalarField::Internal> calculateTransientRRO2()
        {
            tmp<volScalarField::Internal> tRRO2
            (
                volScalarField::Internal::New
                (
                    "tR_CO2",
                    mesh_,
                    dimensionedScalar("Ri", dimensionSet(1, -3, -1, 0, 0, 0, 0), 0)
                )
            );
            
            if (!this->chemistry_)
            {
                return tRRO2;
            }

            volScalarField::Internal& RRO2 = tRRO2.ref();
            volScalarField ssArea
            (
                "ssArea",
                mag(fvc::grad(coke_))*(4.0*coke_*(1-coke_))
            );
            const scalarField& T=thermo_.T();
            const scalarField& rho=thermo_.rho();
     
            scalar aki=0;
            forAll(T, celli)
            {
                const scalar rhoi = rho[celli];  
                const scalar epsi = eps_[celli];
                for(label i=0;i<nGasSpecies_;i++)
                {
                    c0_[i]=epsi*rhoi*Y_[i][celli]/compositions_.Wi(i);
                }
                aki=ssArea[celli]*A_*std::exp(-Ta_/T[celli]);
                RRO2[celli]=-aki*c0_[O2Index_]*compositions_.Wi(O2Index_);
            }
            return tRRO2;
        }

};

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
    Info<<"fluid porosity: "<<coke.field()<<", solid fraction: "<<(1-eps)->field()<<endl;
    Info<<"coke fraction: "<<coke.field()<<", rock fraction: "<<rock.field()<<endl;

    
 

    cokeCombustion combustion(mesh,thermo,cokeThermo,rockThermo);
    scalar flowDt=runTime.deltaTValue();
    combustion.solve(flowDt);

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

    Info<<"Test Qdot..."<<endl;
    Info<<"Qdot: "<<combustion.Qdot()->field()<<endl;


    return 0;
}
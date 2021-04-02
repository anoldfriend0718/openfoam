/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     | Website:  https://openfoam.org
    \\  /    A nd           | Copyright (C) 2011-2019 OpenFOAM Foundation
     \\/     M anipulation  |
-------------------------------------------------------------------------------
License
    This file is part of OpenFOAM.

    OpenFOAM is free software: you can redistribute it and/or modify it
    under the terms of the GNU General Public License as published by
    the Free Software Foundation, either version 3 of the License, or
    (at your option) any later version.

    OpenFOAM is distributed in the hope that it will be useful, but WITHOUT
    ANY WARRANTY; without even the implied warranty of MERCHANTABILITY or
    FITNESS FOR A PARTICULAR PURPOSE.  See the GNU General Public License
    for more details.

    You should have received a copy of the GNU General Public License
    along with OpenFOAM.  If not, see <http://www.gnu.org/licenses/>.

\*---------------------------------------------------------------------------*/

#include "cokeCombustion.H"
#include "constSolidThermo.H"
#include "UniformField.H"
#include "dictionary.H"
#include "fvc.H"
#include "addToRunTimeSelectionTable.H"

namespace Foam
{
    defineTypeNameAndDebug(cokeCombustion,1);
}

Foam::cokeCombustion::cokeCombustion(const Foam::fvMesh& mesh, 
                                    const Foam::rhoReactionThermo& thermo, 
                                    const Foam::constSolidThermo& cokeThermo,
                                    const Foam::constSolidThermo& rockThermo):
    IOdictionary
    (
        IOobject
        (
            "chemistryProperties",
            thermo.db().time().constant(),
            thermo.db(),
            IOobject::MUST_READ_IF_MODIFIED,
            IOobject::NO_WRITE
        )
    ),
    chemistry_(lookup("chemistry")),
    integrateReactionRate_(lookupOrDefault("integrateReactionRate", true)),
    chemicalTimeCoeff_(lookupOrDefault("chemicalTimeCoeff", 1.0)),
    deltaTChemIni_(readScalar(lookup("initialChemicalTimeStep"))),
    deltaTChemMax_(lookupOrDefault("maxChemicalTimeStep", great)),
    deltaTChem_
    (
        IOobject
        (
            "deltaTChem",
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
    odeSolver_(cokeCombustionDict_.lookupOrDefault("odeSolver",word("4thRKFull"))),
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
    if (integrateReactionRate_)
    {
        Info<< "using integrated reaction rate" << endl;
        Info<<"coke combustion ode solver: "<<odeSolver_<<endl;
        Info<<"chemical time coefficient: "<<chemicalTimeCoeff_<<endl;
    }
    else
    {
        Info<< "using instantaneous reaction rate" << endl;
    } 
    
}

void Foam::cokeCombustion::init()
{
    sortIndex();
}

void Foam::cokeCombustion::sortIndex()
{
    cokeIndex_=nTotalSpecies_-2;
    rockIndex_=nTotalSpecies_-1;

    O2Index_=-1;
    CO2Index_=-1;
    forAll(Y_, i)
    {
        if(debug)
        {
            Info<<Y_[i].name()<<" index: "<<i<<endl;
        }
        
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

void Foam::cokeCombustion::correct()
{
    const Time& runtime=mesh_.time();
    scalar flowDeltaT=runtime.deltaTValue();

    if(debug)
    {
        Info<<endl;
        Info<<"@@@@@@@@@@@@@@@@@ TimeIndex="<<runtime.timeIndex()<<" @@@@@@@@@@@@@@@@@"<<endl;
        Info<<"about to correct coke combustion at Time "
            <<runtime.timeName()
            <<", with time step: "
            <<flowDeltaT
            <<endl;
    }
    
    if(integrateReactionRate_)
    {
        this->solve(flowDeltaT);
    }
    else
    {
        this->calculate();
    }
    
    if(debug)
    {
        Info<<"complete correct coke combustion at Time "
            <<runtime.timeName()<<endl;
        Info<<"$$$$$$$$$$$$$$$$$ TimeIndex="<<runtime.timeIndex()<<" $$$$$$$$$$$$$$$$$"<<endl;
        Info<<endl;
    }
}

void Foam::cokeCombustion::solve(const scalar deltaTValue)
{
    if (!this->chemistry_)
    {
        return;
    }

    UniformField<scalar> deltaT(deltaTValue);
    // Info<<"Flow Delta Time: "<<deltaTValue<<endl;

    const scalarField& rho=thermo_.rho()();
    const scalarField& T=thermo_.T();
    const scalarField& p=thermo_.p();
    volScalarField cokeGrad //without filter depth
    (
        "cokeGrad",
        mag(fvc::grad(coke_))
    );

    if(debug>1)
    {
        Info<<"rho: "<<rho<<endl;
        Info<<"T: "<<T<<endl;
        Info<<"p: "<<p<<endl;
        Info<<"Y O2: "<<Y_[O2Index_].field()<<endl;
        Info<<"eps: "<<eps_.field()<<endl;
        Info<<"coke: "<<coke_.field()<<endl;
        Info<<"coke gradient: "
            <<cokeGrad.field()<<endl;
        // Info<<"filter depth: "
        //     <<(4*coke_*(1-coke_))->field()<<endl;
    }
    
    forAll(rho,celli)
    {
            // not reference here, donot update the temperature in T[celli], but update Ti
        scalar Ti=T[celli];
        //  donot update the coke fraction in coke[celli], but update cokei
        scalar cokei=coke_[celli];
        c_[cokeIndex_]=(cokei*cokeThermo_.density/cokeThermo_.molWeight).value();
        const scalar rocki=rock_[celli];
        c_[rockIndex_]=(rocki*rockThermo_.density/rockThermo_.molWeight).value();

        const scalar rhoi = rho[celli];  
        const scalar pi=p[celli];
        const scalar cokeGradi=cokeGrad[celli];
        const scalar epsi = eps_[celli];
        for(label i=0;i<nGasSpecies_;i++)
        {
            c_[i]=epsi*rhoi*Y_[i][celli]/compositions_.Wi(i);
        }

        for(label i=0;i<nTotalSpecies_;i++)
        {
            c0_[i]=c_[i];
        }

        if(debug>1)
        {
            Info<<"cell ["<<celli<<"]: "
                <<"c0: "<<c0_<<endl;
        }
        
        // Initialise time progress
        scalar timeLeft = deltaT[celli]; //flow time step
        scalar& subDeltaT=deltaTChem_[celli]; //Latest estimation of integration step, it will be updated in the finer chemical step

    // Calculate the chemical source terms
        label nChemicalTimeStep=0;
        while (timeLeft > small) //2.22045e-16
        {
            if(debug>1)
            {
                Info<<"cell ["<<celli<<"]: "
                    <<"timeLeft: "<<timeLeft<<endl;
                Info<<"cell ["<<celli<<"]: "
                    <<">>>>>>>>>>>>>>>>>>>>>>> Chemical Time Index= "
                    <<nChemicalTimeStep
                    <<" >>>>>>>>>>>>>>>>>>>>>>>"
                    <<endl;   
            }

            scalar dt = timeLeft;
            solvei(c_, Ti, cokei, pi, cokeGradi, dt, subDeltaT);
            timeLeft -= dt;
            
            if(debug>1)
            {
                Info<<"cell ["<<celli<<"]: "
                    <<"======================  Chemical Time Index= "
                    <<nChemicalTimeStep
                    <<" ======================"
                    <<endl
                    <<endl;
            }
            nChemicalTimeStep +=1;
             
        }

        deltaTChem_[celli] =min(deltaTChem_[celli], deltaTChemMax_);

        RRO2_[celli]=(c_[O2Index_]-c0_[O2Index_])*compositions_.Wi(O2Index_)/deltaT[celli];
        RRCO2_[celli]=(c_[CO2Index_]-c0_[CO2Index_])*compositions_.Wi(CO2Index_)/deltaT[celli];
        RRCoke_[celli]=(c_[cokeIndex_]-c0_[cokeIndex_])*cokeThermo_.molWeight.value()/deltaT[celli];
    }
}


void Foam::cokeCombustion::solvei(scalarField& c,scalar& Ti,scalar& cokei,
                                const scalar pi, const scalar cokeGradi,
                                scalar& dt,scalar& subDeltaT)
{
    scalar Cpf=c[0]*compositions_.Wi(0)*compositions_.Cp(0, pi, Ti);
    for(label i=1;i<nGasSpecies_;i++)
    {
        Cpf+=c[i]*compositions_.Wi(i)*compositions_.Cp(i, pi, Ti);
    }
    scalar Cps=c[cokeIndex_]*cokeThermo_.molWeight.value()*cokeThermo_.Cp.value()+
                c[rockIndex_]*rockThermo_.molWeight.value()*rockThermo_.Cp.value();
    scalar ha=(Cpf+Cps)*Ti;

    if(debug>1)
    {
        Info<<"Cpf: "<<Cpf<<", Cps: "<<Cps<<", Ti: "<<Ti<<", Ha: "<<ha<<endl;
    }

    // const scalar filterDepth=4.0*cokei*(1-cokei);
    const scalar filterDepth=(1-(1-cokei)*(1-cokei)); //either coke or eps in one cell
    const scalar effssi=2.0*cokeGradi*filterDepth; // need a factor of 2 to correct the species surface area
    const scalar ak=effssi*A_*std::exp(-Ta_/Ti);
    const scalar epsi=(1-cokei);
    const scalar cokeReactionRate=ak*(c[O2Index_]/epsi); //use the bulk oxygen concentration
    if(debug>1)
    {
        Info<<"coke reaction rate: "<<cokeReactionRate<<endl;
    }

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
    subDeltaT=tMin;  //Latest estimation of integration step was updated
    //compare the left flow time and the latest estimation of integration step to get the chemical time step
    dt=min(dt,subDeltaT); 

    if(debug>1)
    {
        Info<<"chemical time step: "<<dt
            <<", latest estimation of integration step: "<<subDeltaT<<endl;
    }


    //solve for the new composition
    //4th R-K method 
    scalar c_O2_new=small;
    if(odeSolver_=="EulerImplicit") 
    {
        c_O2_new=solveODEByEulerImplicit(ak,c[O2Index_],dt,epsi);
    }
    else if(odeSolver_=="2ndRK")
    {
        c_O2_new=solveODEBy2ndRK(ak,c[O2Index_],dt,epsi);
    }
    else if(odeSolver_=="4thRK")
    {
        c_O2_new=solveODEBy4thRK(ak,c[O2Index_],dt,epsi);
    }
    else 
    {
        c_O2_new=solveODEBy4thRKFull(effssi,Ti,c[O2Index_],Cps,Cpf,dt,epsi);
    }
    
    scalar deltaC_O2=c_O2_new-c[O2Index_];
    if(debug>1)
    {
        Info<<"c_O2_new: "<<c_O2_new<<", deltaC_O2: "<<deltaC_O2<<endl;
    }
    
    c[O2Index_]+=deltaC_O2;
    c[CO2Index_]-=deltaC_O2;
    c[cokeIndex_]+=deltaC_O2;
    // Limit the composition
    for (label i=0; i<nTotalSpecies_; i++)
    {
        c[i] = max(0, c[i]);
    }
    if(debug>1)
    {
        Info<<"new c updated: "<<c<<endl;
    }
    
    //solve for new coke fraction
    cokei=c[cokeIndex_]*cokeThermo_.molWeight.value()/cokeThermo_.density.value();
    if(debug>1)
    {
        Info<<"new coke fraction updated: "<<cokei<<endl;
    }
    
    // solve for the new temperature
    Cpf=c[0]*compositions_.Wi(0)*compositions_.Cp(0, pi, Ti); //Species Cp is not related to p and T, but mixture depend on c
    for(label i=1;i<nGasSpecies_;i++)
    {
        Cpf+=c[i]*compositions_.Wi(i)*compositions_.Cp(i, pi, Ti);
    }
    Cps=c[cokeIndex_]*cokeThermo_.molWeight.value()*cokeThermo_.Cp.value()+
        c[rockIndex_]*rockThermo_.molWeight.value()*rockThermo_.Cp.value();

    Ti=(ha-deltaC_O2*hr_)/(Cpf+Cps);
    if(debug>1)
    {
        Info<<"new Ti updated: "<<Ti<<endl;
    }
    
}

inline Foam::scalar Foam::cokeCombustion::solveODEByEulerImplicit(const scalar ak,const scalar ci,const scalar dt,const scalar epsi) const
{
    scalar c_O2_new=(ci/dt)/(ak/epsi+1/dt);
    return c_O2_new;
}

inline Foam::scalar Foam::cokeCombustion::solveODEBy2ndRK(const scalar ak,const scalar ci,const scalar dt,const scalar epsi) const
{
    scalar f1=-ak*(ci/epsi);
    scalar f2=-ak*((ci+dt*f1)/epsi);
    scalar c_O2_new=ci+0.5*dt*(f1+f2);
    return c_O2_new;
}

inline Foam::scalar Foam::cokeCombustion::solveODEBy4thRK(const scalar ak,const scalar ci,const scalar dt,const scalar epsi) const
{
    scalar f1=-ak*(ci/epsi);
    scalar f2=-ak*((ci+dt/2.0*f1)/epsi);
    scalar f3=-ak*((ci+dt/2.0*f2)/epsi);
    scalar f4=-ak*((ci+dt*f3)/epsi);
    scalar c_O2_new=ci+dt/6.0*(f1+2.0*f2+2.0*f3+f4);
    return c_O2_new;
}

inline Foam::scalar Foam::cokeCombustion::solveODEBy4thRKFull(const scalar effssi,const scalar Ti0,const scalar ci0,
    const scalar Cps0,const scalar Cpf0,
    const scalar dt,const scalar epsi) const
{
    scalar ak1=effssi*A_*exp(-Ta_/Ti0);
    scalar f1=-ak1*(ci0/epsi);
    if(debug>1)
    {
        Info<<"4th RK full method [1]: "
            <<" Ti1: "<<Ti0
            <<", ak1: "<<ak1
            <<", f1: "<<f1
            <<endl;
    }

    scalar Ti2=Ti0-(hr_*f1*dt/2.0)/(Cps0+Cpf0);
    scalar ak2=effssi*A_*exp(-Ta_/Ti2);
    scalar f2=-ak2*((ci0+dt/2.0*f1)/epsi);
    if(debug>1)
    {
        Info<<"4th RK full method [2]: "
            <<" Ti2: "<<Ti2
            <<", ak2: "<<ak2
            <<", f2: "<<f2
            <<endl;
    }

    scalar Ti3=Ti0-(hr_*f2*dt/2.0)/(Cps0+Cpf0);
    scalar ak3=effssi*A_*exp(-Ta_/Ti3);
    scalar f3=-ak3*((ci0+dt/2.0*f2)/epsi);
    if(debug>1)
    {
        Info<<"4th RK full method [3]: "
            <<" Ti3: "<<Ti3
            <<", ak3: "<<ak3
            <<", f3: "<<f3
            <<endl;
    }

    scalar Ti4=Ti0-(hr_*f3*dt)/(Cps0+Cpf0);
    scalar ak4=effssi*A_*exp(-Ta_/Ti4);
    scalar f4=-ak4*((ci0+dt*f3)/epsi);
    if(debug>1)
    {
        Info<<"4th RK full method [4]: "
            <<" Ti4: "<<Ti4
            <<", ak4: "<<ak4
            <<", f4: "<<f4
            <<endl;
    }

    scalar c_O2_new=ci0+dt/6.0*(f1+2.0*f2+2.0*f3+f4);

    return c_O2_new;
}




Foam::tmp<Foam::fvScalarMatrix> Foam::cokeCombustion::R(volScalarField& Y) const
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

Foam::tmp<Foam::fvScalarMatrix> Foam::cokeCombustion::Rs(volScalarField& frac) const
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

Foam::tmp<Foam::volScalarField> Foam::cokeCombustion::Qdot() const
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
        Qdot[celli] = -hr_*RRCoke_[celli]/cokeThermo_.molWeight.value();
    }
    return tQdot;
}

void Foam::cokeCombustion::calculate() 
{    
    if (!this->chemistry_)
    {
        return ;
    }

    volScalarField ssArea
    (
        "ssArea",
        2.0*mag(fvc::grad(coke_))*(1-(1-coke_)*(1-coke_))
    );

    // Info<<"mag(fvc::grad(coke_)): "<<mag(fvc::grad(coke_))->field()<<endl;
    // Info<<"ssArea: "<<ssArea.field()<<endl;

    const scalarField& T=thermo_.T();
    const scalarField& rho=thermo_.rho();

    scalar aki=0;
    forAll(T, celli)
    {
        const scalar rhoi = rho[celli];  
        // const scalar epsi = eps_[celli];
        for(label i=0;i<nGasSpecies_;i++)
        {
            // c0_[i]=epsi*rhoi*Y_[i][celli]/compositions_.Wi(i);
            c0_[i]=rhoi*Y_[i][celli]/compositions_.Wi(i);
        }
        aki=ssArea[celli]*A_*std::exp(-Ta_/T[celli]);
        scalar dcdt=aki*c0_[O2Index_];
        RRO2_[celli]=-dcdt*compositions_.Wi(O2Index_);
        RRCO2_[celli]=dcdt*compositions_.Wi(CO2Index_);
        RRCoke_[celli]=-dcdt*cokeThermo_.molWeight.value();        
    }
}


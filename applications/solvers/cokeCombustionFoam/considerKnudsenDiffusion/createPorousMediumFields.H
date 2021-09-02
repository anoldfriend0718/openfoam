    Info<<"Reading the porosity and coke fraction"<<nl<<endl;
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

    volScalarField rEps
    {
        "rEps",
        1.0/(eps+SMALL)
    };

    surfaceScalarField rEpsf
    {
        "rEpsf",
        fvc::interpolate(rEps)
    };

    Info<< "calculating face flux field phi in porous medium"<< nl << endl;
    surfaceScalarField phiByEpsf
    {
        "phiByEpsf",
        phi*rEpsf
    };

    volScalarField rhoByEps("rhoByEps",rho*rEps);


    Info<< "Reading porous medium transportProperties"<< nl << endl;
    IOdictionary porousTransportProperties
    (
        IOobject
        (
            "porousTransportProperties",
            runTime.constant(),
            mesh,
            IOobject::MUST_READ_IF_MODIFIED,
            IOobject::NO_WRITE
        )
    );
   

    Info<< "Createing permeability field"<< nl << endl;

    scalar dp=readScalar(porousTransportProperties.lookup("dp")); 
    
    dimensionedScalar K0 //rock 
    {
        "K0",
        dimensionSet(0, 2, 0, 0, 0, 0, 0),
        porousTransportProperties.lookup("K0")
    };

    dimensionedScalar rK0 //rock
    {
        "rK0",
        dimensionSet(0, -2, 0, 0, 0, 0, 0),
        (1.0/K0).value()
    };

    dimensionedScalar K1 //coke
    {
        "K0",
        dimensionSet(0, 2, 0, 0, 0, 0, 0),
        porousTransportProperties.lookup("K1")
    };

    dimensionedScalar rK1 //coke
    {
        "rK1",
        dimensionSet(0, -2, 0, 0, 0, 0, 0),
        (1.0/K1).value()
    };

    Info<<"coke pore diameter: "<<dp<<endl
        <<"coke permeability factor K0: "<<K1.value()<<endl
        <<"rock permeability factor K0: "<<K0.value() <<endl;

    volScalarField rK
    (
        "rK",
        rK0*(1.0-eps)*(1.0-eps)/max((eps*eps*eps),SMALL)
    );

    //Definition of coke region indicator
    volScalarField cokeRegion
    (
        IOobject
        (
            "cokeRegion",
            runTime.timeName(),
            mesh,
            IOobject::NO_READ,
            IOobject::NO_WRITE
        ),
        mesh,
        dimensionedScalar ("cokeRegion", dimensionSet(0,0,0,0,0,0,0), 0) 
    );

    //label coke regions 
    forAll(coke,celli)
    {
        if(coke[celli]>SMALL)
        {
            cokeRegion[celli]=1.0;
        }
        else
        {
            cokeRegion[celli]=0.0;
        }
    }

    volScalarField::Boundary& cokeBf = coke.boundaryFieldRef();
    volScalarField::Boundary& cokeRegionBf=cokeRegion.boundaryFieldRef();
    forAll(cokeBf,patchi)
    {
        forAll(cokeBf[patchi],facei)
        {
            if(cokeBf[patchi][facei]>SMALL)
            {
                cokeRegionBf[patchi][facei]=1.0;
            }
            else
            {
                cokeRegionBf[patchi][facei]=0.0;
            }
        }
    }

    forAll(rK,celli)
    {
        if(cokeRegion[celli]>1.0-SMALL) //==1.0
        {
            scalar epsi=eps[celli];
            rK[celli]=(1.0-epsi)*(1.0-epsi)/max(epsi*epsi*epsi,SMALL)*rK1.value();
        }
    }


    volScalarField rKnano
    (
        IOobject
        (
            "rKnano",
            runTime.timeName(),
            mesh,
            IOobject::NO_READ,
            IOobject::NO_WRITE
        ),
        mesh,
        dimensionedScalar(dimless,0.0)
    );

    {
        const volScalarField& mu=thermo.mu()();
        const volScalarField& T=thermo.T();
        const scalar MN2=composition.Wi(N2Index);
        const scalar MCO2=composition.Wi(CO2Index);
        const scalar M=2*MN2*MCO2/(MN2+MCO2);
                    
        forAll(cokeRegion,i)
        {
            if(cokeRegion[i]>1-SMALL)
            {
                scalar mui=mu[i];
                scalar pi=p[i];
                scalar Ti=T[i];
                scalar kni=2*(mui/pi*std::sqrt(PI*Rg1000*Ti/(2*M)))/dp;
                scalar alphai=128/(15*PI2)*std::atan(4*std::pow(kni,0.4));
                scalar rKnanoi=1./((1+alphai*kni)*(1+4*kni/(1+kni)));
                rKnano[i]=rKnanoi;
            }
            else
            {
                rKnano[i]=1.0;
            }
        }
    }

    Info<<"Creating the porous medium drag force field"<< nl <<endl;
    volScalarField drag
    {
        "drag",
        fvc::average(mu*rK*rKnano)
    };

    // Info<<"Label the cells that have porous medium"<<endl;
    // Definition of Solid Indicator
    volScalarField solid
    (
        IOobject
        (
            "solid",
            runTime.timeName(),
            mesh,
            IOobject::NO_READ,
            IOobject::NO_WRITE
        ),
        mesh,
        dimensionedScalar ("solid", dimensionSet(0,0,0,0,0,0,0), 0) 
    );


    //label solid regions 
    forAll(eps,celli)
    {
        if(eps[celli]>0.99)
        {
            solid[celli]=0.0;
        }
        else
        {
            solid[celli]=1.0;
        }
    }

    forAll(drag,celli)
    {
        if(solid[celli]<small) //==0
        {
            drag[celli]=0.0;
        }
    }
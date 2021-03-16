Info<<"Reading the porosity"<<nl<<endl;
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
    1-eps
);

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

dimensionedScalar K0
{
    "K0",
    dimensionSet(0, 2, 0, 0, 0, 0, 0),
    porousTransportProperties.lookup("K0")
};

dimensionedScalar rK0
{
    "rK0",
    dimensionSet(0, -2, 0, 0, 0, 0, 0),
    (1.0/K0).value()
};

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


Info<<"Label the cells that have porous medium"<<endl;
// Definition of Solid Indicator
volScalarField solid
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
    dimensionedScalar ("solid", dimensionSet(0,0,0,0,0,0,0), 0) 
);

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

volScalarField::Boundary& epsBf = eps.boundaryFieldRef();
volScalarField::Boundary& solidBf=solid.boundaryFieldRef();
forAll(epsBf,patchi)
{
    forAll(epsBf[patchi],facei)
    {
        if(epsBf[patchi][facei]>0.99)
        {
            solidBf[patchi][facei]=0.0;
        }
        else
        {
            solidBf[patchi][facei]=1.0;
        }
    }
}

Info<< "Createing permeability field"<< nl << endl;
volScalarField rK
(
    "rK",
    rK0*(1.0-eps)*(1.0-eps)/max((eps*eps*eps),SMALL)
);

Info<< "calculating face flux field phi in porous medium"<< nl << endl;
surfaceScalarField phiByEpsf
{
    "phiByEpsf",
    phi*rEpsf
};

volScalarField rhoByEps("rhoByEps",rho*rEps);

Info<<"Creating the porous medium drag force field"<< nl <<endl;
volScalarField drag
{
    "drag",
    fvc::average(mu*rK)
};


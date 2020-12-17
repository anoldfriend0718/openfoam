Info<< "\nReading field porosity" << endl;
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

volScalarField tEps(eps+SMALL);

//Permeability Models
Info<< "Reading intrinsic permeability field field K\n" << endl;
dimensionedScalar K0
{
    "K0",
    dimensionSet(0, 2, 0, 0, 0, 0, 0),
    transportProperties.lookup("K0")
};

dimensionedScalar rK0
{
    "rK0",
    dimensionSet(0, -2, 0, 0, 0, 0, 0),
    (1.0/K0).value()
};

Info<< "create permeability field\n" << endl;

volScalarField rK
(
    "rK",
    rK0*(max((1-eps)*(1-eps),SMALL))/(tEps*tEps*tEps)
);

// Foam::volScalarField K
// (
//     "K",
//     1./rK
// );

// K.write();

surfaceScalarField repsf
(
    "repsf",
    fvc::interpolate(1/tEps)
);

surfaceScalarField phiByEpsf
(
    "phiByEpsf",
     phi*repsf
);

Info<< "Creating Solid Indicator\n" << endl;
// Defenition of Solid Indicator
volScalarField Solid
(
    IOobject
    (
        "Solid",
        runTime.timeName(),
        mesh,
        IOobject::NO_READ,
        IOobject::AUTO_WRITE
    ),
    mesh,
    dimensionedScalar ("solid", dimensionSet(0,0,0,0,0,0,0), 0) 
);


//Labeling the places that have porous media
forAll(eps,celli) 
{
	if(eps[celli]>0.99) 
	{	
		Solid[celli]=0; 
	}
	else
	{
		Solid[celli]=1; 
	}
}

forAll(mesh.boundary(), patchi) 
{
	forAll(Solid.boundaryField()[patchi],facei)
	{
		if(eps.boundaryField()[patchi][facei]>0.99) 
		{
			Solid.boundaryFieldRef()[patchi][facei]=0;
		}
		else
		{
			Solid.boundaryFieldRef()[patchi][facei]=1;
		}
	}
}
Solid.write();
surfaceScalarField Solidf(fvc::interpolate(Solid));

//Eliminating intermediate values of Solid Indicator Function
forAll(Solidf,facei)
{
	if(Solidf[facei]<1) // if <1 ==0.  
	{	
		Solidf[facei]=0; 
	}
}

//Drag Coefficient Calculation
volScalarField drag ("drag", nu*rK);

forAll(drag,celli)
{
	if (Solid[celli] == 0)
	{
		drag[celli]= 0;
	}
}



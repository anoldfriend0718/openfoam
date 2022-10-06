Info<< "Creating coke combustion reaction model"<<nl << endl;
cokeCombustion reaction(mesh,thermo,cokeThermo,rockThermo);

volScalarField rhoCoke
(
    IOobject
    (
        "rhoCoke",
        runTime.timeName(),
        mesh,
        IOobject::NO_READ,
        IOobject::NO_WRITE
    ),
    mesh,
    cokeThermo.density
);

//just for outputing coke reaction rate (kg/m3/s)
volScalarField cokeRectionRate
(
    IOobject
    (
        "cokeRectionRate",
        runTime.timeName(),
        mesh,
        IOobject::READ_IF_PRESENT,
        IOobject::AUTO_WRITE
    ),
    mesh,
    dimensionedScalar("Ri", dimensionSet(1, -3, -1, 0, 0, 0, 0), 0)
);

//just for outputing combustion heat rate (kg/m3/s)
volScalarField Qdot
(
    IOobject
    (
        "Qdot",
        runTime.timeName(),
        mesh,
        IOobject::NO_READ,
        IOobject::AUTO_WRITE
    ),
    mesh,
    dimensionedScalar("Qdot", dimEnergy/dimVolume/dimTime, 0)
);


// const tmp<volScalarField>& trRhoCoke=1.0/rhoCoke;
// const volScalarField& rRhoCoke=trRhoCoke.ref();

// volScalarField RRg
// (
//     IOobject
//     (
//         "RRg",
//         runTime.timeName(),
//         mesh,
//         IOobject::NO_READ,
//         IOobject::NO_WRITE
//     ),
//     mesh,
//     dimensionedScalar("Ri", dimensionSet(1, -3, -1, 0, 0, 0, 0), 0)
// );

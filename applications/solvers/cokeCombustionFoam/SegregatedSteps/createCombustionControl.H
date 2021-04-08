
IOdictionary combustionControl
(
    IOobject
    (
        "combustionControl",
        runTime.constant(),
        mesh,
        IOobject::MUST_READ_IF_MODIFIED,
        IOobject::NO_WRITE
    )
);

scalar targetResidualCokeFraction=combustionControl.lookupOrDefault("targetResidualCokeFraction",1.0);
scalar residualCokeFrac=(fvc::domainIntegrate(coke)/gSum(mesh.V())).value();
Info<<"Target residual coke fraction: "<<targetResidualCokeFraction<<endl;
Info<<"Initial coke fraction: "<<residualCokeFrac<<endl;
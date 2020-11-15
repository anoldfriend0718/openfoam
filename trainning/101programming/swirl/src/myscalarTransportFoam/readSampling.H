Info<<"Reading sampling properties"<<endl;
Foam::IOdictionary logProperties
(
    Foam::IOobject
    (
        "samplingProperties",
        runTime.system(),
        mesh,
        Foam::IOobject::MUST_READ_IF_MODIFIED,
        Foam::IOobject::NO_WRITE
    )
);

const word fieldName {logProperties.lookup("field")};
const label cellIndex=Foam::readLabel(logProperties.lookup("cellIndex"));

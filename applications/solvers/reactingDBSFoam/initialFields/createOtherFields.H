Info<<"==Create extra property Fields=="<<nl<<endl;
//heat capacity ratio: solid Cp / fluid Cp pr Cv based on specific energy model 
volScalarField CpsByCpvf
{
    IOobject
    (
        "CpsByCpvf",
        runTime.timeName(),
        mesh,
        IOobject::NO_READ,
        IOobject::NO_WRITE
    ),
    Cps/thermo.Cpv()
};
// CpsByCpvf.write();


// construct effective thermal diffusivity for energy by 
//alphaEff=(1-eps)*alphas/Cpvf+eps*alphahe
//alphahe: fluid Thermal diffusivity for energy 
//alphahe=alpha*CpByCpv= ladal/Cpf * (Cpf/Cpvf)
volScalarField alphaEff
{
    IOobject
    (
        "alphaEff",
        runTime.timeName(),
        mesh,
        IOobject::NO_READ,
        IOobject::NO_WRITE
    ),
    eps*thermo.alphahe()+(1-eps)*ladals/thermo.Cpv()
};

Info<<"create rhofEps "<<endl;
Foam::volScalarField rhofEps
{
    IOobject
    (
        "rhofEps",
        runTime.timeName(),
        mesh,
        IOobject::NO_READ,
        IOobject::NO_WRITE
    ),
    eps*rhof
};
Info<<"store rhofEps old field "<<endl;
rhofEps.oldTime();

Info<<"create effRhosCpsByCpvf"<<endl;
Foam::volScalarField effRhosCpsByCpvf
{
    IOobject
    (
        "effRhosCpsByCpvf",
        runTime.timeName(),
        mesh,
        IOobject::NO_READ,
        IOobject::NO_WRITE
    ),
    rhos*(1-eps)*CpsByCpvf
};
Info<<"store effRhosCpsByCpvf old field "<<endl;
effRhosCpsByCpvf.oldTime();


Info<<"==Complete creating extra property Fields=="<<nl<<endl;
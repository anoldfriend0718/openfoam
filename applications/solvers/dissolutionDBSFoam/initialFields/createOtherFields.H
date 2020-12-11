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


Info<<"create rhosCpsByCpvf"<<endl;
Foam::volScalarField rhosCpsByCpvf
{
    IOobject
    (
        "rhosCpsByCpvf",
        runTime.timeName(),
        mesh,
        IOobject::NO_READ,
        IOobject::NO_WRITE
    ),
    rhos*CpsByCpvf
};
Info<<"store rhosCpsByCpvf old field "<<endl;
rhosCpsByCpvf.oldTime();


Info<<"==Complete creating extra property Fields=="<<nl<<endl;
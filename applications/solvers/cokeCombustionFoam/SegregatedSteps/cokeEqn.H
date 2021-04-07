
//solving coke volume evolution equation
fvScalarMatrix cokeEqn
(
    fvm::ddt(rhoCoke,coke)
  ==
    reaction.Rs(coke) //note signs
  + fvOptions(rhoCoke,coke)
);
cokeEqn.relax();
fvOptions.constrain(cokeEqn);

cokeEqn.solve();
fvOptions.correct(coke);
coke=max(min(coke,1.0),0.0);

Info<<"updating the porous medium and related fields"<<endl;
eps=1-coke-rock;
eps.max(0.0);
rEps=1.0/(eps+SMALL);
rEpsf=fvc::interpolate(rEps);
phiByEpsf=phi*rEpsf;

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

rK=rK0*(1.0-eps)*(1.0-eps)/max((eps*eps*eps),SMALL);
drag=fvc::average(mu*rK);
forAll(drag,celli)
{
    if(solid[celli]<small) //==0
    {
        drag[celli]=0.0;
    }
}

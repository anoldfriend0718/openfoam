{
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

    forAll(coke,celli)
    {
        //remove the residual coke,the threshold 1e-6 defined in the specialSurfaceArea class
        if(coke[celli]<=1e-6) 
        {
            coke[celli]=0.0;
        }
    }

    Info<<"updating the porous medium and related fields"<<endl;
    eps=1-coke-rock;
    // eps.max(0.0); //bug in HPC
    rEps=1.0/(eps+SMALL);
    rEpsf=fvc::interpolate(rEps);
    phiByEpsf=phi*rEpsf;

    //reupdate the coke region indication
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
    //calculate the rK
    rK=rK0*(1.0-eps)*(1.0-eps)/max((eps*eps*eps),SMALL);
    forAll(cokeRegion,celli)
    {
        if(cokeRegion[celli]>1.0-SMALL) //==1.0
        {
            scalar epsi=eps[celli];
            rK[celli]=(1.0-epsi)*(1.0-epsi)/max(epsi*epsi*epsi,SMALL)*rK1.value();
        }

    }

    // //calculate the rknano
    // const volScalarField& mu=thermo.mu()();
    // const volScalarField& T=thermo.T();
    // const scalar MN2=composition.Wi(N2Index);
    // const scalar MCO2=composition.Wi(CO2Index);
    // const scalar M=2*MN2*MCO2/(MN2+MCO2);
    // forAll(cokeRegion,i)
    // {
    //     if(cokeRegion[i]>1-SMALL)
    //     {
    //         scalar mui=mu[i];
    //         scalar pi=p[i];
    //         scalar Ti=T[i];
    //         scalar kni=2*(mui/pi*std::sqrt(PI*Rg1000*Ti/(2*M)))/dp;
    //         scalar alphai=128/(15*PI2)*std::atan(4*std::pow(kni,0.4));
    //         scalar rKnanoi=1./((1+alphai*kni)*(1+4*kni/(1+kni)));
    //         rKnano[i]=rKnanoi;
    //     }
    //     else
    //     {
    //         rKnano[i]=1.0;
    //     }
    // }
    // drag=fvc::average(mu*rK*rKnano);

    drag=fvc::average(mu*rK);
    
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

}
{
    volScalarField& he = thermo.he();
    //update the solid Cp/fluid Cpv
    CpsByCpvf= Cps/thermo.Cpv();
    //update the effective thermal diffusivity for energy 
    alphaEff=eps*thermo.alphahe()+(1-eps)*ladals/thermo.Cpv();

    //store the old rhofEps and effRhosCpsByCpvf due to ddt computations
    rhofEps.storeOldTimes();
    effRhosCpsByCpvf.storeOldTimes();
    //update the rhofEps and effRhosCpsByCpvf
    rhofEps=rhof*eps;
    effRhosCpsByCpvf=rhos*(1-eps)*CpsByCpvf;
    mdotSolid=ssArea*reactionRate*rhof*(Y-Yeq)*beta;

    fvScalarMatrix EEqn
    (
       fvm::ddt(rhofEps, he) + fvm::div(phi, he)
         //fvm::ddt(eps,rhof, he) + fvm::div(phi, he)
      + fvm::ddt(effRhosCpsByCpvf, he)
      + fvc::ddt(rhof, K) + fvc::div(phi, K)
      + (
            he.name() == "e"
          ? fvc::div
            (
                fvc::absolute(phi/fvc::interpolate(rhof), U),
                p,
                "div(phiv,p)"
            )
          : -dpdt
        )
	    - fvm::laplacian(alphaEff, he)
     ==
        rhof*(U&g)
      + fvm::Su(mdotSolid*reactionHeat,he)
      + radiation->Sh(thermo, he)
      + fvOptions(rhof, he)
    );

    EEqn.relax();

    fvOptions.constrain(EEqn);

    EEqn.solve();

    fvOptions.correct(he);

    thermo.correct();
    radiation->correct();
}

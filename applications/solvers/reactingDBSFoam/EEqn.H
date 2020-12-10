{
    volScalarField& he = thermo.he();
    //update the solid Cp/fluid Cpv
    CpsByCpvf= Cps/thermo.Cpv();
    //update the effective thermal diffusivity for energy 
    alphaEff=eps*thermo.alphahe()+epsSolid*ladals/thermo.Cpv();

    //store the old rhosCpsByCpvf due to ddt computations
    rhosCpsByCpvf.storeOldTimes();
    //update  rhosCpsByCpvf
    rhosCpsByCpvf=rhos*CpsByCpvf;

    mdotSolid=ssArea*reactionRate*rhof*(Y-Yeq)*beta;

    fvScalarMatrix EEqn
    (
       fvm::ddt(eps,rhof, he) + fvm::div(phi, he)
      + fvm::ddt(epsSolid,rhosCpsByCpvf, he)
      // + fvc::ddt(rhof, K) + fvc::div(phi, K)
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

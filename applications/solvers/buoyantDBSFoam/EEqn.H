{
    volScalarField& he = thermo.he();
    //update the solid Cp/fluid Cpv
    CpsByCpvf= Cps/thermo.Cpv();
    //update the effective thermal diffusivity for energy 
    alphaEff=eps*thermo.alphahe()+eps_solid*ladals/thermo.Cpv();

    //store the old rhofEps and effRhosCpsByCpvf due to ddt computations
    rhosCpsByCpvf.storeOldTimes();
    //update the rhofEps and effRhosCpsByCpvf
    rhosCpsByCpvf=rhos*CpsByCpvf;

    fvScalarMatrix EEqn
    (
        fvm::ddt(eps,rhof, he) + fvm::div(phi, he)
      + fvm::ddt(eps_solid,rhosCpsByCpvf, he)
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

{
    volScalarField& he = thermo.he();
    //update the solid Cp/fluid Cpv
    CpsByCpvf= Cps/thermo.Cpv();
    //update the effective thermal diffusivity for energy 
    alphaEff=eps*thermo.alphahe()+(1-eps)*ladals/thermo.Cpv();
    
    fvScalarMatrix EEqn
    (
        fvm::ddt(rhof*eps, he) + fvm::div(phi, he)
      + fvc::ddt(rhof, K) + fvc::div(phi, K)
      + fvm::ddt(rhos*(1-eps)*CpsByCpvf, he)
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

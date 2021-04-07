{
    volScalarField& he = thermo.he();

    cokeRhoCpByCpvf.storeOldTimes();
    cokeRhoCpByCpvf=cokeThermo.density*cokeThermo.Cp/thermo.Cpv();
    
    rockRhoCpByCpvf.storeOldTimes();
    rockRhoCpByCpvf=rockThermo.density*rockThermo.Cp/thermo.Cpv();

    tmp<volScalarField> talphaEff = eps * thermo.alphahe()
                                    + coke* cokeThermo.kappa/thermo.Cpv()
                                    + rock* rockThermo.kappa/thermo.Cpv();
    const volScalarField& alphaEff = talphaEff();    
    
    //only suitable for the const cpv    
    fvScalarMatrix EEqn
    (
          fvm::ddt(eps, rho,            he) 
        + fvm::ddt(coke,cokeRhoCpByCpvf,he)
        + fvm::ddt(rock,rockRhoCpByCpvf,he)
        // + fvm::div(phi, he)
        + mvConvection->fvmDiv(phi, he)
        + fvc::ddt(rho, K) + fvc::div(phi, K)
        + (
                he.name() == "e"
            ? fvc::div
                (
                    fvc::absolute(phi/fvc::interpolate(rho), U),
                    p,
                    "div(phiv,p)"
                )
            : -dpdt
            )
        - fvm::laplacian(alphaEff, he)
        ==
          rho*(U&g)
        + fvm::Su(reaction.Qdot(),he)
        + fvOptions(rho, he)
    );   

    EEqn.relax();

    fvOptions.constrain(EEqn);

    EEqn.solve();

    fvOptions.correct(he);

    thermo.correct();

    const volScalarField& T=thermo.T();

    // Info<< "min/max fvm::ddt(eps, rho,            he)  = "
    //     << min(fvm::ddt(eps, rho,            he)&he).value() << ", " << max(fvm::ddt(eps, rho,            he)&he).value() << endl;    

    // Info<< "min/max Qdot = "
    //     << min(reaction.Qdot()).value() << ", " << max(reaction.Qdot()).value() << endl;    

    Info<< "min/max(T) = "
        << min(T).value() << ", " << max(T).value() << endl;    
}
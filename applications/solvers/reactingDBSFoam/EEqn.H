{
    volScalarField& he = thermo.he();
    //update the solid Cp/fluid Cpv
    CpsByCpvf= Cps/thermo.Cpv();
    //update the effective thermal diffusivity for energy 
    alphaEff=eps*thermo.alphahe()+(1-eps)*ladals/thermo.Cpv();

    // rhofEps.storeOldTimes();
    rhofEps=rhof*eps;

    // if(runTime.timeName()=="1")
    // {
    //     Info<< "Time = " << runTime.timeName() << endl;

    //     Info << "[QXU mark]: eps"<< endl; 
    //     Info << eps <<endl;
        
    //     Info << "[QXU mark]: he"<< endl; 
    //     Info << he <<endl; 


    //     Info << "[QXU mark]: rhof"<< endl; 
    //     Info << rhof <<endl;


    //     Info << "[QXU mark]: rhofEps"<< endl; 
    //     Info << rhofEps <<endl;

    //     Info << "[QXU mark]: ddtRhof"<< endl;
    //     fvScalarMatrix ddtRhof(fvm::ddt(rhof, he));
    //     Info<<"ddtRhof = "<<ddtRhof<<endl;
    //     Info<<"ddtRhof.lduAddr().lowerAddr() = "<<ddtRhof.lduAddr().lowerAddr()<<endl;
    //     Info<<"ddtRhof.lduAddr().upperAddr() = "<<ddtRhof.lduAddr().upperAddr()<<endl;


    //     Info << "[QXU mark]: ddtRhofEps1"<< endl;
    //     fvScalarMatrix ddtRhofEps(fvm::ddt(rhof*eps, he));
    //     Info<<"ddtRhofEps1 = "<<ddtRhofEps<<endl;
    //     Info<<"ddtRhofEps1.lduAddr().lowerAddr() = "<<ddtRhofEps.lduAddr().lowerAddr()<<endl;
    //     Info<<"ddtRhofEps1.lduAddr().upperAddr() = "<<ddtRhofEps.lduAddr().upperAddr()<<endl;


    //     Info << "[QXU mark]: ddtRhofEps2"<< endl;
    //     fvScalarMatrix ddtRhofEps2(fvm::ddt(rhofEps, he));
    //     Info<<"ddtRhofEps2 = "<<ddtRhofEps2<<endl;
    //     Info<<"ddtRhofEps2.lduAddr().lowerAddr() = "<<ddtRhofEps2.lduAddr().lowerAddr()<<endl;
    //     Info<<"ddtRhofEps2.lduAddr().upperAddr() = "<<ddtRhofEps2.lduAddr().upperAddr()<<endl;

    //     Info << "[QXU mark]: ddtRhofEps3"<< endl;
    //     fvScalarMatrix ddtRhofEps3(fvm::ddt(eps,rhof, he));
    //     Info<<"ddtRhofEps3 = "<<ddtRhofEps3<<endl;
    //     Info<<"ddtRhofEps3.lduAddr().lowerAddr() = "<<ddtRhofEps3.lduAddr().lowerAddr()<<endl;
    //     Info<<"ddtRhofEps3.lduAddr().upperAddr() = "<<ddtRhofEps3.lduAddr().upperAddr()<<endl;

    // }


    fvScalarMatrix EEqn
    (
       fvm::ddt(rhofEps, he) + fvm::div(phi, he)
         //fvm::ddt(eps,rhof, he) + fvm::div(phi, he)
      +  fvm::ddt(rhos*(1-eps)*CpsByCpvf, he)
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

{
    PtrList<volScalarField> rhoDeffs(Y.size());
    rhoDeffs.set(O2Index,new volScalarField("rhoEps2DO2",rho*eps*eps*DO2));
    rhoDeffs.set(CO2Index,new volScalarField("rhoEps2DCO2",rho*eps*eps*DCO2));

    volScalarField Yt(0.0*Y[0]);

    //solving O2 advection-diffusion equation
    forAll(Y,i)
    {        
        if(Y[i].name()!="N2")
        {
            volScalarField& Yi = Y[i];
            volScalarField& rhoDeff=rhoDeffs[i];

            fvScalarMatrix YiEqn
            (
                  fvm::ddt(eps,rho, Yi)
                + mvConvection->fvmDiv(phi, Yi)
                - fvm::laplacian(rhoDeff, Yi)
              ==
                  reaction.R(Yi)
                + fvOptions(rho, Yi)
            );

            YiEqn.relax();

            fvOptions.constrain(YiEqn);

            YiEqn.solve("Yi");

            fvOptions.correct(Yi);

            Yi.max(0.0);
            Yi.min(1.0);

            Yt += Yi;

            if(Y[i].name()=="O2")
            {
                // Info<<"update the reaction rate after solving O2 transfer equation"<<endl;
                reaction.correct();
                cokeRectionRate=reaction.Rs(coke) & coke;
                Qdot=reaction.Qdot().ref();
            }

        }
    }

    Y[N2Index] = scalar(1) - Yt;
    Y[N2Index].max(0.0);
}
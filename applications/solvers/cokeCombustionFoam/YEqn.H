tmp<fv::convectionScheme<scalar>> mvConvection
(
    fv::convectionScheme<scalar>::New
    (
        mesh,
        fields,
        phi,
        mesh.divScheme("div(phi,Yi_h)")
    )
);

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
        }
    }

    Y[N2Index] = scalar(1) - Yt;
    Y[N2Index].max(0.0);
}
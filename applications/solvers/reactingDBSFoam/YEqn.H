{
    eps2Df=eps*eps*Df;
    mdotY=ssArea*reactionRate*rhof;

    fvScalarMatrix YEqn
    (
        fvm::ddt(eps,rhof, Y)
        + fvm::div(phi,Y)
        - fvm::laplacian(eps2Df,Y)
        ==
        - fvm::Sp(mdotY,Y)
        + fvm::Su(mdotY*Yeq,Y)
        + fvOptions(rhof, Y)
    );

    YEqn.relax();

    fvOptions.constrain(YEqn);

    YEqn.solve();

    fvOptions.correct(Y);

    Y.max(0.0);
    Y.min(1.0);
}
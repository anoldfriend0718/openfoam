// {
//     volScalarField& he = thermo.he();

//     CpsByCpvf= Cps/thermo.Cpv();

//     fvScalarMatrix EEqn
//     (
//         fvm::ddt(rho*eps, he) + fvm::div(phi, he)
//       + fvc::ddt(rho, K) + fvc::div(phi, K)
//       + 
//       + (
//             he.name() == "e"
//           ? fvc::div
//             (
//                 fvc::absolute(phi/fvc::interpolate(rho), U),
//                 p,
//                 "div(phiv,p)"
//             )
//           : -dpdt
//         )
//       - fvm::laplacian(turbulence->alphaEff(), he)
//      ==
//         rho*(U&g)
//       + radiation->Sh(thermo, he)
//       + fvOptions(rho, he)
//     );

//     EEqn.relax();

//     fvOptions.constrain(EEqn);

//     EEqn.solve();

//     fvOptions.correct(he);

//     thermo.correct();
//     radiation->correct();
// }

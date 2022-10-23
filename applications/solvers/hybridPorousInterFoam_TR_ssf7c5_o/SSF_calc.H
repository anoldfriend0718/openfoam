// C_sh read from fvsolution>PIMPLE>Cpc
// Surface force
// Cpc should be 0.5 for dynamic problems, > 0.9 for static problems
    //scalar Cpc (readScalar      (     alpha1.mesh().solutionDict().subDict("PIMPLE").lookup("Cpc")        )      );
    //scalar Cpc(transportProperties.lookupOrDefault("interfaceSharpenCoeff",0.5)); 
    // scalar gradAlphaThreshold(transportProperties.lookupOrDefault("gradAlphaThreshold",1.0));
    //- Sharp force coefficient (put 0.98-0.99 for static problems, 0.4-0.5 for dynamic)
//     if (CorrGradAlpha == 1) 
//     {
//         //Fluidf = (max(1.0-Solidf,SMALL));
//         alpha1_i = fvc::average(fvc::interpolate(max(alpha1,SMALL),"Fluid_c2f")*Fluidf)/fvc::average(Fluidf);
//         alpha1_i = (1-Solid)*alpha1 + Solid*alpha1_i;
//         //alpha1_i2 = 0.5*alpha1_i + 0.5*alpha1;
//         //alpha1_w = fvc::average(fvc::interpolate(max(alpha1,SMALL),"Fluid_c2f")*Fluidwf)/fvc::average(Fluidwf);
//         // forAll(alpha1_i,celli) 
//         // {
//         //     if(alpha1_i[celli]<0.05) 
//         //     {	
//         //         alpha1_i[celli]=0; 
//         //     }
//         // }
//         // alpha1 = (1-Solidw)*alpha1 + Solidw*alpha1_w;
//     }
//     else
//     {
//         alpha1_i = alpha1;
//     }

//     scalar cPc
//     (  readScalar ( alpha1.mesh().solverDict(alpha1.name()).lookup("cPc") ) );

// // eq 23
// // Sharpen interface function
//     // volScalarField alpha_pc = (1.0/(1.0-cPc))*(min( max(alpha1,cPc/2.0), (1.0 -cPc/2.0) ) - cPc/2.0);
//     volScalarField alpha_pc = (1.0/(1.0-cPc))*(min( max(alpha1_i,cPc/2.0), (1.0 -cPc/2.0) ) - cPc/2.0);
    // surfaceScalarField deltasf = fvc::snGrad(alpha_pc);   //fvc::interpolate(Solidv);

    // // surfaceScalarField deltasf = fvc::snGrad(alpha_pc);

    // forAll(deltasf,facei) 
    //     {
    //         if(deltasf[facei] <= gradAlphaThreshold_)   //will cause non-physics fc
    //         {	
    //             deltasf[facei]=0.0; 
    //         }
    //     }

    alpha1_i =  mixture.alphapc();

//eq 22
    //surfaceScalarField fcf = mixture.sigmaKSSF()*deltasf;
    fc = mixture.fc();
    // fc2 = mixture.fc0();
    // thetafc = mixture.thetafc();
    // thetafc0 = mixture.thetafc0();
    Klg = mixture.Klg();
    // Klg0 = mixture.Klg0();
    nlg = mixture.nlg();

    // Knx = mixture.knx();
    // thetafc0 = mixture.thetafc0();
// //Produce fc on cell centers
//     fc = fvc::average(fcf*Solidf0*mesh.Sf()/mesh.magSf());   //(1.-Solidf)

    gradAlpha1 = fvc::grad(alpha1, "nHat"); 
    gradAlphapc1 = fvc::reconstruct(mixture.deltasf()*mesh.magSf());

    // volScalarField gradAlpha_mag(mag(gradAlphapc1));  

    // volScalarField gradAlpha_nfilt(0.0*alpha1);  


    // forAll(gradAlpha_mag,celli) 
    // {
    //     if(gradAlpha_mag[celli] <= 1.0)   //will cause non-physics fc
    //     {	
    //         gradAlpha_nfilt[celli]=0.0; 
    //     }
    //     else
    //     {
    //         gradAlpha_nfilt[celli]=1.0; 
    //     }
  
    // }
    // gradAlphapc2 *= gradAlpha_nfilt;

//     gradAlphapc2 = fvc::reconstruct(fvc::snGrad(alpha_pc)*mesh.magSf());
//     fc2 = fvc::reconstruct(fcf*Solidf0* mesh.magSf());

//     scalar SmootherCount = 0;
//     if (SmoothingMethod == 1)
//     {
//         volVectorField gradSolid(fvc::grad(Solid, "nHatSolid"));
//         // surfaceVectorField gradSolidf(fvc::interpolate(gradSolid));
//         nHatSolidfv0 = (gradSolid/(mag(gradSolid) + deltaN));
//         volScalarField nHatSolidfvMag0(mag(nHatSolidfv0));
//         volScalarField nHatSolidfvMag(max(nHatSolidfvMag0,deltaN.value()));
//         nHatSolidfv1 = nHatSolidfv0;
//         nHatSolidfv2 = nHatSolidfv0;
//         while (SmootherCount <= SmoothingCycles)
//         {
//         // Smoothing the normal of the solid surface
//             nHatSolidfv1 = fvc::average(fvc::interpolate(nHatSolidfv1*nHatSolidfvMag))/fvc::average(fvc::interpolate(nHatSolidfvMag));
//             nHatSolidfv2 = nHatSolidfvMag0*fvc::average(fvc::interpolate(nHatSolidfv2*nHatSolidfvMag))/fvc::average(fvc::interpolate(nHatSolidfvMag));
//             nHatSolidfv2 = nHatSolidfv2/(mag(nHatSolidfv2) + deltaN.value());
//             SmootherCount +=1; 
//         }

//     }


/*
    alpha_pc = 1.0/(1.0-cPc_)*(min( max(alpha1_,cPc_/2.0), (1.0-cPc_/2.0) ) - cPc_/2.0);

    alpha_pc.correctBoundaryConditions();

    surfaceScalarField deltasf = fvc::snGrad(alpha_pc);

    //surface tension force
    surfaceScalarField stf = fvc::interpolate(sigmaK())*deltasf;


    //surface tension force flux
    phic_ = stf*magSf;
*/

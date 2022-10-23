    Info<< "solve CEqn" << endl;
    // volScalarField eps2D (eps*eps*Df);

    //volScalarField rhoEpsDa (rho*eps*eps*Da1*Da2/(alpha1*Da2 + alpha2*Da1));
    //volScalarField DaCST ("DaCST", -rhoEpsDa*C*(rho2*Ha-rho1)/(rho1*Ha*alpha1+rho2*alpha2));    // - _v2
    //volScalarField suCST ("suCST", fvc::laplacian(DaCST,alpha1)); 

    //volScalarField Solidv ("Solidv",max((1.0-Solid),SMALL));
    //Info<< "flag1" << endl;
    rhoEpsDa = (eps*eps*Da1*Da2/(alpha1*Da2 + alpha2*Da1))*Solidv; //(max((1.0-Solid),SMALL));
    //DaCST =  -rhoEpsDa*C*(rho2*Ha-rho1)/(rho1*Ha*alpha1+rho2*alpha2);    // - _v2
    //Info<< "flag2" << endl;
    DaCST =  - rhoEpsDa*(1.0-Ha)/(alpha1 + Ha*alpha2); 
    //Info<< "flag3" << endl;
    suCST = fvc::laplacian(DaCST*C,alpha1);   //for debug
    //gradAlpha = fvc::grad(alpha1);
    //Info<< "flag4" << endl;
    //CST source term using snGrad(Alpha1)
    PhiGradAlphaf=mesh.magSf()*fvc::snGrad(alpha1);  //surfaceScalarField
    //PhiGradAlphaf2=PhiGradAlphaf;  //surfaceScalarField
    //Info<< "flag5" << endl;
/*
    // filter phi near solid
    forAll(PhiGradAlphaf2,facei) 
    {
	if(Solidf0[facei]>0) 
	{	
	PhiGradAlphaf2[facei]=0; 
	}
    }
*/
    suCST_ga=fvc::div(PhiGradAlphaf,DaCST*C);  //for debug
    // Info<< "flag6" << endl;
    suCST_err=suCST_ga - suCST;   //for debug
    //Info<< "flag7" << endl;
    //surfaceScalarField PhiGradAlphaCSTf ("PhiGradAlphaCSTf",fvc::interpolate(Solidv)*fvc::interpolate(DaCST)*PhiGradAlphaf);   // fvc::interpolate(some_vsf, "harmonic")
    PhiGradAlphaCSTf = Solidf0*fvc::interpolate(DaCST)*PhiGradAlphaf;  // fvc::interpolate(Solidv)
    //Info<< "flag8" << endl;
    suCST_ga2=fvc::div(PhiGradAlphaCSTf,C);
    //Info<< "flag9" << endl;
    fvScalarMatrix CEqn
    (
        fvm::ddt(eps,C) + fvm::div(phi, C)
     ==
        fvm::laplacian(rhoEpsDa,C) 
      //+ fvm::SuSp(suCST,C)
      //- fvm::SuSp(-suCST,C) 
      //+ fvm::Su(suCST_ga2,C) 
        + fvm::div(PhiGradAlphaCSTf, C)
        - fvm::SuSp(-mdotci*(C-Ceq),C)      //reaction Sp  

        // + fvm::Sp(mdotci,C)      //reaction Sp  
        // - fvm::Su(mdotci*Ceq,C)  //reaction Su
    );

    CEqn.relax();
    
    CEqn.solve();
    // C = max(C,0.0);
    Info << "CEqn" << " Min(C) = " << gMin(C) << " Max(C) = " << gMax(C) << endl;
    C = max(C,Cmin);
    Info << "limit: CEqn" << " Min(C) = " << gMin(C) << " Max(C) = " << gMax(C) << endl;

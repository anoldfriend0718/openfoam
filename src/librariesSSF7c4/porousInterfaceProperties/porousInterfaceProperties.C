/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     |
    \\  /    A nd           | Copyright (C) 2011-2017 OpenFOAM Foundation
     \\/     M anipulation  |
-------------------------------------------------------------------------------
License
    This file is part of OpenFOAM.

    OpenFOAM is free software: you can redistribute it and/or modify it
    under the terms of the GNU General Public License as published by
    the Free Software Foundation, either version 3 of the License, or
    (at your option) any later version.

    OpenFOAM is distributed in the hope that it will be useful, but WITHOUT
    ANY WARRANTY; without even the implied warranty of MERCHANTABILITY or
    FITNESS FOR A PARTICULAR PURPOSE.  See the GNU General Public License
    for more details.

    You should have received a copy of the GNU General Public License
    along with OpenFOAM.  If not, see <http://www.gnu.org/licenses/>.

\*---------------------------------------------------------------------------*/

// #include "dimensionSets.H"
#include "porousInterfaceProperties.H"
#include "alphaContactAngleFvPatchScalarField.H"
#include "mathematicalConstants.H"
#include "surfaceInterpolate.H"
#include "fvcDiv.H"
#include "fvcGrad.H"
#include "fvcSnGrad.H"
#include "fvcAverage.H"
#include "fvCFD.H"


// * * * * * * * * * * * * * * * Static Member Data  * * * * * * * * * * * * //

const Foam::scalar Foam::interfaceProperties::convertToRad =
    Foam::constant::mathematical::pi/180.0;


// * * * * * * * * * * * * * Private Member Functions  * * * * * * * * * * * //

// Correction for the boundary condition on the unit normal nHat on
// walls to produce the correct contact angle.

// The dynamic contact angle is calculated from the component of the
// velocity on the direction of the interface, parallel to the wall.
void Foam::interfaceProperties::correctContactAngle
(
    surfaceVectorField::Boundary& nHatb,
    const surfaceVectorField::Boundary& gradAlphaf
) const
{
    const fvMesh& mesh = alpha1_.mesh();
    const volScalarField::Boundary& abf = alpha1_.boundaryField();

    const fvBoundaryMesh& boundary = mesh.boundary();

    forAll(boundary, patchi)
    {
        if (isA<alphaContactAngleFvPatchScalarField>(abf[patchi]))
        {
            alphaContactAngleFvPatchScalarField& acap =
                const_cast<alphaContactAngleFvPatchScalarField&>
                (
                    refCast<const alphaContactAngleFvPatchScalarField>
                    (
                        abf[patchi]
                    )
                );

            fvsPatchVectorField& nHatp = nHatb[patchi];
            const scalarField theta
            (
                convertToRad*acap.theta(U_.boundaryField()[patchi], nHatp)
            );

            const vectorField nf
            (
                boundary[patchi].nf()
            );

            // Reset nHatp to correspond to the contact angle

            const scalarField a12(nHatp & nf); 

            const scalarField b1(cos(theta));

            scalarField b2(nHatp.size());
            forAll(b2, facei)
            {
                b2[facei] = cos(acos(a12[facei]) - theta[facei]);
            }

            const scalarField det(1.0 - a12*a12);

            scalarField a((b1 - a12*b2)/det);
            scalarField b((b2 - a12*b1)/det);

            nHatp = a*nf + b*nHatp;
            nHatp /= (mag(nHatp) + deltaN_.value());

            acap.gradient() = (nf & nHatp)*mag(gradAlphaf[patchi]);
            acap.evaluate();
        }
    }
}


// The dynamic contact angle is calculated from the component of the
// velocity on the direction of the interface, parallel to the porosity wall
void Foam::interfaceProperties::correctContactAngleDB
(
    surfaceVectorField& nHatb, 
    surfaceVectorField& gradAlphaf
) const
{
    // Reading Values From Main Simulation
    volScalarField Solid_ = alpha1_.mesh().lookupObject<volScalarField>("Solid"); 

    scalar SmoothingCycles_(transportPropertiesDict_.lookupOrDefault("SmoothingCycles",0));
    scalar SmoothingMethod_(transportPropertiesDict_.lookupOrDefault("SmoothingMethod",1));
    scalar SmootherCount = 0;

    // Defining the Solid Surface Interface
    const surfaceScalarField SolidSurface(fvc::interpolate(Solid_));
    //forAll(SolidSurface, facei) //setting all face B values to 1 at solid interface and zero otherwise
    //        {
    //            if (SolidSurface[facei]>0) //>0 ==1
    //            {SolidSurface[facei]=1;}
    //        }
  
    //Calculating surface gradient at the face and smoothing
    const volVectorField gradSolid(fvc::grad(Solid_, "nHatSolid"));
    surfaceVectorField gradSolidf(fvc::interpolate(gradSolid));
    //surfaceVectorField nHatSolidfv0(gradSolidf/(mag(gradSolidf) + deltaN_));
    surfaceVectorField nHatSolidfv(gradSolidf/(mag(gradSolidf) + deltaN_));

    if (SmoothingMethod_ == 0)
    {
        while (SmootherCount <= SmoothingCycles_)
        {
        if (SmootherCount < SmoothingCycles_)
            {gradSolidf = SolidSurface*fvc::interpolate(fvc::average(gradSolidf));}
        else 
            {gradSolidf = fvc::interpolate(fvc::average(gradSolidf));}
        SmootherCount +=1;  
        }
        nHatSolidfv = gradSolidf/(mag(gradSolidf) + deltaN_);
    }

    //Calculating the normal of the solid surface
    
    if (SmoothingMethod_ == 1)
    {
        //surfaceVectorField nHatSolidfv(gradSolidf/(mag(gradSolidf) + deltaN_));
        surfaceScalarField nHatSolidfvMag0(mag(nHatSolidfv));
        surfaceScalarField nHatSolidfvMag(max(nHatSolidfvMag0,deltaN_.value()));

        // Smoothing the normal of the solid surface
        while (SmootherCount <= SmoothingCycles_)
        {
            //nHatSolidfv = nHatSolidfvMag*fvc::interpolate(fvc::average(nHatSolidfv*nHatSolidfvMag)/(fvc::average(nHatSolidfvMag)));
            nHatSolidfv = nHatSolidfvMag0*fvc::interpolate(fvc::average(nHatSolidfv*nHatSolidfvMag))/fvc::interpolate(fvc::average(nHatSolidfvMag));
            // nHatSolidfv = fvc::interpolate(fvc::average(nHatSolidfv*nHatSolidfvMag)/(fvc::average(nHatSolidfvMag)));
            nHatSolidfv = nHatSolidfv/(mag(nHatSolidfv)+deltaN_.value());
            SmootherCount +=1;  
        }

    }

    if (SmoothingMethod_ == 2)
    {
        volVectorField nHatSolidv(gradSolid/(mag(gradSolid) + deltaN_));
        volScalarField nHatSolidvMag0(mag(nHatSolidv));
        volScalarField nHatSolidvMag(max(nHatSolidvMag0,deltaN_.value()));

        // Smoothing the normal of the solid surface
        while (SmootherCount <= SmoothingCycles_)
        {
            // if (SmootherCount == 0)
            // {   nHatSolidfv = fvc::interpolate(fvc::average(nHatSolidfv*nHatSolidfvMag))/fvc::interpolate(fvc::average(nHatSolidfvMag));  }
            // else
            // { 
                //nHatSolidfv = nHatSolidfvMag*fvc::interpolate(fvc::average(nHatSolidfv*nHatSolidfvMag)/(fvc::average(nHatSolidfvMag)));
                nHatSolidv = nHatSolidvMag0*fvc::average(fvc::interpolate(nHatSolidv*nHatSolidvMag))/fvc::average(fvc::interpolate(nHatSolidvMag));
                // nHatSolidfv = fvc::interpolate(fvc::average(nHatSolidfv*nHatSolidfvMag)/(fvc::average(nHatSolidfvMag)));
                nHatSolidv = nHatSolidv/(mag(nHatSolidv)+deltaN_.value());
            // }
            SmootherCount +=1;  
        }
        //surfaceVectorField nHatSolidfv(fvc::interpolate(nHatSolidv));
        nHatSolidfv = fvc::interpolate(nHatSolidv);
    }

    if (SmoothingMethod_ == 3)
    {
        volVectorField nHatSolidv(gradSolid/(mag(gradSolid) + deltaN_));
        volScalarField nHatSolidvMag0(mag(nHatSolidv));
        volScalarField nHatSolidvMag(max(nHatSolidvMag0,deltaN_.value()));

        // Smoothing the normal of the solid surface
        while (SmootherCount <= SmoothingCycles_)
        {
            // if (SmootherCount == 0)
            // {   nHatSolidfv = fvc::interpolate(fvc::average(nHatSolidfv*nHatSolidfvMag))/fvc::interpolate(fvc::average(nHatSolidfvMag));  }
            // else
            // { 
                //nHatSolidfv = nHatSolidfvMag*fvc::interpolate(fvc::average(nHatSolidfv*nHatSolidfvMag)/(fvc::average(nHatSolidfvMag)));
                nHatSolidv = nHatSolidvMag0*fvc::average(fvc::interpolate(nHatSolidv*nHatSolidvMag))/fvc::average(fvc::interpolate(nHatSolidvMag));
                // nHatSolidfv = fvc::interpolate(fvc::average(nHatSolidfv*nHatSolidfvMag)/(fvc::average(nHatSolidfvMag)));
                nHatSolidv = nHatSolidv/(mag(nHatSolidv)+deltaN_.value());
            // }
            SmootherCount +=1;  
        }
        //surfaceVectorField nHatSolidfv(fvc::interpolate(nHatSolidv));
        nHatSolidfv = fvc::interpolate(nHatSolidv);
        nHatSolidfv /= (mag(nHatSolidfv) + deltaN_.value());
    }

    //Converting Input Contact Angle to Radians
    scalar theta0Rad_ = convertToRad*theta0_;
 
	    //dot product between both normals = cos(thetaI)
            const surfaceScalarField a12(nHatb & nHatSolidfv);
	
	    //cos(theta0)
            const scalar b1(cos(theta0Rad_));

	    //cos(thetaDiff)
            surfaceScalarField b2(cos(acos(a12) - theta0Rad_));

	    //1-cos^2(thetaI)
            const surfaceScalarField det(1.0 - a12*a12);

	    //cos(theta0) - cos(thetaI)*cos(thetaDiff)/ ( 1-cos^2(thetaI) )
            surfaceScalarField a((b1 - a12*b2)/det);

	    //cos(thetaDiff) - cos(thetaI)*cos(theta0)/ ( 1-cos^2(thetaI) )
            surfaceScalarField b((b2 - a12*b1)/det);

	    // correction of the interface and re-normalizing
	    nHatb = a*nHatSolidfv + b*nHatb;
            nHatb /= (mag(nHatb) + deltaN_.value());

        surfaceScalarField nHatbmag(mag(nHatb));
        forAll(nHatbmag,facei)
        {
            if(nHatbmag[facei]<0.5) // if <1 ==0.  
            {	
            nHatbmag[facei]=0; 
            }
            if(nHatbmag[facei]>0.5) // if <1 ==0.  
            {	
            nHatbmag[facei]=1; 
            }
        }
        nHatb *= nHatbmag;
}

/////Calculating the curvature/////   ADD Smoothing by LZY

void Foam::interfaceProperties::calculateK()
{
    const fvMesh& mesh = alpha1_.mesh();
    const surfaceVectorField& Sf = mesh.Sf();
    volScalarField Solid_ = alpha1_.mesh().lookupObject<volScalarField>("Solid"); 
    surfaceScalarField Fluidf_ = alpha1_.mesh().lookupObject<surfaceScalarField>("Fluidf");   // 1216
    surfaceScalarField Fluidf2_ = alpha1_.mesh().lookupObject<surfaceScalarField>("Fluidf2");   // 1216
    scalar CorrGradAlpha_(transportPropertiesDict_.lookupOrDefault("CorrGradAlpha",1));
    scalar cSK2_(transportPropertiesDict_.lookupOrDefault("cSK2",0.9));
    scalar nSK2_(transportPropertiesDict_.lookupOrDefault("nSK2",2));
    //scalar alphaThreshold_(transportPropertiesDict_.lookupOrDefault("alphaThreshold",4.0e-5));
    scalar alphaThresholdK_(transportPropertiesDict_.lookupOrDefault("alphaThresholdK",0.0));
    // scalar gradAlphaThreshold_(transportPropertiesDict_.lookupOrDefault("gradAlphaThreshold",1.0));
        dimensionedScalar small("small", dimless, 1e-10);
    // dimensionedScalar gradAlphaThreshold_("gradAlphaThreshold_",alphaThreshold_/pow(average(mesh.V()), 1.0/3.0));
    volScalarField talpha2(1.0-alpha1_);


        scalar xfcset(transportPropertiesDict_.lookupOrDefault("xfcset",1.0));
    //volScalarField alpha1_ave (alpha1_);

    if (CorrGradAlpha_ == 1)  //alpha_s_ini = 1
    {
        // volScalarField alpha1_ave (fvc::average(fvc::interpolate(max(alpha1_,SMALL),"Fluid_c2f")*Fluidf_)/fvc::average(Fluidf_));
        alpha1_ave = (fvc::average(fvc::interpolate(max(alpha1_,small),"Fluid_c2f")*Fluidf_)/fvc::average(Fluidf_));
        alpha1_ave = (1.0-Solid_)*alpha1_ + Solid_*alpha1_ave;
    }
    if (CorrGradAlpha_ == 2)
    {
        // volScalarField alpha1_ave (fvc::average(fvc::interpolate(max(alpha1_,SMALL),"Fluid_c2f")*Fluidf_)/fvc::average(Fluidf_));
        alpha1_ave = fvc::average(fvc::interpolate(max(alpha1_,SMALL),"Fluid_c2f")*Fluidf_)/fvc::average(Fluidf_);
        alpha1_ave = (1-Solid_)*alpha1_ + Solid_*alpha1_ave;
        alpha1_ave = (fvc::average(fvc::interpolate(max(alpha1_ave,small),"Fluid_c2f")*Fluidf2_)/fvc::average(Fluidf2_));
        alpha1_ave = (1.0-Solid_)*alpha1_ + Solid_*alpha1_ave;
    }

    if (CorrGradAlpha_ == 11)  //alpha_s_ini = 0
    {
        // volScalarField alpha1_ave (fvc::average(fvc::interpolate(max(alpha1_,SMALL),"Fluid_c2f")*Fluidf_)/fvc::average(Fluidf_));
        alpha1_ave = (fvc::average(fvc::interpolate(max(talpha2,small),"Fluid_c2f")*Fluidf_)/fvc::average(Fluidf_));
        alpha1_ave = (1.0-Solid_)*talpha2 + Solid_*alpha1_ave;
        alpha1_ave = 1.0 - alpha1_ave;
    }
    if (CorrGradAlpha_ == 12)
    {
        // volScalarField alpha1_ave (fvc::average(fvc::interpolate(max(alpha1_,SMALL),"Fluid_c2f")*Fluidf_)/fvc::average(Fluidf_));
        alpha1_ave = fvc::average(fvc::interpolate(max(talpha2,SMALL),"Fluid_c2f")*Fluidf_)/fvc::average(Fluidf_);
        alpha1_ave = (1-Solid_)*talpha2 + Solid_*alpha1_ave;
        alpha1_ave = (fvc::average(fvc::interpolate(max(alpha1_ave,small),"Fluid_c2f")*Fluidf2_)/fvc::average(Fluidf2_));
        alpha1_ave = (1.0-Solid_)*talpha2 + Solid_*alpha1_ave;
        alpha1_ave = 1.0 - alpha1_ave;
    }

     
    // else
    // {
    //     volScalarField alpha1_ave (alpha1_);
    // }


    /* -----------------------original hybridInterFoam---------------------------------
    //scalar C_(transportPropertiesDict_.lookupOrDefault("alphaSmoothingCoeff",1.));

    //Smoothing Alpha ()
    volScalarField alpha1Smoothed (C_*alpha1_ +(1-C_)*fvc::average(fvc::interpolate(alpha1_)));

    // Cell gradient of alpha
    //const volVectorField gradAlpha(fvc::grad(alpha1_, "nHat"));
    const volVectorField gradAlpha(fvc::grad(alpha1Smoothed, "nHat")); 

    // Interpolated face-gradient of alpha
    surfaceVectorField gradAlphaf(fvc::interpolate(gradAlpha));

    // Face unit interface normal of the fluid-fluid interface
    surfaceVectorField nHatfv(gradAlphaf/(mag(gradAlphaf) + deltaN_));

    if (activatePorousContactAngle_==1)
    {correctContactAngleDB(nHatfv,gradAlphaf);}
    correctContactAngle(nHatfv.boundaryFieldRef(), gradAlphaf.boundaryField());

    // Face unit interface normal flux
    nHatf_ = nHatfv & Sf; 

    // Simple expression for curvature
    K_ = -fvc::div(nHatf_);
    K_.correctBoundaryConditions();
    -------------------------------------------------------------------------*/

    //-----------------------Smoothing Alpha using Raeini's SSF--------------------------------
    // smooth alpha1_, eq18a-c 
    //scalar CSK(transportPropertiesDict_.lookupOrDefault("curvatureSmoothCoeff",0.5)); 
    //volScalarField alpha1s_ = CSK * (fvc::average(fvc::interpolate(alpha1_))) + (1.0 - CSK) * alpha1_;
    //volScalarField alpha2s_ = CSK * fvc::average(fvc::interpolate(alpha1s_)) + (1.0 - CSK) * alpha1s_;
    //volScalarField alpha3s_ = CSK * fvc::average(fvc::interpolate(alpha2s_)) + (1.0 - CSK) * alpha2s_;

    // indicator of Fluid region
    //volScalarField Fluid_ (1.0-Solid1_); 
    //surfaceScalarField Fluidf_ (fvc::interpolate(Fluid_,"Fluid_c2f")); 



    // init alpha smoothed    from interGCFoam
    // volScalarField alpha1s = alpha1_ave;
    volScalarField alpha1s ( min(max((1.0+2.0*alphaThresholdK_)*(alpha1_ave-0.5)+0.5,0.0),1.0) );  // clip //??
    // forAll(alpha1s,celli) 
    //     {
    //         if(alpha1s[celli] <= alphaThreshold_)   //will cause non-physics fc
    //         {	
    //             alpha1s[celli]=0; 
    //         }
    //         if(alpha1s[celli] >= (1-alphaThreshold_)) 
    //         {	
    //             alpha1s[celli]=1; 
    //         }
    //     }
    volScalarField alpha1s_0 = alpha1s;    
    for (int i=0;i<nSK_;i++)
    {
        alpha1s = Solid_*alpha1_ave + (1.0-Solid_)*(cSK_ * fvc::average(fvc::interpolate(alpha1s*(1.0-Solid_)))/fvc::average(fvc::interpolate(max(1.0-Solid_,small))) + (1.0 - cSK_) * alpha1s);
    }

    for (int i=0;i<nSK2_;i++)
    {
        alpha1s_0 = Solid_*alpha1_ave + (1.0-Solid_)*(cSK2_ * fvc::average(fvc::interpolate(alpha1s_0*(1.0-Solid_)))/fvc::average(fvc::interpolate(max(1.0-Solid_,small))) + (1.0 - cSK2_) * alpha1s_0);
    }
    
        //Calculating the normal of the solid surface
    // surfaceVectorField nHatSolidfv(gradSolidf/(mag(gradSolidf) + deltaN_));
            //const volVectorField gradAlpha(fvc::grad(alpha1_, "nHat"));
     const volVectorField gradAlpha(fvc::grad(alpha1s, "nHat")); 
     const volVectorField gradAlpha_0(fvc::grad(alpha1s_0, "nHat"));
    // volVectorField gradAlpha(fvc::grad(alpha1s, "nHat")); 
    // volVectorField gradAlpha_0(fvc::grad(alpha1s_0, "nHat"));
    // volScalarField gradAlpha_mag(mag(gradAlpha));  
    // volScalarField gradAlpha_mag0(mag(gradAlpha_0));
    // volScalarField gradAlpha_nfilt(0.0*alpha1_);  
    // volScalarField gradAlpha_nfilt0(0.0*alpha1_);
    // forAll(gradAlpha_mag,celli) 
    //     {
    //         if(gradAlpha_mag[celli] <= gradAlphaThreshold_)   //will cause non-physics fc
    //         {	
    //             gradAlpha_nfilt[celli]=0.0; 
    //         }
    //         else
    //         {
    //             gradAlpha_nfilt[celli]=1.0; 
    //         }
    //         if(gradAlpha_mag0[celli] <= gradAlphaThreshold_) 
    //         {	
    //             gradAlpha_nfilt0[celli]=0.0; 
    //         }
    //         else
    //         {
    //             gradAlpha_nfilt0[celli]=1.0;
    //         }
        

    //     }
    // gradAlpha *= gradAlpha_nfilt;
    // gradAlpha_0 *= gradAlpha_nfilt0;

    
        // const volVectorField gradAlpha((1.0-Solid_)*fvc::grad(alpha1s, "nHat")); 
    //volVectorField gradAlphas (fvc::average(fvc::interpolate(gradAlpha*Fluid_)));
    //update interface vector at cell center, used for filtering
    nI_ = gradAlpha/(Foam::mag(gradAlpha) + deltaN_);
    
    //nI_ = gradAlpha/(Foam::mag(gradAlpha) + deltaN_)*Fluid_;
    // nI_ = fvc::average(fvc::interpolate(nI_));   // unit normal can not interpolate, mag may less than 1 !!!

    // Interpolated face-gradient of alpha
    //surfaceVectorField gradAlphaf(fvc::interpolate(gradAlpha)*Fluidf_);  // should be fvc::fvc::interpolate(gradAlpha)*fvc::interpolate(solidv,"solidvf")
    surfaceVectorField gradAlphaf(fvc::interpolate(gradAlpha)); 
    // Face unit interface normal of the fluid-fluid interface
    //surfaceVectorField nHatfv_g(gradAlphaf/(mag(gradAlphaf) + deltaN_));    //eq 19
    //surfaceVectorField nHatfv(fvc::interpolate(nI_)*Fluidf_);  // should be fvc::fvc::interpolate(nI_)*fvc::interpolate(solidv,"solidvf")
    // surfaceVectorField nHatfv(fvc::interpolate(nI_)); 
    //volVectorField nI_s (nI_);
    volVectorField nI_0 (gradAlpha_0/(Foam::mag(gradAlpha_0) + deltaN_));
    // nHatfv /= (mag(nHatfv)+deltaN_.value());

    
    // if (nSN_ > 0)  // smooth normal of the interface
    // {
    //     //volVectorField nI_s (nI_);
    //     volScalarField nI_s_Mag0(mag(nI_s));
    //     volScalarField nI_s_Mag(max(nI_s_Mag0,deltaN_.value()));

    //     // Smoothing the normal of the interface  [default nsK 1 cSK 1]
    //     for (int i=0;i<nSN_;i++)
    //     {
    //         //nI_ = mag(nI0_)*fvc::average((fvc::interpolate(nI_*nI0_Mag))/(fvc::interpolate(nI0_Mag)));
    //         nI_s = nI_s_Mag0*fvc::average(fvc::interpolate(nI_s*nI_s_Mag))/fvc::average(fvc::interpolate(nI_s_Mag));
    //         nI_s = nI_s/(mag(nI_s) + deltaN_.value());
    //     }
    // }

    surfaceVectorField nHatfv(fvc::interpolate(nI_));   // the key to reduce Usc. BUT WHY?
    // Cell gradient of alpha3s_
    nIsi_ = fvc::average(nHatfv);
    surfaceVectorField nHatfv_0 (fvc::interpolate(nI_0)*Fluidf2_); 
    surfaceVectorField nHatfv_g (gradAlphaf/(mag(gradAlphaf) + deltaN_)); 
    //volVectorField ns = fvc::grad(Solid_)/(mag(fvc::grad(Solid_)) + deltaN_);   //test
    //volVectorField ns = fvc::grad(Solid_)/(mag(fvc::grad(Solid_)) + dimensionedScalar(dimless/dimLength, small));  //test 
    //volVectorField gradAlphac(gradAlpha - ns*mag(gradAlpha));    //test
    //surfaceVectorField gradAlphaf(fvc::interpolate(gradAlphac));      //test

    surfaceVectorField nHatfv_ave (nHatfv*xfcset+nHatfv_g*(1.0-xfcset));

    if (activatePorousContactAngle_==1)
    {
        correctContactAngleDB(nHatfv_ave,gradAlphaf);
        // correctContactAngleDB(nHatfv_g,gradAlphaf);
        }
    correctContactAngle(nHatfv_ave.boundaryFieldRef(), gradAlphaf.boundaryField());
    //correctContactAngle(nHatfv_g.boundaryFieldRef(), gradAlphaf.boundaryField());

    // Face unit interface normal flux
    nHatf_ = nHatfv_ave & Sf; 
    
    // Simple expression for curvature
    // volScalarField K1_ = -fvc::div(nHatf_);    // eq 20

    // Complex expression for curvature.
    // Correction is formally zero but numerically non-zero.
    //  Kcorr=0 simple, Kcorr=1 Complex (default simple)
    K_ = -fvc::div(nHatf_) + (nI_ & fvc::grad(nHatfv) & nI_)*Kcorr_;

    K_.correctBoundaryConditions();

    // nHatf_0 = nHatfv_0 & Sf;
    // K_0 = -fvc::div(nHatf_0);   //for calculate theta K<0, theta = 180-theta

    // nHatf_g = nHatfv_g & Sf;
    // K_g = -fvc::div(nHatf_g);   //for calculate theta K<0, theta = 180-theta

    if (nSC_>0)
    {
        // apha_c
        volScalarField alpha1c_ = Foam::min(1.0, Foam::max(alpha1_ave, 0.0));  //alpha1_ave  //alpha1_
        // eq 21b(w)  
        volScalarField w = Foam::sqrt(alpha1c_*(1.0 - alpha1c_) + 1e-6);
        volScalarField factor = 2.0 * Foam::sqrt(alpha1c_*(1.0 - alpha1c_));
        volScalarField Kstar = fvc::average(fvc::interpolate(K_*w*(1.0-Solid_)))/fvc::average(fvc::interpolate(w*max(1.0-Solid_,small)));
        volScalarField Ks = factor * K_ + (1.0 - factor) * Kstar;

        for (int i=1;i<nSC_;i++)
        {
            Kstar = fvc::average(fvc::interpolate(Ks*w*(1.0-Solid_)))/fvc::average(fvc::interpolate(w*max(1.0-Solid_,small)));
            Ks = factor * K_ + (1.0 - factor) * Kstar;
        }

        Kf_ = fvc::interpolate(w*Ks*(1.0-Solid_))/fvc::interpolate(w*max(1.0-Solid_,small));
    }
    // if (CorrGradAlpha_ > 0)
    //     alpha_pc = 1.0/(1.0-cPc_)*(min( max(alpha1_ave,cPc_/2.0), (1.0-cPc_/2.0) ) - cPc_/2.0);
    // else
    //     alpha_pc = 1.0/(1.0-cPc_)*(min( max(alpha1_,cPc_/2.0), (1.0-cPc_/2.0) ) - cPc_/2.0);

    // forAll(alpha_pc,celli) 
    //     {
    //         if(alpha_pc[celli] <= alphaThreshold_)   //will cause non-physics fc
    //         {	
    //             alpha_pc[celli]=0.0; 
    //         }
    //         else if(alpha_pc[celli] >= (1.0-alphaThreshold_)) 
    //         {
    //             alpha_pc[celli]=1.0; 
    //         }
    //     }


    // alpha_pc.correctBoundaryConditions();

}


void Foam::interfaceProperties::calculateFc() //const
{
    const fvMesh& mesh = alpha1_.mesh();
    const surfaceScalarField& magSf = mesh.magSf();
    surfaceScalarField Solidf0_ = alpha1_.mesh().lookupObject<surfaceScalarField>("Solidf0");   // 1216

    scalar CorrGradAlpha_(transportPropertiesDict_.lookupOrDefault("CorrGradAlpha",1));
    scalar alphaThresholdFc_(transportPropertiesDict_.lookupOrDefault("alphaThresholdFc",0.001));
    // scalar gradAlphaThreshold_(transportPropertiesDict_.lookupOrDefault("gradAlphaThreshold",1.0));  //gradAlphaThreshold
    scalar SoomthAlphapc_(transportPropertiesDict_.lookupOrDefault("SoomthAlphapc",0.0));
    // Info << "gradAlphaThreshold = " << gradAlphaThreshold_ << endl;
    // scalar thetaLimit_(transportPropertiesDict_.lookupOrDefault("thetaLimit",80));
    // scalar alphalsLimit_(transportPropertiesDict_.lookupOrDefault("alphalsLimit",0.8));
    // //scalar xfcset(transportPropertiesDict_.lookupOrDefault("xfcset",1.0));
    // scalar thetaLimit2_(transportPropertiesDict_.lookupOrDefault("thetaLimit2",90));
        //volScalarField alpha1_ave (alpha1_);

    // if (CorrGradAlpha_ == 1)
    // {
    //     // volScalarField alpha1_ave (fvc::average(fvc::interpolate(max(alpha1_,SMALL),"Fluid_c2f")*Fluidf_)/fvc::average(Fluidf_));
    //     alpha1_ave = (fvc::average(fvc::interpolate(max(alpha1_,small),"Fluid_c2f")*Fluidf_)/fvc::average(Fluidf_));
    //     alpha1_ave = (1.0-Solid_)*alpha1_ + Solid_*alpha1_ave;
    // }

    // volScalarField alpha_pc
    // (
    //     IOobject 
    //     (
    //         "alpha_pc",
    //         alpha1_.time().timeName(),
    //         alpha1_.mesh()
    //     ),
    //  alpha1_
    //  //pc_.boundaryField().types()
    // );

    // Sharpen interface function
    // Raeini's thesis, equation (2.22) and (2.23)
    // if (CorrGradAlpha_ > 0)
    //     alpha_pc = 1.0/(1.0-cPc_)*(min( max(alpha1_ave,cPc_/2.0), (1.0-cPc_/2.0) ) - cPc_/2.0);
    // else
    //     alpha_pc = 1.0/(1.0-cPc_)*(min( max(alpha1_,cPc_/2.0), (1.0-cPc_/2.0) ) - cPc_/2.0);

    // forAll(alpha_pc,celli) 
    //     {
    //         if(alpha_pc[celli] <= alphaThreshold_)   //will cause non-physics fc
    //         {	
    //             alpha_pc[celli]=0.0; 
    //         }
    //         else if(alpha_pc[celli] >= (1.0-alphaThreshold_)) 
    //         {
    //             alpha_pc[celli]=1.0; 
    //         }
    //     }


    // alpha_pc.correctBoundaryConditions();

    
    if (CorrGradAlpha_ > 0)
    {
        alpha_pc = min(max((1.0+2.0*alphaThresholdFc_)*(alpha1_ave-0.5)+0.5,0.0),1.0);  // clip
        alpha_pc = 1.0/(1.0-cPc_)*(min( max(alpha_pc,cPc_/2.0), (1.0-cPc_/2.0) ) - cPc_/2.0);
        //alpha_pc = alpha1_ave;
    }
    else
    {
        //alpha_pc = 1.0/(1.0-cPc_)*(min( max(alpha1_,cPc_/2.0), (1.0-cPc_/2.0) ) - cPc_/2.0);
        alpha_pc = min(max((1.0+2.0*alphaThresholdFc_)*(alpha1_-0.5)+0.5,0.0),1.0);  // clip
        alpha_pc = 1.0/(1.0-cPc_)*(min( max(alpha_pc,cPc_/2.0), (1.0-cPc_/2.0) ) - cPc_/2.0);
        //alpha_pc = 1.0*alpha1_;
    }

    if (SoomthAlphapc_ > 0.0001)  // smooth after sharp
    {
        alpha_pc = SoomthAlphapc_*fvc::average(alpha_pc) + (1.0-SoomthAlphapc_)*alpha_pc;
    }

    alpha_pc.correctBoundaryConditions();

    deltasf0_ = fvc::snGrad(alpha_pc);

    // volScalarField magGradAlpha (mag(fvc::Grad(alpha_pc)));

    // Info << "max deltasf0_ = " << Foam::max(deltasf0_) << endl;

    // volScalarField magGradAlpha (mag(fvc::grad(alpha_pc)));
    // volScalarField filterGradApha (alpha_pc);

    // forAll(magGradAlpha,celli) 
    // {
    //     if(magGradAlpha[celli] <= gradAlphaThreshold_)   //will cause non-physics fc
    //     {	
    //         filterGradApha[celli]=1.0;
    //     }
    //     else
    //     {
    //         filterGradApha[celli]=0.0;
    //     }
    // }
    // surfaceScalarField filterGradAphaf (fvc::interpolate(filterGradApha));

    //surface tension force
    stf0_ = sigmaKSSF()*deltasf0_;
    // surfaceScalarField stf_g = fvc::interpolate( (sigmaPtr_->sigma())*K_g)*deltasf;
    // fcg_ = fvc::reconstruct(stf_g*Solidf0_*mesh.magSf());
    
    //volVectorField nIsi (fvc::average(nHatfv));
    fc0_ = fvc::reconstruct(stf0_*Solidf0_*mesh.magSf());

    // forAll(deltasf,facei) 
    //     {
    //         if(mag(deltasf[facei]) <= gradAlphaThreshold_)   //will cause non-physics fc
    //         {	
    //             deltasf[facei]=0.0;
    //             stf[facei]=0.0;
    //         }
    //     }
    stf_ =  stf0_;
    deltasf_ = deltasf0_;
    fc_ = fc0_;
    //stf = sigmaKSSF()*deltasf;

    // forAll(stf_,facei) 
    // {
    //     if(filterGradAphaf[facei] > 0.1)   //will cause non-physics fc
    //     {	
    //             // deltasf[facei]=0.0;
    //         stf_[facei]=0.0;
    //         deltasf_[facei]=0.0;
    //     }
    // }
    // // surfaceScalarField stf_g = fvc::interpolate( (sigmaPtr_->sigma())*K_g)*deltasf;
    // // fcg_ = fvc::reconstruct(stf_g*Solidf0_*mesh.magSf());
    
    // //volVectorField nIsi (fvc::average(nHatfv));
    // fc_ = fvc::reconstruct(stf_*Solidf0_*mesh.magSf());    



    // dimensionedScalar dimfc ("dimfc",fc0_.dimensions(),1);
    // volVectorField nfc2 (fc0_/(mag(fc0_) + dimfc*deltaN_.value()));
    // volScalarField afc2 (nfc2 & nIsi_);
    // // volScalarField theta_fc2 (acos(afc2)/convertToRad);
    // theta_fc1 = (acos(afc2)/convertToRad);   //before correction
    // theta_fc2 = theta_fc1;

    // volVectorField nfc1 (fc0g_/(mag(fc0g_) + dimfc*deltaN_.value()));
    // volScalarField afc1 (nfc1 & nIsi_);
    // // volScalarField theta_fc2 (acos(afc2)/convertToRad);
    // theta_fc1 = (acos(afc1)/convertToRad);
    // volScalarField xfc (0.0*alpha1_);
    // alpha_ls = mag((alpha1_ave-0.5)*2.0);
    // forAll(theta_fc1,celli)
    // {
    //     if (K_0[celli]<0)
    //     {
    //         theta_fc2[celli] = 180.0 - theta_fc1[celli];
    //     }
    //     else
    //     {
    //         theta_fc2[celli] = theta_fc1[celli];
    //     }
    //     if(theta_fc2[celli]>thetaLimit_ && alpha_ls[celli]>alphalsLimit_) // if <1 ==0.   60-180 
    //     {	
    //        xfc[celli] = 0.0; 
    //     }
    //     else
    //     {
    //        xfc[celli] = 1.0; 
    //     }
    // }
    // fc_ = xfc*fc0_;
    // dimensionedScalar theta_sum (Foam::sum(theta_fc2-90.0));  // <0 nothing, >0 180-theta_fc
    // if (theta_sum.value()>0)
    // {
    //     theta_fc2 = 180.0 - theta_fc2;
    //     theta_fc1 = 180.0 - theta_fc1;
    // }

    // volScalarField xfc2(0.0*alpha1_);//volScalarField afc1 (nfc1 & nI1b);
    // volScalarField xfc1(0.0*alpha1_);

    // forAll(theta_fc2,celli)
    // {
    //     if (K_0[celli]<0)
    //     {
    //         theta_fc2[celli] = 180.0 - theta_fc1[celli];
    //     }
    //     else
    //     {
    //         theta_fc2[celli] = theta_fc1[celli];
    //     }

    //     if(theta_fc2[celli]>thetaLimit2_ && theta_fc1[celli]>thetaLimit2_) // if <1 ==0.   60-180 
    //     {	
    //        xfc2[celli] = 0.0; 
    //        xfc1[celli] = 0.0; 
    //     }
    //     else if(theta_fc2[celli]>thetaLimit_) // if <1 ==0.   20-60
    //     {	
    //        if (theta_fc2[celli] <= theta_fc1[celli])
    //        {
    //            xfc2[celli] = xfcset;   //xfc
    //            xfc1[celli] = 1.0-xfcset; 
    //        }
    //        else
    //        {
    //            xfc2[celli] = 1.0-xfcset; 
    //            xfc1[celli] = xfcset; 
    //        }
    //     }
    //     else    // 0-20
    //     {
    //        xfc2[celli] = xfcset; 
    //        xfc1[celli] = 1.0-xfcset; 
    //     }
    // }
    // fc_ = xfc2*fc0_ + xfc1*fcg_;
    // fc_ = fc0_;

    // stfc_ = (fvc::interpolate(fc_) & mesh.Sf())/mesh.magSf();
    // volScalarField xfc(alpha1_);//volScalarField afc1 (nfc1 & nI1b);
    // forAll(theta_fc2,celli)
    // {
    //     if(theta_fc2[celli]>thetaLimit_) // if <1 ==0.  
    //     {	
    //        xfc[celli] = 0.0; 
    //     }
    //     else
    //     {
    //        xfc[celli] = 1.0; 
    //     }
    // }
    // // volVectorField fc_ (xfc*fc0_);
    // fc_ = xfc*fc0_;

}
 
/*
    // apha_c
    volScalarField alpha1c_ = Foam::min(1.0, Foam::max(alpha1_, 0.0));
    const dictionary& MULEScontrols = alpha1_.mesh().solverDict(alpha1_.name());
    scalar maxUnboundedness
    (
        readScalar(MULEScontrols.lookup("maxUnboundedness"))
    );
// eq 21b(w)  
    volScalarField w = Foam::sqrt(alpha1c_*(1.0 - alpha1c_) + maxUnboundedness);
    //volScalarField factor = 2.0*Foam::sqrt(alpha1_*(1.0 - alpha1_)+1e-6) 
    //  * pos(alpha1_-1e-6) * pos(0.9999999-alpha1_);
    volScalarField factor = 2.0 * Foam::sqrt(alpha1c_*(1.0 - alpha1c_));
    //volScalarField factor = 2.0 * w; // alternative to clipping alpha1 - use w
    // obtain smoothed curvature, eq. 15 & 16
// eq 21b(Ks_star)  
    volScalarField Ks_star_ = 
        fvc::average(fvc::interpolate(K1_*w))/fvc::average(fvc::interpolate(w));
// eq 21a
    volScalarField Ks1_ = factor * K1_ + (1.0 - factor) * Ks_star_;
// eq 21c
    Ks_star_ = 
        fvc::average(fvc::interpolate(Ks1_*w))/fvc::average(fvc::interpolate(w));

    volScalarField Ks2_ = factor * K1_ + (1.0 - factor) * Ks_star_;

// eq 21d
    K_ = fvc::interpolate(w*Ks2_)/fvc::interpolate(w);

*/
    //K_ = -fvc::div(nHatf_);              
    //K_.correctBoundaryConditions();




// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

Foam::interfaceProperties::interfaceProperties
(
    const volScalarField& alpha1,
    const volVectorField& U,
    const IOdictionary& dict
)
:
    transportPropertiesDict_(dict),
    cAlpha_
    (
        readScalar
        (
            alpha1.mesh().solverDict(alpha1.name()).lookup("cAlpha")
        )
    ),
    //- smoothing coefficient for alphaS (generally 0.5)
    cSK_
    (
        readScalar
        (
            alpha1.mesh().solverDict(alpha1.name()).lookup("cSK")
        )
    ),
    //-number of smoothing cycle for alphaS
    nSK_
    (
        readScalar
        (
            alpha1.mesh().solverDict(alpha1.name()).lookup("nSK")
        )
    ),

    //-number of smoothing cycle for curvature K_
    nSC_
    (
        readScalar
        (
            alpha1.mesh().solverDict(alpha1.name()).lookup("nSC")
        )
    ),

    //-number of smoothing cycle for curvature K_
    Kcorr_
    (
        readScalar
        (
            alpha1.mesh().solverDict(alpha1.name()).lookup("Kcorr")
        )
    ),

    //- Sharp force coefficient (put 0.98-0.99 for static problems, 0.4-0.5 for dynamic)
    cPc_
    (
       readScalar
       (
           alpha1.mesh().solverDict(alpha1.name()).lookup("cPc")
       )
    ),


    theta0_ 
    (
           transportPropertiesDict_.lookupOrDefault("theta0",90)
    ),

    activatePorousContactAngle_
    (
           transportPropertiesDict_.lookupOrDefault("activatePorousContactAngle",0)
    ),

    sigmaPtr_(surfaceTensionModel::New(dict, alpha1.mesh())),

    deltaN_
    (
        "deltaN",
        1e-8/pow(average(alpha1.mesh().V()), 1.0/3.0)
    ), 

    alpha1_(alpha1),
    alpha1_ave
    (
        IOobject
        (
            "alpha1_ave",
            alpha1_.time().timeName(),
            alpha1_.mesh(),
            IOobject::NO_READ,
            IOobject::AUTO_WRITE
        ),
        alpha1_
    ),
    alpha_pc
    (
        IOobject
        (
            "alpha_pc",
            alpha1_.time().timeName(),
            alpha1_.mesh(),
            IOobject::NO_READ,
            IOobject::AUTO_WRITE
        ),
        alpha1_
    ),
    alpha_ls
    (
        IOobject
        (
            "alpha_ls",
            alpha1_.time().timeName(),
            alpha1_.mesh(),
            IOobject::NO_READ,
            IOobject::AUTO_WRITE
        ),
        alpha1_ave
    ),
    

    U_(U),

    nHatf_
    (
        IOobject
        (
            "nHatf",
            alpha1_.time().timeName(),
            alpha1_.mesh()
        ),
        alpha1_.mesh(),
        dimensionedScalar("nHatf", dimArea, 0.0)
    ),

    nHatf_0
    (
        IOobject
        (
            "nHatf0",
            alpha1_.time().timeName(),
            alpha1_.mesh()
        ),
        alpha1_.mesh(),
        dimensionedScalar("nHatf0", dimArea, 0.0)
    ),

    nHatf_g
    (
        IOobject
        (
            "nHatf_g",
            alpha1_.time().timeName(),
            alpha1_.mesh()
        ),
        alpha1_.mesh(),
        dimensionedScalar("nHatf_g", dimArea, 0.0)
    ),

    //interface normal vector at cell center
    nI_
    (
        IOobject
        (
            "nI",
            alpha1_.time().timeName(),
            alpha1_.mesh(),
            IOobject::NO_READ,
            IOobject::NO_WRITE
        ),
        alpha1_.mesh(),
        dimensionedVector(dimless, vector(0.0,0.0,0.0)) //(0.0 0.0 0.0)
    ),
        nI_0
    (
        IOobject
        (
            "nI0",
            alpha1_.time().timeName(),
            alpha1_.mesh(),
            IOobject::NO_READ,
            IOobject::NO_WRITE
        ),
        alpha1_.mesh(),
        dimensionedVector(dimless, vector(0.0,0.0,0.0)) //(0.0 0.0 0.0)
    ),

        //interface normal vector at cell center
    nIsi_
    (
        IOobject
        (
            "nIsi",
            alpha1_.time().timeName(),
            alpha1_.mesh(),
            IOobject::NO_READ,
            IOobject::NO_WRITE
        ),
        alpha1_.mesh(),
        dimensionedVector(dimless, vector(0.0,0.0,0.0)) //(0.0 0.0 0.0)
    ),

    K_
    (
        IOobject
        (
            "interfaceProperties:K",
            alpha1_.time().timeName(),
            alpha1_.mesh()
        ),
        alpha1_.mesh(),
        dimensionedScalar("K", dimless/dimLength, 0.0)
    ),
        K_0
    (
        IOobject
        (
            "interfaceProperties:K0",
            alpha1_.time().timeName(),
            alpha1_.mesh()
        ),
        alpha1_.mesh(),
        dimensionedScalar("K0", dimless/dimLength, 0.0)
    ),
        K_g
    (
        IOobject
        (
            "interfaceProperties:Kg",
            alpha1_.time().timeName(),
            alpha1_.mesh()
        ),
        alpha1_.mesh(),
        dimensionedScalar("Kg", dimless/dimLength, 0.0)
    ),
    deltasf0_
    (
        IOobject
        (
            "deltasf0",
            alpha1_.time().timeName(),
            alpha1_.mesh()
        ),
        alpha1_.mesh(),
        dimensionedScalar(dimless/dimLength, 0.0)
    ),
    deltasf_
    (
        IOobject
        (
            "deltasf",
            alpha1_.time().timeName(),
            alpha1_.mesh()
        ),
        alpha1_.mesh(),
        dimensionedScalar( dimless/dimLength, 0.0)
    ),
    stf0_
    (
        IOobject
        (
            "stf0",
            alpha1_.time().timeName(),
            alpha1_.mesh()
        ),
        alpha1_.mesh(),
        dimensionedScalar(dimMass/(dimLength*dimLength*dimTime*dimTime), 0.0)
    ),
    stf_
    (
        IOobject
        (
            "stf",
            alpha1_.time().timeName(),
            alpha1_.mesh()
        ),
        alpha1_.mesh(),
        dimensionedScalar( dimMass/(dimLength*dimLength*dimTime*dimTime), 0.0)
    ),


    Kf_
    (
        IOobject
        (
            "interfaceProperties:Kf",
            alpha1_.time().timeName(),
            alpha1_.mesh()
        ),
        alpha1_.mesh(),
        dimensionedScalar("Kf", dimless/dimLength, 0.0)
    ),

    // stfc_
    // (
    //     IOobject
    //     (
    //         "stfc",
    //         alpha1_.time().timeName(),
    //         alpha1_.mesh()
    //     ),
    //     alpha1_.mesh(),
    //     dimensionedScalar("stfc", Foam::dimPressure/dimLength, 0.0)
    // ),

    // theta_fc2
   theta_fc2
    (
        IOobject
        (
            "theta_fc2",
            alpha1_.time().timeName(),
            alpha1_.mesh()
        ),
        alpha1_.mesh(),
        dimensionedScalar(dimless, 0.0)
    ),

    // theta_fc2
   theta_fc1
    (
        IOobject
        (
            "theta_fc1",
            alpha1_.time().timeName(),
            alpha1_.mesh()
        ),
        alpha1_.mesh(),
        dimensionedScalar(dimless, 0.0)
    ),

    fc_
    (
        IOobject
        (
            "fc_",
            alpha1_.time().timeName(),
            alpha1_.mesh(),
            IOobject::NO_READ,
            IOobject::NO_WRITE
        ),
        alpha1_.mesh(),
        dimensionedVector(Foam::dimPressure/dimLength, vector(0.0,0.0,0.0)) //(0.0 0.0 0.0)
    ),
    //     fc0g_
    // (
    //     IOobject
    //     (
    //         "fc0g_",
    //         alpha1_.time().timeName(),
    //         alpha1_.mesh(),
    //         IOobject::NO_READ,
    //         IOobject::NO_WRITE
    //     ),
    //     alpha1_.mesh(),
    //     dimensionedVector(Foam::dimPressure/dimLength, vector(0.0,0.0,0.0)) //(0.0 0.0 0.0)
    // ),
    fc0_
    (
        IOobject
        (
            "fc0_",
            alpha1_.time().timeName(),
            alpha1_.mesh(),
            IOobject::NO_READ,
            IOobject::NO_WRITE
        ),
        alpha1_.mesh(),
        dimensionedVector(Foam::dimPressure/dimLength, vector(0.0,0.0,0.0)) //(0.0 0.0 0.0)
    ),
        fcg_
    (
        IOobject
        (
            "fcg_",
            alpha1_.time().timeName(),
            alpha1_.mesh(),
            IOobject::NO_READ,
            IOobject::NO_WRITE
        ),
        alpha1_.mesh(),
        dimensionedVector(Foam::dimPressure/dimLength, vector(0.0,0.0,0.0)) //(0.0 0.0 0.0)
    )
    //fvc::reconstruct(fcf*Solidf0*mesh.magSf())

{
    calculateK();
    //calculateFc();
}


// * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * * //

//Foam::tmp<Foam::volScalarField>
//Foam::interfaceProperties::sigmaK() const
//{
//    return sigmaPtr_->sigma()*K_;
//}

//Foam::tmp<Foam::surfaceScalarField>
//Foam::interfaceProperties::surfaceTensionForce() const
//{
//   return fvc::interpolate(sigmaK())*fvc::snGrad(alpha1_); 
//}

Foam::tmp<Foam::surfaceScalarField>
Foam::interfaceProperties::sigmaKSSF() const
{
    if (nSC_>0)
        {return fvc::interpolate( sigmaPtr_->sigma() )*Kf_;}    // smooth K_
    else
        {return fvc::interpolate( (sigmaPtr_->sigma())*K_);}  //don't smooth K_ 
}

Foam::tmp<Foam::surfaceScalarField>
Foam::interfaceProperties::surfaceTensionForce() const
{
   return stf0_; // sigmaKSSF()*fvc::snGrad(alpha_pc); 
   //return stfc_;
   //return fvc::interpolate(sigmaK())*fvc::snGrad(alpha_pc); 
   // stf0_ = sigmaKSSF()*deltasf0_;
}

Foam::tmp<Foam::surfaceScalarField>
Foam::interfaceProperties::stf() const
{
   return stf_; 
   //return stfc_;
   //return fvc::interpolate(sigmaK())*fvc::snGrad(alpha_pc); 
}
Foam::tmp<Foam::surfaceScalarField>
Foam::interfaceProperties::stf0() const
{
   return stf0_; 
   //return stfc_;
   //return fvc::interpolate(sigmaK())*fvc::snGrad(alpha_pc); 
}
Foam::tmp<Foam::surfaceScalarField>
Foam::interfaceProperties::deltasf() const
{
   return deltasf_; 
   //return stfc_;
   //return fvc::interpolate(sigmaK())*fvc::snGrad(alpha_pc); 
}
Foam::tmp<Foam::surfaceScalarField>
Foam::interfaceProperties::deltasf0() const
{
   return deltasf0_; 
   //return stfc_;
   //return fvc::interpolate(sigmaK())*fvc::snGrad(alpha_pc); 
}

Foam::tmp<Foam::volScalarField>
Foam::interfaceProperties::Klg() const
{
    return K_;
}

Foam::tmp<Foam::volVectorField>
Foam::interfaceProperties::nlg() const
{
    return nI_;
}

Foam::tmp<Foam::volVectorField>
Foam::interfaceProperties::fc() const
{
    return fc_;
}

Foam::tmp<Foam::volVectorField>
Foam::interfaceProperties::fc0() const
{
    return fc0_;
}

Foam::tmp<Foam::volScalarField>
Foam::interfaceProperties::thetafc() const
{
    return theta_fc2;
}

Foam::tmp<Foam::volScalarField>
Foam::interfaceProperties::thetafc0() const
{
    return theta_fc1;
}

Foam::tmp<Foam::volScalarField>
Foam::interfaceProperties::alphals() const
{
    return alpha_ls;
}

Foam::tmp<Foam::volScalarField>
Foam::interfaceProperties::alphapc() const
{
    return alpha_pc;
}

Foam::tmp<Foam::volScalarField>
Foam::interfaceProperties::nearInterface() const
{
    return pos(alpha1_ - 0.01)*pos(0.99 - alpha1_);
}


void Foam::interfaceProperties::correct()
{
    calculateK();
    calculateFc();
}


bool Foam::interfaceProperties::read()
{
    alpha1_.mesh().solverDict(alpha1_.name()).lookup("cAlpha") >> cAlpha_;
    sigmaPtr_->readDict(transportPropertiesDict_);

    return true;
}


// ************************************************************************* //

## Limit
```cpp
template
<
    class RdeltaTType,
    class RhoType,
    class SpType,
    class SuType,
    class PsiMaxType,
    class PsiMinType
>
void Foam::MULES::limiterCorr
(
    scalarField& allLambda,
    const RdeltaTType& rDeltaT,
    const RhoType& rho,
    const volScalarField& psi,
    const surfaceScalarField& phi,
    const surfaceScalarField& phiCorr,
    const SpType& Sp,
    const SuType& Su,
    const PsiMaxType& psiMax,
    const PsiMinType& psiMin
)
{
    const scalarField& psiIf = psi;
    const volScalarField::Boundary& psiBf = psi.boundaryField();

    const fvMesh& mesh = psi.mesh();

    const dictionary& MULEScontrols = mesh.solverDict(psi.name());

    const label nLimiterIter
    (
        readLabel(MULEScontrols.lookup("nLimiterIter"))
    );

    const scalar smoothLimiter
    (
        MULEScontrols.lookupOrDefault<scalar>("smoothLimiter", 0)
    );

    const scalar extremaCoeff
    (
        MULEScontrols.lookupOrDefault<scalar>("extremaCoeff", 0)
    );

    const scalar boundaryExtremaCoeff
    (
        MULEScontrols.lookupOrDefault<scalar>
        (
            "boundaryExtremaCoeff",
            extremaCoeff
        )
    );

    const scalar boundaryDeltaExtremaCoeff
    (
        max(boundaryExtremaCoeff - extremaCoeff, 0)
    );

    const labelUList& owner = mesh.owner();
    const labelUList& neighb = mesh.neighbour();
    tmp<volScalarField::Internal> tVsc = mesh.Vsc();
    const scalarField& V = tVsc();

    const surfaceScalarField::Boundary& phiBf =
        phi.boundaryField();

    const scalarField& phiCorrIf = phiCorr;
    const surfaceScalarField::Boundary& phiCorrBf =
        phiCorr.boundaryField();

    slicedSurfaceScalarField lambda
    (
        IOobject
        (
            "lambda",
            mesh.time().timeName(),
            mesh,
            IOobject::NO_READ,
            IOobject::NO_WRITE,
            false
        ),
        mesh,
        dimless,
        allLambda, //1.0
        false   // Use slices for the couples
    );

    scalarField& lambdaIf = lambda;
    surfaceScalarField::Boundary& lambdaBf =
        lambda.boundaryFieldRef();

    scalarField psiMaxn(psiIf.size()); //psiIf： 低阶精度的通量
    scalarField psiMinn(psiIf.size()); //

    psiMaxn = psiMin; //0.0
    psiMinn = psiMax; //1.0

    scalarField sumPhip(psiIf.size(), 0.0);
    scalarField mSumPhim(psiIf.size(), 0.0);

    forAll(phiCorrIf, facei)
    {
        const label own = owner[facei];
        const label nei = neighb[facei];

        psiMaxn[own] = max(psiMaxn[own], psiIf[nei]);
        psiMinn[own] = min(psiMinn[own], psiIf[nei]);

        psiMaxn[nei] = max(psiMaxn[nei], psiIf[own]);
        psiMinn[nei] = min(psiMinn[nei], psiIf[own]);

        const scalar phiCorrf = phiCorrIf[facei];

        if (phiCorrf > 0) //流出为正
        {
            sumPhip[own] += phiCorrf; // sum the outflow A (anti-diffusion flux) to the owner 
            mSumPhim[nei] += phiCorrf; //sum the inflow A to the neiborring cell 
        }
        else
        {
            mSumPhim[own] -= phiCorrf; // change to the positive sign, sum the inflow A (anti-diffusion flux) to the owner 
            sumPhip[nei] -= phiCorrf; // change to the positive sign, sum the outflow A (anti-diffusion flux) to the owner 
        }
    }

    forAll(phiCorrBf, patchi)
    {
        const fvPatchScalarField& psiPf = psiBf[patchi];
        const scalarField& phiCorrPf = phiCorrBf[patchi];

        const labelList& pFaceCells = mesh.boundary()[patchi].faceCells();

        if (psiPf.coupled())
        {
            const scalarField psiPNf(psiPf.patchNeighbourField());

            forAll(phiCorrPf, pFacei)
            {
                label pfCelli = pFaceCells[pFacei];

                psiMaxn[pfCelli] = max(psiMaxn[pfCelli], psiPNf[pFacei]);
                psiMinn[pfCelli] = min(psiMinn[pfCelli], psiPNf[pFacei]);
            }
        }
        else if (psiPf.fixesValue())
        {
            forAll(phiCorrPf, pFacei)
            {
                const label pfCelli = pFaceCells[pFacei];

                psiMaxn[pfCelli] = max(psiMaxn[pfCelli], psiPf[pFacei]);
                psiMinn[pfCelli] = min(psiMinn[pfCelli], psiPf[pFacei]);
            }
        }
        else
        {
            // Add the optional additional allowed boundary extrema
            if (boundaryDeltaExtremaCoeff > 0)
            {
                forAll(phiCorrPf, pFacei)
                {
                    const label pfCelli = pFaceCells[pFacei];

                    const scalar extrema =
                        boundaryDeltaExtremaCoeff
                       *(psiMax[pfCelli] - psiMin[pfCelli]);

                    psiMaxn[pfCelli] += extrema;
                    psiMinn[pfCelli] -= extrema;
                }
            }
        }

        forAll(phiCorrPf, pFacei)
        {
            const label pfCelli = pFaceCells[pFacei];

            const scalar phiCorrf = phiCorrPf[pFacei];

            if (phiCorrf > 0)
            {
                sumPhip[pfCelli] += phiCorrf;
            }
            else
            {
                mSumPhim[pfCelli] -= phiCorrf;
            }
        }
    }

    psiMaxn = min(psiMaxn + extremaCoeff*(psiMax - psiMin), psiMax);
    psiMinn = max(psiMinn - extremaCoeff*(psiMax - psiMin), psiMin);

    if (smoothLimiter > small)
    {
        psiMaxn =
            min(smoothLimiter*psiIf + (1.0 - smoothLimiter)*psiMaxn, psiMax); //相当于做了个松弛
        psiMinn =
            max(smoothLimiter*psiIf + (1.0 - smoothLimiter)*psiMinn, psiMin);
    }

    psiMaxn =
        V
       *(
           (rho.field()*rDeltaT - Sp.field())*psiMaxn
         - Su.field()
         - rho.field()*psi.primitiveField()*rDeltaT
        );

    psiMinn =
        V
       *(
           Su.field()
         - (rho.field()*rDeltaT - Sp.field())*psiMinn
         + rho.field()*psi.primitiveField()*rDeltaT
        );

    scalarField sumlPhip(psiIf.size());
    scalarField mSumlPhim(psiIf.size());

    for (int j=0; j<nLimiterIter; j++)
    {
        sumlPhip = 0;
        mSumlPhim = 0;

        forAll(lambdaIf, facei)
        {
            const label own = owner[facei];
            const label nei = neighb[facei];

            const scalar lambdaPhiCorrf = lambdaIf[facei]*phiCorrIf[facei];

            if (lambdaPhiCorrf > 0)
            {
                sumlPhip[own] += lambdaPhiCorrf; //outflow away from the owner due to the corrected anti-diffusion flow
                mSumlPhim[nei] += lambdaPhiCorrf;//inflow into the neighboring due to the corrected anti-diffusion flow
            }
            else
            {
                mSumlPhim[own] -= lambdaPhiCorrf; //inflow into the owner due to the corrected anti-diffusion flow
                sumlPhip[nei] -= lambdaPhiCorrf;//outflow away from the neighboring due to the corrected anti-diffusion flow
            }
        }

        forAll(lambdaBf, patchi)
        {
            scalarField& lambdaPf = lambdaBf[patchi];
            const scalarField& phiCorrfPf = phiCorrBf[patchi];

            const labelList& pFaceCells = mesh.boundary()[patchi].faceCells();

            forAll(lambdaPf, pFacei)
            {
                label pfCelli = pFaceCells[pFacei];

                scalar lambdaPhiCorrf = lambdaPf[pFacei]*phiCorrfPf[pFacei];

                if (lambdaPhiCorrf > 0)
                {
                    sumlPhip[pfCelli] += lambdaPhiCorrf;
                }
                else
                {
                    mSumlPhim[pfCelli] -= lambdaPhiCorrf;
                }
            }
        }

        forAll(sumlPhip, celli)
        {
            sumlPhip[celli] = 
                max(min
                (
                    (sumlPhip[celli] + psiMaxn[celli]) //sumlPhip: outflow by the corrected anti-diffusion (last step), psiMaxn: maximum allowed inflow
                   /(mSumPhim[celli] + rootVSmall), //mSumPhim: inflow by the anti-diffusion
                    1.0), 0.0
                );

            mSumlPhim[celli] =
                max(min
                (
                    (mSumlPhim[celli] + psiMinn[celli]) //mSumlPhim: inflow by the corrected anti-diffusion(last step)), psiMinn: maximumm allowed outflow
                   /(sumPhip[celli] + rootVSmall), //sumPhip: outflow by the anti-diffusion
                    1.0), 0.0
                );
        }

        const scalarField& lambdam = sumlPhip; //for inflow by the anti-diffusion
        const scalarField& lambdap = mSumlPhim; //for outflow by the anti-diffusion

        forAll(lambdaIf, facei)
        {
            if (phiCorrIf[facei] > 0) //anti-diffusion is outflow
            {
                lambdaIf[facei] = min
                (
                    lambdaIf[facei],
                    min(lambdap[owner[facei]], lambdam[neighb[facei]]) //owner used the lamda for outflow, nei used the lambda for inflow
                );
            }
            else
            {
                lambdaIf[facei] = min
                (
                    lambdaIf[facei],
                    min(lambdam[owner[facei]], lambdap[neighb[facei]])
                );
            }
        }


        forAll(lambdaBf, patchi)
        {
            fvsPatchScalarField& lambdaPf = lambdaBf[patchi];
            const scalarField& phiCorrfPf = phiCorrBf[patchi];
            const fvPatchScalarField& psiPf = psiBf[patchi];

            if (isA<wedgeFvPatch>(mesh.boundary()[patchi]))
            {
                lambdaPf = 0;
            }
            else if (psiPf.coupled())
            {
                const labelList& pFaceCells =
                    mesh.boundary()[patchi].faceCells();

                forAll(lambdaPf, pFacei)
                {
                    const label pfCelli = pFaceCells[pFacei];

                    if (phiCorrfPf[pFacei] > 0)
                    {
                        lambdaPf[pFacei] =
                            min(lambdaPf[pFacei], lambdap[pfCelli]);
                    }
                    else
                    {
                        lambdaPf[pFacei] =
                            min(lambdaPf[pFacei], lambdam[pfCelli]);
                    }
                }
            }
            else
            {
                const labelList& pFaceCells =
                    mesh.boundary()[patchi].faceCells();
                const scalarField& phiPf = phiBf[patchi];

                forAll(lambdaPf, pFacei)
                {
                    // Limit outlet faces only
                    if ((phiPf[pFacei] + phiCorrfPf[pFacei]) > small*small)
                    {
                        const label pfCelli = pFaceCells[pFacei];

                        if (phiCorrfPf[pFacei] > 0)
                        {
                            lambdaPf[pFacei] =
                                min(lambdaPf[pFacei], lambdap[pfCelli]);
                        }
                        else
                        {
                            lambdaPf[pFacei] =
                                min(lambdaPf[pFacei], lambdam[pfCelli]);
                        }
                    }
                }
            }
        }

        syncTools::syncFaceList(mesh, allLambda, minEqOp<scalar>());
    }
}
```
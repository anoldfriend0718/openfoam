{
    const scalar MO2=composition.Wi(O2Index);
    const scalar MCO2=composition.Wi(CO2Index);
    const scalar MN2=composition.Wi(N2Index);
    const volScalarField& T=thermo.T();
    forAll(p,i)
    {
        scalar Ti=T[i];
        scalar pi=p[i]/atmpressure; //atm

        scalar TstartiO2N2=Ti/epsO2N2Kb;
        scalar omegaO2N2=1.6036/std::pow(TstartiO2N2,0.15610) + 0.193/std::exp(0.47635*TstartiO2N2)
                            + 1.03587/std::exp(1.52996*TstartiO2N2)+1.76474/std::exp(3.89411*TstartiO2N2);
        scalar DO2N2i=1.8583e-7*std::sqrt(Ti*Ti*Ti*(1/MO2+1/MN2))/(pi*sigmaO2N2*sigmaO2N2*omegaO2N2);

        scalar TstartiO2CO2=Ti/epsO2CO2Kb;
        scalar omegaO2CO2=1.6036/std::pow(TstartiO2CO2,0.15610) + 0.193/std::exp(0.47635*TstartiO2CO2)
                            + 1.03587/std::exp(1.52996*TstartiO2CO2)+1.76474/std::exp(3.89411*TstartiO2CO2);
        scalar DO2CO2i=1.8583e-7*std::sqrt(Ti*Ti*Ti*(1/MO2+1/MCO2))/(pi*sigmaO2CO2*sigmaO2CO2*omegaO2CO2);

        scalar TstartiCO2N2=Ti/epsCO2N2Kb;
        scalar omegaCO2N2=1.6036/std::pow(TstartiCO2N2,0.15610) + 0.193/std::exp(0.47635*TstartiCO2N2)
                            + 1.03587/std::exp(1.52996*TstartiCO2N2)+1.76474/std::exp(3.89411*TstartiCO2N2);
        scalar DCO2N2i=1.8583e-7*std::sqrt(Ti*Ti*Ti*(1/MCO2+1/MN2))/(pi*sigmaCO2N2*sigmaCO2N2*omegaCO2N2);

        scalar YO2i=Y[O2Index][i];
        scalar YCO2i=Y[CO2Index][i];
        scalar YN2i=Y[N2Index][i];
        scalar totalMolei=YCO2i/MCO2+YO2i/MO2+YN2i/MN2;
        scalar xO2i=YO2i/MO2/totalMolei;
        scalar xCO2i=YCO2i/MCO2/totalMolei;
        scalar xN2i=YN2i/MN2/totalMolei;


        scalar DO2bi=(1-xO2i)/(xCO2i/DO2CO2i+xN2i/DO2N2i);
        scalar DCO2bi=(1-xCO2i)/(xO2i/DO2CO2i+xN2i/DCO2N2i);

        scalar DO2Ki=2*dp/3.0*std::sqrt(2*Rg1000*Ti/PI/MO2);
        scalar DCO2Ki=2*dp/3.0*std::sqrt(2*Rg1000*Ti/PI/MCO2);

        if(cokeRegion[i]>1.0-SMALL) //==1.0
        {
            DeffO2[i]=eps[i]*DO2bi*DO2Ki/(DO2bi+DO2Ki);
            DeffCO2[i]=eps[i]*DCO2bi*DCO2Ki/(DCO2bi+DCO2Ki);
        }
        else
        {
            DeffO2[i]=eps[i]*DO2bi;
            DeffCO2[i]=eps[i]*DCO2bi;
        }                    
    }

    const volScalarField::Boundary& TBf = thermo.T().boundaryFieldRef();
    const volScalarField::Boundary& pBf = p.boundaryFieldRef();
    const volScalarField::Boundary& YO2Bf = Y[O2Index].boundaryFieldRef();
    const volScalarField::Boundary& YCO2Bf = Y[CO2Index].boundaryFieldRef();
    const volScalarField::Boundary& YN2Bf = Y[N2Index].boundaryFieldRef();
    const volScalarField::Boundary& cokeRegionBf = cokeRegion.boundaryFieldRef();
    const volScalarField::Boundary& epsBf = eps.boundaryFieldRef();
    volScalarField::Boundary& DeffO2Bf = DeffO2.boundaryFieldRef();
    volScalarField::Boundary& DeffCO2Bf = DeffCO2.boundaryFieldRef();
    
    forAll(TBf,patchi)
    {
        forAll(TBf[patchi],facei)
        {
            scalar Ti=TBf[patchi][facei];
            scalar pi=pBf[patchi][facei]/101325; //atm

            scalar TstartiO2N2=Ti/epsO2N2Kb;
            scalar omegaO2N2=1.6036/std::pow(TstartiO2N2,0.15610) + 0.193/std::exp(0.47635*TstartiO2N2)
                            + 1.03587/std::exp(1.52996*TstartiO2N2)+1.76474/std::exp(3.89411*TstartiO2N2);
            scalar DO2N2i=1.8583e-7*std::sqrt(Ti*Ti*Ti*(1/MO2+1/MN2))/(pi*sigmaO2N2*sigmaO2N2*omegaO2N2);

            scalar TstartiO2CO2=Ti/epsO2CO2Kb;
            scalar omegaO2CO2=1.6036/std::pow(TstartiO2CO2,0.15610) + 0.193/std::exp(0.47635*TstartiO2CO2)
                            + 1.03587/std::exp(1.52996*TstartiO2CO2)+1.76474/std::exp(3.89411*TstartiO2CO2);
            scalar DO2CO2i=1.8583e-7*std::sqrt(Ti*Ti*Ti*(1/MO2+1/MCO2))/(pi*sigmaO2CO2*sigmaO2CO2*omegaO2CO2);

            scalar TstartiCO2N2=Ti/epsCO2N2Kb;
            scalar omegaCO2N2=1.6036/std::pow(TstartiCO2N2,0.15610) + 0.193/std::exp(0.47635*TstartiCO2N2)
                            + 1.03587/std::exp(1.52996*TstartiCO2N2)+1.76474/std::exp(3.89411*TstartiCO2N2);
            scalar DCO2N2i=1.8583e-7*std::sqrt(Ti*Ti*Ti*(1/MCO2+1/MN2))/(pi*sigmaCO2N2*sigmaCO2N2*omegaCO2N2);

            scalar YO2i=YO2Bf[patchi][facei];
            scalar YCO2i=YCO2Bf[patchi][facei];
            scalar YN2i=YN2Bf[patchi][facei];
            scalar totalMolei=YCO2i/MCO2+YO2i/MO2+YN2i/MN2;
            scalar xO2i=YO2i/MO2/totalMolei;
            scalar xCO2i=YCO2i/MCO2/totalMolei;
            scalar xN2i=YN2i/MN2/totalMolei;

            scalar DO2bi=(1-xO2i)/(xCO2i/DO2CO2i+xN2i/DO2N2i);
            scalar DCO2bi=(1-xCO2i)/(xO2i/DO2CO2i+xN2i/DCO2N2i);

            scalar DO2Ki=2*dp/3.0*std::sqrt(2*Rg1000*Ti/PI/MO2);
            scalar DCO2Ki=2*dp/3.0*std::sqrt(2*Rg1000*Ti/PI/MCO2);

            if(cokeRegionBf[patchi][facei]>1.0-SMALL) //==1.0
            {
                DeffO2Bf[patchi][facei]=epsBf[patchi][facei]*DO2bi*DO2Ki/(DO2bi+DO2Ki);
                DeffCO2Bf[patchi][facei]=epsBf[patchi][facei]*DCO2bi*DCO2Ki/(DCO2bi+DCO2Ki);
            }
            else
            {
                DeffO2Bf[patchi][facei]=epsBf[patchi][facei]*DO2bi;
                DeffCO2Bf[patchi][facei]=epsBf[patchi][facei]*DCO2bi;
            }
        }
    }

}
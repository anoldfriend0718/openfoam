    const scalar PI=3.1415926;
    const scalar PI2=3.1415926*3.1415926;
    const scalar Rg=8.314;
    const scalar Rg1000=Rg*1000;
    const scalar epsO2Kb=113;
    const scalar epsN2Kb=99.8;
    const scalar epsCO2Kb=190;
    const scalar epsO2N2Kb=std::sqrt(epsO2Kb*epsN2Kb);
    const scalar epsCO2N2Kb=std::sqrt(epsCO2Kb*epsN2Kb);
    const scalar epsO2CO2Kb=std::sqrt(epsO2Kb*epsCO2Kb);

    const scalar sigmaO2=3.433;
    const scalar sigmaCO2=3.996;
    const scalar sigmaN2=3.667;

    const scalar sigmaO2CO2=(sigmaO2+sigmaCO2)/2.0;
    const scalar sigmaO2N2=(sigmaO2+sigmaN2)/2.0;
    const scalar sigmaCO2N2=(sigmaCO2+sigmaN2)/2.0;
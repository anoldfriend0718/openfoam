# Workflow for each cell
1. get the deltaT, mass fractions of fluid species,fluid density, coke fraction, temperature 
2. calculate the molar concentrations for fluid species and coke ,solid
3. back up the initial molar concentrations
4. time-split 
4.1 one chemical calculation
- compose the gas_mixture Cpf: gas_mixture Cpf += (this->specieThermo_[i].W()*c[i])*this->specieThermo_[i].Cp(p,T) by loop N2,CO2,CO
- calculate the solid mixture Cps: coke molar concentration* coke molar weight*cp,coke + rock molar concentration* rock molar weight*cp,rock
- calculate the current ha: (Cps + Cpf)*T
- calculate the transient reaction rate, RR, with the current temperature, oxygen molar concentration and fluid density
- calculate the tMin=(current c_O2+small)/RR, and get the chemical time step: deltaT = min(deltaT, tMin);
- solve the c_O2 of next step by implicit way: (av*k + 1/deltaT)* (c_O2 )= c_02_0/deltaT; delta_c_O2=(c_O2-c_02_0) [less than zero]
- update the c_CO2, coke molar concentration : c_CO2 -=delta_c_O2 [increase], c_coke +=delta_c_O2 [decrease]
- update the coke fraction: coke_fraction -=(delta_c_coke*W_coke)/coke mass density 
- update the gas mixture
- get the heat capacity of gas mixture, Cpf, and coke 
- calculate the adabatic temperature: ((1-eps)*coke molar concentration* molar weight*cp +eps* gas_mixture Cpf )* T= ha
- return updated delataT, molar concentration, adabatic temperature
4.2 increae the time by the deltaT, get the timeleft, and go to 4.1
4.3 till chemistry time is equal to the transport time

5. calculate the effective chemical reaction rate 


// CEqn
    // HCl reaction rate
    mdotci = -mag(fvc::grad(eps))*(4.0*eps*(1.0-eps))*alpha1*r_si/C; //Foam::max(C,Cmin);  //fvm::Sp(mdotci,C); -fvm::Su(mdotci*Ceq,C)
    mdotc = mdotci*(C-Ceq)*C;  
// SEqn
    // solid dissolution rate
    mdots = mdotc*gama_si*Mc/rho_si;  //fvm::Su(mdots,eps_si)
// alphaEqn
    mdotle = mdotc*gama_l*Mc/rho1; //explicit mdotl/rhol
    mdotli =  -mdotc*(gama_si/rho_si+gama_l/rho1+gama_g/rho2)*Mc;
    divU = fvc::div(fvc::absolute(phi, U));  //fvc::div(phi,U);
// pEqn
   mdotp = -mdotli;
   mdotp1=mdotp;   //debug



//Permeability Models
k0 = K*(eps*eps*eps)/(max((1-eps)*(1-eps),SMALL)); //Uncomment for Kozeny-Carman
k0f = fvc::interpolate(k0,"k0");

//Labeling the places that have porous media
forAll(eps,celli) 
{
	if(eps[celli]>0.99) 
	{	
	Solid[celli]=0; 
	}
	else
	{
	Solid[celli]=1; 
	}
	if(eps[celli]<=0.01) 
	{	
	Solidw[celli]=1; 
	}
	else
	{
	Solidw[celli]=0; 
	}
}

forAll(mesh.boundary(), patchi) 
{
	forAll(Solid.boundaryField()[patchi],facei)
	{
		if(eps.boundaryField()[patchi][facei]>0.99) 
		{Solid.boundaryFieldRef()[patchi][facei]=0;}
		else
		{Solid.boundaryFieldRef()[patchi][facei]=1;}
	}
}

Solidf = fvc::interpolate(Solid);
Solidwf = fvc::interpolate(Solidw);
Solidwf0 = Solidwf;
// Solidf0 = Solidf;
// Solid1 = fvc::average(fvc::interpolate(Solid));
//Eliminating intermediate values of Solid Indicator Function
forAll(Solidf,facei)
{
	if(Solidf[facei]<1) // if <1 ==0.  
	{	
	Solidf[facei]=0; 
	}
}

forAll(Solidwf,facei)
{
	if(Solidwf[facei]>0) // if >0 ==1.  
	{	
	Solidwf[facei]=1; 
	}
}

forAll(Solidwf0,facei)
{
	if(Solidwf0[facei]<1) // if <1 ==0.  
	{	
	Solidwf0[facei]=0; 
	}
}


// forAll(Solid1,celli)
// {
// 	if(Solid1[celli]>0) // if <1 ==0.  
// 	{	
// 	Solid1[celli]=1; 
// 	}
// }

Solidv = max((1.0-Solid),SMALL);
Solidf0 = fvc::interpolate(Solidv);
forAll(Solidf0,facei)
{
	if(Solidf0[facei]<1) // if <1 ==0.  
	{	
	Solidf0[facei]=0; 
	}
}
Fluidf = max(1.0-Solidf,SMALL);
Fluidwf = max(1.0-Solidwf0,SMALL);

Fluidf2 = fvc::interpolate(fvc::average(Fluidf));
forAll(Fluidf2,facei)
{
	if(Fluidf2[facei]<1) // if <1 ==0.  
	{	
	Fluidf2[facei]=1; 
	}
}





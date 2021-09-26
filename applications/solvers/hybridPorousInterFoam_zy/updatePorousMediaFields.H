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
Solidf0 = Solidf;

//Eliminating intermediate values of Solid Indicator Function
forAll(Solidf,facei)
{
	if(Solidf[facei]<1) // if <1 ==0.  
	{	
	Solidf[facei]=0; 
	}
}

Solidv = max((1.0-Solidw),SMALL);



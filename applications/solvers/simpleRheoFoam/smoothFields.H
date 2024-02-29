scalar wDiv = 0.1;
for (label i = 0; i < nGamma; i++)
{
    edgeScalarField divUsf(fac::interpolate(divUs));
    divUs = wDiv*fac::average(divUsf) + (1. - wDiv)*divUs;
    edgeScalarField gammae(fac::interpolate(gamma));
    gamma = wGamma*fac::average(gammae) + (1. - wGamma)*gamma;
}

Info<< "t_res: " << runTime.timeName() 
    << " gamma: " << residuals.gamma() 
    << " Us: " << residuals.Us() 
    << " iters: " << looping.index() << endl;

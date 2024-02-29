for (label i = 0; i < nGamma; i++)
{
    edgeScalarField gammae(fac::interpolate(gamma));
    gamma = wGamma*fac::average(gammae) + (1. - wGamma)*gamma;
}

//Info<< "t_res: " << runTime.timeName()
//    << " " << residuals.gamma() 
//    << " " << residuals.Us() << endl;

dictionary solDict = aMesh.solutionDict().subDict("solution");
scalar tol = solDict.lookupOrDefault<scalar>("tolerance", 1e-06);
bool coupledToBulk = solDict.lookupOrDefault<bool>("coupledToBulk", true);
label nGamma = solDict.lookupOrDefault<label>("nGamma", 0);
scalar wGamma = solDict.lookupOrDefault<scalar>("wGamma", 0.0);

interfaceModel->coupledToBulk(coupledToBulk);

rheoResiduals residuals(tol);

// Check under-relaxation factor
dictionary relaxDict = aMesh.solutionDict().subDict("relaxationFactors");
scalar alphaUs = 
    relaxDict.subDict("equations").lookupOrDefault<scalar>("Us", 1);
scalar alphaUsRheo = 1;
scalar theta = interfaceModel->Theta().value();
if (theta > 0)
{
    alphaUsRheo = 1./theta;
}
else
{
    alphaUsRheo = interfaceModel->Bq().value()/interfaceModel->Bqk().value();
}

if (alphaUs > alphaUsRheo)
{
    Info<<"Warning: alphaUs > " << alphaUsRheo 
        << " may lead to convergence issues." << endl;
}


dictionary solDict = aMesh.solutionDict().subDict("solution");
scalar tol = solDict.lookupOrDefault<scalar>("tolerance", 1e-06);
bool coupledToBulk = solDict.lookupOrDefault<bool>("coupledToBulk", true);
label nGamma = solDict.lookupOrDefault<label>("nGamma", 0);
scalar wGamma = solDict.lookupOrDefault<scalar>("wGamma", 0.0);

interfaceModel->coupledToBulk(coupledToBulk);

rheoResiduals residuals(tol);

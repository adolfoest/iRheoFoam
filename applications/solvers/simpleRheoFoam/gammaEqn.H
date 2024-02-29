faScalarMatrix gammaEqn
(
    fam::ddt(gamma)
  + fam::div(phis, gamma)
  - fam::laplacian(1./interfaceModel->Pes(), gamma)
);

gammaEqn.solve();
gamma.relax();
gamma.correctBoundaryConditions();

Info<< "max(gamma): " << gMax(gamma.primitiveField()) 
    << " min(gamma): " << gMin(gamma.primitiveField()) 
    << endl;

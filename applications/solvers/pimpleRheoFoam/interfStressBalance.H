tmp<faVectorMatrix> tUsEqn
(
    interfaceModel->divTaus(Us)
  + interfaceModel->marangoni(gamma)
 ==
    interfaceModel->nTaut(U, Us)
);
faVectorMatrix& UsEqn = tUsEqn.ref();

UsEqn.relax();
UsEqn.solve();

// Calculate phis
Us.correctBoundaryConditions();
if (isDynMesh)
{
    mfam.correctBC(mesh, Us);
}
Usf = fac::interpolate(Us);
phis = Usf & aMesh.Le();

//correct Usf
if (isDynMesh)
{
    mfam.correct(Usf, phis);
}

// Correct pressure gradient at interface
if (coupledToBulk)
{
    if
    (
        p.boundaryField()[interfaceModel->patchID()].type()
    == fixedGradientFvPatchField<scalar>::typeName
    )
    {
        fixedGradientFvPatchField<scalar>& pb =
            refCast<fixedGradientFvPatchField<scalar>>
            (
                p.boundaryFieldRef()[interfaceModel->patchID()]
            ); 
        pb.gradient() = interfaceModel->snGradTau(U);
    }
}

// Calculate div(Us)
divUs = fac::div(phis);
scalar maxDivUs = max(mag(min(divUs.primitiveField())), max(divUs.primitiveField()));
scalar avDivUs = gAverage(divUs.primitiveField());
reduce(maxDivUs, maxOp<scalar>());
Info<< "max(div(Us)): " << maxDivUs << " average(div(Us)): " << avDivUs << endl;

if (maxDivUs > 200)
{
    FatalErrorInFunction << "Divergence error: max(div(Us)) =  "
        << maxDivUs
        << exit(FatalError);
}
// Make the fluxes relative to the mesh motion
if (isDynMesh)
{
    mfam.makeRelative(mesh, phis);
}

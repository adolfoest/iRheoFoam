if (runTime.writeTime())
{
    // vsm.mapToVolume(divUs, divU.boundaryFieldRef());
    const label& patchID = interfaceModel->patchID();
    forAll(mesh.boundary()[patchID], i)
    {
        const label& bCell = mesh.boundary()[patchID].faceCells()[i];
        gammaVol[bCell] = gamma[i];
        gammaVol.boundaryFieldRef()[patchID][i] = gamma[i];  
        UsVol[bCell] = Us[i];
        UsVol.boundaryFieldRef()[patchID][i] = Us[i];
        
        // Rest of face boundaries next to the interface
        forAll(mesh.cells()[bCell], j)
        {
            const label& faceJ = mesh.cells()[bCell][j];
            if (!mesh.isInternalFace(faceJ))
            {
                const label& pf = mesh.boundaryMesh().whichPatch(faceJ);
                if (pf != patchID)
                {
                    const label& k = mesh.boundaryMesh()[pf].whichFace(faceJ);
                    gammaVol.boundaryFieldRef()[pf][k] = gamma[i];
                    UsVol.boundaryFieldRef()[pf][k] = Us[i];
                }
            }
        }
    }
}

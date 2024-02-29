//#include "createControl.H"
pimpleControlRheo pimple(mesh);
#include "createTimeControls.H"
 
bool correctPhi
(
    pimple.dict().getOrDefault("correctPhi", mesh.dynamic())
);
 
bool checkMeshCourantNo
(
    pimple.dict().getOrDefault("checkMeshCourantNo", false)
);
 
bool moveMeshOuterCorrectors
(
    pimple.dict().getOrDefault("moveMeshOuterCorrectors", false)
);

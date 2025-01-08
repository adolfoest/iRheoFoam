{
    const scalar& deltaT = runTime.deltaTValue();
    const scalar& Ma = interfaceModel->Ma().value();
    const scalar& Pes = interfaceModel->Pes().value();

    const scalar&  maxDeltaTRheo(min(min(1/Ma,Pes),1));
    if (deltaT > maxDeltaTRheo)
    {
        Info<< "Warning: deltaT > min(min(1/Ma,Pes),1) "
            << "may affect stability of calculations." << endl;
    }

}

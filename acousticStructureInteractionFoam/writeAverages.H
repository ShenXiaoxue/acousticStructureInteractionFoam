    fieldAcousticEnergyFile << runTime.timeName() << " "
    << fvc::domainIntegrate(p * p / (rhoAcoustic * c02)).value() << " "
    << fvc::domainIntegrate(mag(p * (-fvc::grad(phi)))/c0).value()
        << " " << endl;

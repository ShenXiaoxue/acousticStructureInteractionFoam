    if (runTime.outputTime())
    {
        volSymmTensorField sigma
        (
            IOobject
            (
                "sigma",
                runTime.timeName(),
                stressMesh,
                IOobject::NO_READ,
                IOobject::AUTO_WRITE
            ),
            rho*sigmaU
        );

        if (thermalStress)
        {
            const volScalarField& T = Tptr();
            sigma = sigma - I*(rho*threeKalpha*T);
        }

        volScalarField sigmaEq
        (
            IOobject
            (
                "sigmaEq",
                runTime.timeName(),
                stressMesh,
                IOobject::NO_READ,
                IOobject::AUTO_WRITE
            ),
            sqrt((3.0/2.0)*magSqr(dev(sigma)))
        );

        Info<< "Max sigmaEq = " << max(sigmaEq).value()
            << endl;

        runTime.write();
    }

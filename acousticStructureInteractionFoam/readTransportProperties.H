    Info<< "Reading transportProperties\n" << endl;
    IOdictionary transportProperties
    (
        IOobject
        (
            "transportProperties",
            runTime.constant(),
            mesh,
            IOobject::MUST_READ,
            IOobject::NO_WRITE
        )
    );

    Info<< "Reading speed of sound c0\n" << endl;
    dimensionedScalar c0
    (
        transportProperties.lookup("c0")
    );

    Info<< "Reading density of medium rho\n" << endl;
    dimensionedScalar rhoAcoustic
    (
        transportProperties.lookup("rhoAcoustic")
    );

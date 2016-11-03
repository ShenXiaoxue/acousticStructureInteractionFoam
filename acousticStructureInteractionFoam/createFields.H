    Info<< "Reading field phi\n" << endl;
    volScalarField phi
    (
        IOobject
        (
            "phi",
            runTime.timeName(),
            mesh,
            IOobject::MUST_READ,
            IOobject::AUTO_WRITE
        ),
        mesh
    );

#if WRITE_PROBES
    Info<< "Creating field p\n" << endl;
//    volScalarField p("p", fvc::ddt(phi) * rhoAcoustic);

    volScalarField p
    (
        IOobject
        (
            "p",
            runTime.timeName(),
            mesh,
            IOobject::NO_READ,
            IOobject::AUTO_WRITE
        ),
	fvc::ddt(phi) * rhoAcoustic     
    );

#endif



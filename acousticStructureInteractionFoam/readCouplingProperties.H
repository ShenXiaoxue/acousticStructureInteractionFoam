    Info << "Reading coupling properties" << endl;
    IOdictionary couplingProperties
    (
        IOobject
        (
            "couplingProperties",
            runTime.constant(),
            mesh,
            IOobject::MUST_READ,
            IOobject::NO_WRITE
        )
    );

    // Read solid patch data
    word solidPatchName(couplingProperties.lookup("solidPatch"));

    label solidPatchID =
        stressMesh.boundaryMesh().findPatchID(solidPatchName);


    // Read fluid patch data
    word acousticPatchName(couplingProperties.lookup("acousticPatch"));

    label acousticPatchID =
        mesh.boundaryMesh().findPatchID(acousticPatchName);


    if (solidPatchID < 0 || acousticPatchID < 0)
    {
        FatalErrorIn(args.executable())
            << "Problem with patch interpolation definition"
            << abort(FatalError);
    }


    // Create interpolators
    patchToPatchInterpolation interpolatorAcousticSolid
    (
        mesh.boundaryMesh()[acousticPatchID],
        stressMesh.boundaryMesh()[solidPatchID]
    );

    patchToPatchInterpolation interpolatorSolidAcoustic
    (
        stressMesh.boundaryMesh()[solidPatchID],
        mesh.boundaryMesh()[acousticPatchID]
    );


    // Grab solid patch field
    tractionDisplacementFvPatchVectorField& tForce =
        refCast<tractionDisplacementFvPatchVectorField>
        (
            Usolid.boundaryField()[solidPatchID]
        );


    // Grab motion field

    // Read fluid patch data
    word movingRegionName(couplingProperties.lookup("movingRegion"));

    const fvMesh& motionMesh =
        runTime.objectRegistry::lookupObject<fvMesh>(movingRegionName);

    tetPointVectorField& motionU =
        const_cast<tetPointVectorField&>
        (
            motionMesh.objectRegistry::lookupObject<tetPointVectorField>
            (
                "motionU"
            )
        );


    fixedValueTetPolyPatchVectorField& motionUAcousticPatch =
        refCast<fixedValueTetPolyPatchVectorField>
        (
            motionU.boundaryField()[acousticPatchID]
        );

    tetPolyPatchInterpolation tppi
    (
        refCast<const faceTetPolyPatch>(motionUAcousticPatch.patch())
    );

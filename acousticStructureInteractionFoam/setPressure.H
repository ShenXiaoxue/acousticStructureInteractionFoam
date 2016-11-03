{
    // Setting pressure on solid patch
    Info << "Setting pressure" << endl;

    scalarField solidPatchPressure =
        interpolatorAcousticSolid.faceInterpolate
        (
            p.boundaryField()[acousticPatchID]
        );

    solidPatchPressure *= rhoAcoustic.value();

    tForce.pressure() = solidPatchPressure;


    vector totalPressureForce =
        sum
        (
            p.boundaryField()[acousticPatchID]*
            mesh.Sf().boundaryField()[acousticPatchID]
        );


    Info << "Total pressure force = " << totalPressureForce << endl;
}

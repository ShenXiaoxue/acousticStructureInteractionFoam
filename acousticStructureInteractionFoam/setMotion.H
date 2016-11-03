{
    // Setting mesh motion

    pointVectorField solidPointsDispl =
        cpi.interpolate(Usolid - Usolid.oldTime());

    vectorField newPoints =
        stressMesh.points()
      + solidPointsDispl.internalField();

    stressMesh.movePoints(newPoints);


    vectorField acousticPatchPointsDispl =
        interpolatorSolidAcoustic.pointInterpolate
        (
            solidPointsDispl.boundaryField()[solidPatchID].
            patchInternalField()
        );

    motionUAcousticPatch ==
        tppi.pointToPointInterpolate
        (
            acousticPatchPointsDispl/runTime.deltaT().value()
        );

    vibrationInducedPressureFvPatchScalarField& tAccelera =
        refCast<vibrationInducedPressureFvPatchScalarField>
        (
            phi.boundaryField()[acousticPatchID]
        );

/*
    for (int i=0;i<128;i++)
    {
            Info << "acousticPatchPointsDispl  " << i <<"  "<< acousticPatchPointsDispl[i].z() << endl;
    }
*/

    forAll(phi.boundaryField()[acousticPatchID], faceI) 
    { 
            tAccelera.displacement0()[faceI].x() = acousticPatchPointsDispl[faceI].x(); 
            tAccelera.displacement0()[faceI].y() = acousticPatchPointsDispl[faceI].y();
            tAccelera.displacement0()[faceI].z() = acousticPatchPointsDispl[faceI].z();
    } 


//    tAccelera.displacement0() = acousticPatchPointsDispl.z();


    mesh.update();
    /* test code    
    if
    (
        phi.boundaryField()[0].type()
     != virbrationInducedPressureFvPatchField<Type>::typeName
	//!= tractionDisplacementIncrementFvPatchVectorField::typeName
    )
    {
        FatalErrorIn(args.executable())
            << "Bounary condition on " << phi.name()
                <<  " is "
                << phi.boundaryField()[0].type()
                << "for fluid -solid interface patch, instead "
	  //	        << virbrationInducedPressureFvPatchField<Type>::typeName
	        //<< tractionDisplacementIncrementFvPatchVectorField::typeName
                << abort(FatalError);
    }    
 */
 //   virbrationInducedPressureFvPatchField<>& patTest;

#   include "volContinuity.H"

    Info << "Motion magnitude: mean = "
        << average(mag(Usolid.boundaryField()[solidPatchID]))
        << " max = "
        << max(mag(Usolid.boundaryField()[solidPatchID])) << endl;
}

/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | foam-extend: Open Source CFD
   \\    /   O peration     |
    \\  /    A nd           | For copyright notice see file Copyright
     \\/     M anipulation  |
-------------------------------------------------------------------------------
License
    This file is part of foam-extend.

    foam-extend is free software: you can redistribute it and/or modify it
    under the terms of the GNU General Public License as published by the
    Free Software Foundation, either version 3 of the License, or (at your
    option) any later version.

    foam-extend is distributed in the hope that it will be useful, but
    WITHOUT ANY WARRANTY; without even the implied warranty of
    MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the GNU
    General Public License for more details.

    You should have received a copy of the GNU General Public License
    along with foam-extend.  If not, see <http://www.gnu.org/licenses/>.

Application
    acousticStructureInteractionFoam

Description
    Transient solver for interaction between acoustic and structure with
    mesh motion.  Set up as a fake acoustic structure interaction solver

\*---------------------------------------------------------------------------*/
#define WRITE_PROBES 1

#include "fvCFD.H"
#include "OFstream.H"
#include "dynamicFvMesh.H"
#include "patchToPatchInterpolation.H"
#include "tractionDisplacement/tractionDisplacementFvPatchVectorField.H"
#include "vibrationInducedPressure/vibrationInducedPressureFvPatchScalarField.H"
#include "tetFemMatrices.H"
#include "tetPointFields.H"
#include "faceTetPolyPatch.H"
#include "tetPolyPatchInterpolation.H"
#include "fixedValueTetPolyPatchFields.H"

#include "pointMesh.H"
#include "pointFields.H"
#include "volPointInterpolation.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

int main(int argc, char *argv[])
{
#   include "setRootCase.H"
#   include "createTime.H"
#   include "createDynamicFvMesh.H"
#   include "createStressMesh.H"
#   include "readTransportProperties.H"
#   include "createFields.H"
#   include "readStressedFoamControls.H"
#   include "readMechanicalProperties.H"
#   include "readThermalProperties.H"
#   include "createStressFields.H"


#   include "readCouplingProperties.H"
#   include "readTimeControls.H"

/*#   include "createAverages.H" */
//#   include "initContinuityErrs.H"
// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //
    Info<< "\nStarting time loop\n" << endl;

#   include "acousticCourantNo.H"

    dimensionedScalar c02 = c0 * c0;

    while (runTime.run())
    {
//#       include "readPISOControls.H"
//#       include "readTimeControls.H"
//#       include "CourantNo.H"
#       include "setDeltaT.H"

        runTime++;

        Info<< "Time = " << runTime.timeName() << nl << endl;

#       include "setPressure.H"
#       include "solveSolid.H"

#       include "setMotion.H"
#       include "solveAcoustic.H"

        runTime.write();

        Info<< "ExecutionTime = "
            << runTime.elapsedCpuTime()
            << " s\n\n" << endl;
    }

    Info<< "End\n" << endl;

    return(0);
}


// ************************************************************************* //

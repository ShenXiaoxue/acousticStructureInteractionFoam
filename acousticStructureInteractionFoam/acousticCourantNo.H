/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     |
    \\  /    A nd           | Copyright held by original authors.
     \\/     M anipulation  |
-------------------------------------------------------------------------------
License
    OpenFOAM is free software; you can redistribute it and/or modify it
    under the terms of the GNU General Public License as published by the
    Free Software Foundation; either version 2 of the License, or (at your
    option) any later version.

    OpenFOAM is distributed in the hope that it will be useful, but WITHOUT
    ANY WARRANTY; without even the implied warranty of MERCHANTABILITY or
    FITNESS FOR A PARTICULAR PURPOSE.  See the GNU General Public License
    for more details.

    You should have received a copy of the GNU General Public License
    along with OpenFOAM; if not, write to the Free Software Foundation,
    Inc., 51 Franklin St, Fifth Floor, Boston, MA 02110-1301 USA

Global
    acousticCourantNo

Description
    Calculates and outputs the mean and maximum acoustic Courant Numbers.

Copyrights
    (C) 1991-2007 OpenCFD Ltd. (courantNo.H original version)
    (C) 2007 Takuya OSHIMA

\*---------------------------------------------------------------------------*/

scalar CoNum = 0.0;
scalar meanCoNum = 0.0;

if (mesh.nInternalFaces())
{
    CoNum = max(mesh.surfaceInterpolation::deltaCoeffs() * c0)
      .value() * runTime.deltaT().value();

    meanCoNum
      = (sum(mesh.surfaceInterpolation::deltaCoeffs() * c0 * mesh.magSf())
	 /sum(mesh.magSf())).value() * runTime.deltaT().value();
}

Info<< "Mean and max acoustic Courant Numbers = "
    << meanCoNum << " " << CoNum << endl;

// ************************************************************************* //

/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     |
    \\  /    A nd           | Copyright (C) 2007 Takuya OSHIMA
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

Description
    Compatibility symbols across different versions of OpenFOAM.

Author
    Takuya OSHIMA <oshima@eng.niigata-u.ac.jp>

\*---------------------------------------------------------------------------*/

#ifndef convectiveOutletConfig_H
#define convectiveOutletConfig_H

// OF >= 1.4 requires a copy ctor
#if WM_PROJECT_VERSION_MAJOR > 1 || (WM_PROJECT_VERSION_MAJOR == 1 \
    && WM_PROJECT_VERSION_MINOR >= 4)
#define NEEDS_COPY_CTOR 1
#else
#define NEEDS_COPY_CTOR 0
#endif

// OF >= 1.4.1 requires ctor with DimensionedField arguments
#if WM_PROJECT_VERSION_MAJOR == 1 && (WM_PROJECT_VERSION_MINOR < 4 \
    || WM_PROJECT_VERSION_MINOR == 4 && WM_PROJECT_VERSION_PATCH == 0)
#define FIELD_IS_DIMENSIONED 0
#else
#define FIELD_IS_DIMENSIONED 1
#endif

#endif

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

Class
    convectiveOutletFvPatchField<Type>

Description
    Convective outlet boundary condition.

Author
    Takuya OSHIMA <oshima@eng.niigata-u.ac.jp>

SourceFiles
    convectiveOutletFvPatchField.C

\*---------------------------------------------------------------------------*/

#ifndef convectiveOutletFvPatchField_H
#define convectiveOutletFvPatchField_H

#include "convectiveOutletConfig.H"

#include "fvPatchFields.H"
#include "fixedGradientFvPatchFields.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{

/*---------------------------------------------------------------------------*\
                 Class convectiveOutletFvPatch Declaration
\*---------------------------------------------------------------------------*/

template <class Type>
class convectiveOutletFvPatchField
:
    public fixedGradientFvPatchField<Type>
{
    // Private data

        //- Convective velocity
        Field<scalar> convectiveVelocity_;

        //- Time index
        label curTimeIndex_;

	//- Surface-normal directional discretization scheme
	word snGradScheme_;

	//- Temporal discretization scheme
	word ddtScheme_;

	//- Update boundary values in updateCoeffs() or not
	bool updateValue_;

	//- Write boundary values or not
	bool writeValue_;

        //- "field::patch" name of myself
        word fieldPatchName_;

        //- Old time patch-internal gradient
	Field<Type> gi0_;

        //- Old time boundary gradient
	Field<Type> gb0_;

        //- Old time patch-internal field
        Field<Type> pi0_;

        //- Old-old time patch-internal field
        Field<Type> pi00_;

        //- Old time boundary patch field
        Field<Type> pb0_;

        //- Old-old time boundary patch field
        Field<Type> pb00_;

public:

    //- Runtime type information
    TypeName("convectiveOutlet");


    // Constructors

        //- Construct from patch and internal field
        convectiveOutletFvPatchField
        (
            const fvPatch&,
#if FIELD_IS_DIMENSIONED
            const DimensionedField<Type, volMesh>&
#else
            const Field<Type>&
#endif
        );

        //- Construct from patch, internal field and dictionary
        convectiveOutletFvPatchField
        (
            const fvPatch&,
#if FIELD_IS_DIMENSIONED
            const DimensionedField<Type, volMesh>&,
#else
            const Field<Type>&,
#endif
            const dictionary&
        );

        //- Construct by mapping given convectiveOutletFvPatchField
        //  onto a new patch
        convectiveOutletFvPatchField
        (
            const convectiveOutletFvPatchField<Type>&,
            const fvPatch&,
#if FIELD_IS_DIMENSIONED
            const DimensionedField<Type, volMesh>&,
#else
            const Field<Type>&,
#endif
            const fvPatchFieldMapper&
        );

#if NEEDS_COPY_CTOR
        //- Construct as copy
        convectiveOutletFvPatchField
        (
            const convectiveOutletFvPatchField<Type>&
        );
#endif

        //- Construct and return a clone
        virtual tmp<fvPatchField<Type> > clone() const
        {
            return tmp<fvPatchField<Type> >
            (
                new convectiveOutletFvPatchField<Type>(*this)
            );
        }

        //- Construct as copy setting internal field reference
        convectiveOutletFvPatchField
        (
            const convectiveOutletFvPatchField<Type>&,
#if FIELD_IS_DIMENSIONED
            const DimensionedField<Type, volMesh>&
#else
            const Field<Type>&
#endif
        );

        //- Construct and return a clone setting internal field reference
        virtual tmp<fvPatchField<Type> > clone
        (
#if FIELD_IS_DIMENSIONED
            const DimensionedField<Type, volMesh>& iF
#else
            const Field<Type>& iF
#endif
        ) const
        {
            return tmp<fvPatchField<Type> >
            (
                new convectiveOutletFvPatchField<Type>(*this, iF)
            );
        }


    // Member functions

        //- Update the coefficients associated with the patch field
        virtual void updateCoeffs();

        //- Write
        virtual void write(Ostream&) const;


    // Member operators

};

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace Foam

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#ifdef NoRepository
#   include "convectiveOutletFvPatchField.C"
#endif

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#endif

// ************************************************************************* //

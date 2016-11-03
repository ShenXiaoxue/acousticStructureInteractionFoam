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
    vibrationInducedPressureFvPatchScalarField

Description
    vibration boundary condition for acoustic caculation.

Author
    Xiaoxue Shen <shenxiaoxue@csu.edu.cn>

SourceFiles
    vibrationInducedPressureFvPatchScalarField.C

\*---------------------------------------------------------------------------*/

#ifndef vibrationInducedPressureFvPatchScalarField_H
#define vibrationInducedPressureFvPatchScalarField_H

#include "fvPatchFields.H"
#include "fixedGradientFvPatchFields.H"

namespace Foam
{

/*---------------------------------------------------------------------------*\
                 Class vibrationInducedPressureFvPatch Declaration
\*---------------------------------------------------------------------------*/

class vibrationInducedPressureFvPatchScalarField
:
    public fixedGradientFvPatchScalarField
{
  // Private data

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

  //- Old time boundary gradient
  scalarField gb0_;

  //- Old time patch-internal field
  scalarField pi0_;

  //- Old-old time patch-internal field
  scalarField pi00_;

  //- Old time boundary patch field
  scalarField pb0_;

  //- Old-old time boundary patch field
  scalarField pb00_;

  //- Old time displacement of boundary
  vectorField displacement0_;

  //- Old-old time displacement of boundary
  vectorField displacement00_;

public:

  //- Runtime type information
  TypeName("vibrationInducedPressure");



  // Constructors

  //- Construct from patch and internal field
  vibrationInducedPressureFvPatchScalarField
  (
   const fvPatch&,
//#if FIELD_IS_DIMENSIONED
   const DimensionedField<scalar,volMesh>&
//#else
//   const Field&
//#endif
   );

  //- Construct from patch, internal field and dictionary
  vibrationInducedPressureFvPatchScalarField
  (
   const fvPatch&,
//#if FIELD_IS_DIMENSIONED
   const DimensionedField<scalar,volMesh>&,
//#else
//   const Field&,
//#endif
   const dictionary&
   );

  //- Construct by mapping given vibrationInducedPressureFvPatchScalarField
  //  onto a new patch
  vibrationInducedPressureFvPatchScalarField
  (
   const vibrationInducedPressureFvPatchScalarField&,
   const fvPatch&,
//#if FIELD_IS_DIMENSIONED
   const DimensionedField<scalar, volMesh>&,
//#else 
//   const Field&,
//#endif
   const fvPatchFieldMapper&
   );

#if NEEDS_COPY_CTOR
  //- Construct as copy
  vibrationInducedPressureFvPatchScalarField
  (
   const vibrationInducedPressureFvPatchScalarField&
   );
#endif

  //- Construct and return a clone
  virtual tmp<fvPatchField > clone() const
  {
    return tmp<fvPatchField >
      (
       new vibrationInducedPressureFvPatchScalarField(*this)
       );
  }

  //- Construct as copy setting internal field reference
  vibrationInducedPressureFvPatchScalarField
  (
   const vibrationInducedPressureFvPatchScalarField&,
//#if FIELD_IS_DIMENSIONED
   const DimensionedField<scalar, volMesh>&
//#else
//   const Field&
//#endif
   );

  //- Construct and return a clone setting internal field reference
  virtual tmp<fvPatchField > clone
  (
//#if FIELD_IS_DIMENSIONED
   const DimensionedField<scalar, volMesh>& iF
//#else
//   const Field& iF
//#endif
   ) const
  {
    return tmp<fvPatchField >
      (
       	  new vibrationInducedPressureFvPatchScalarField(*this, iF)
       );
  }

// Member functions

//- Update the coefficients associated with the patch field
virtual void updateCoeffs();

virtual const vectorField& displacement0() const
{
    return displacement0_;
}

virtual  vectorField& displacement0()
{
    return displacement0_;
}

//- Write 
virtual void write(Ostream&) const;

//  Member operators

};

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} //End namespace Foam

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //


#endif

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

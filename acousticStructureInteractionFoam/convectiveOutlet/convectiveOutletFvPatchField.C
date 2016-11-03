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
    Convective outlet boundary condition.

Author
    Takuya OSHIMA <oshima@eng.niigata-u.ac.jp>

\*---------------------------------------------------------------------------*/

#include "convectiveOutletFvPatchField.H"

#include "addToRunTimeSelectionTable.H"
#include "fvPatchFieldMapper.H"
#include "volFields.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{

// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

template<class Type>
convectiveOutletFvPatchField<Type>::convectiveOutletFvPatchField
(
    const fvPatch& p,
#if FIELD_IS_DIMENSIONED
    const DimensionedField<Type, volMesh>& iF
#else
    const Field<Type>& iF
#endif
)
:
    fixedGradientFvPatchField<Type>(p, iF),
    convectiveVelocity_(p.size(), 0.0),
    curTimeIndex_(0),
    snGradScheme_("normal"),
    ddtScheme_("CrankNicholson"),
    updateValue_(false),
    writeValue_(false),
    fieldPatchName_("::" + p.name()),
    gi0_(p.size(), pTraits<Type>::zero),
    gb0_(p.size(), pTraits<Type>::zero),
    pi0_(p.size(), pTraits<Type>::zero),
    pi00_(p.size(), pTraits<Type>::zero),
    pb0_(p.size(), pTraits<Type>::zero),
    pb00_(p.size(), pTraits<Type>::zero)
{
  Info << "1wahaha" << nl << endl;
}


template<class Type>
convectiveOutletFvPatchField<Type>::convectiveOutletFvPatchField
(
    const convectiveOutletFvPatchField<Type>& ptf,
    const fvPatch& p,
#if FIELD_IS_DIMENSIONED
    const DimensionedField<Type, volMesh>& iF,
#else
    const Field<Type>& iF,
#endif
    const fvPatchFieldMapper& mapper
)
:
    fixedGradientFvPatchField<Type>(ptf, p, iF, mapper),
    convectiveVelocity_(ptf.convectiveVelocity_, mapper),
    curTimeIndex_(ptf.curTimeIndex_),
    snGradScheme_(ptf.snGradScheme_),
    ddtScheme_(ptf.ddtScheme_),
    updateValue_(ptf.updateValue_),
    writeValue_(ptf.writeValue_),
    fieldPatchName_(ptf.fieldPatchName_),
    gi0_(ptf.gi0_, mapper),
    gb0_(ptf.gb0_, mapper),
    pi0_(ptf.pi0_, mapper),
    pi00_(ptf.pi00_, mapper),
    pb0_(ptf.pb0_, mapper),
    pb00_(ptf.pb00_, mapper)
{
  Info << "2wahaha" << nl << endl;
}


template<class Type>
convectiveOutletFvPatchField<Type>::convectiveOutletFvPatchField
(
    const fvPatch& p,
#if FIELD_IS_DIMENSIONED
    const DimensionedField<Type, volMesh>& iF,
#else
    const Field<Type>& iF,
#endif
    const dictionary& dict
)
:
    fixedGradientFvPatchField<Type>(p, iF),
    convectiveVelocity_("convectiveVelocity", dict, p.size()),
    curTimeIndex_(0),
    snGradScheme_("normal"),
    ddtScheme_("CrankNicholson"),
    updateValue_(false),
    writeValue_(false),
    fieldPatchName_(dict.name().name()),
    gi0_(p.size(), pTraits<Type>::zero),
    gb0_(p.size(), pTraits<Type>::zero),
    pi0_(p.size(), pTraits<Type>::zero),
    pi00_(p.size(), pTraits<Type>::zero),
    pb0_(p.size(), pTraits<Type>::zero),
    pb00_(p.size(), pTraits<Type>::zero)
{
    Info << "3wahaha" << nl << endl;
    if (dict.found("snGradScheme"))
    {
	dict.lookup("snGradScheme") >> snGradScheme_;

	if(snGradScheme_ != "upwind" && snGradScheme_ != "predictorCorrector"
	   && snGradScheme_ != "normal")
	{
	    FatalErrorIn
	      (
	       "convectiveOutletFvPatchField::convectiveOutletFvPatchField()"
	       )   << "    Unsupported surface-normal differencing scheme : "
		   << snGradScheme_
		   << exit(FatalError);
        }
    }

    if (dict.found("ddtScheme"))
    {
	dict.lookup("ddtScheme") >> ddtScheme_;

	if(ddtScheme_ != "Euler" && ddtScheme_ != "backward"
	   && ddtScheme_ != "CrankNicholson")
	  {
	    FatalErrorIn
	      (
	       "convectiveOutletFvPatchField::convectiveOutletFvPatchField()"
	       )   << "    Unsupported temporal differencing scheme : "
		   << ddtScheme_
		   << exit(FatalError);
	  }
    }

    if (dict.found("updateValue"))
    {
	dict.lookup("updateValue") >> updateValue_;
    }

    if (dict.found("writeValue"))
    {
	dict.lookup("writeValue") >> writeValue_;
    }

    pi00_ = this->patchInternalField();
    pi0_ = this->patchInternalField();
    if (dict.found("gradient"))
    {
  	Info << "gradient right!" << nl <<endl;
        this->gradient() = Field<Type>("gradient", dict, p.size());
	pb00_ = pi00_ + this->gradient()/this->patch().deltaCoeffs();
//	Info << this->patch().deltaCoeffs() << nl << endl;
    }
    else
    {
        Info << "gradient r!" << nl <<endl;
	this->gradient() = Field<Type>(this->size(), pTraits<Type>::zero);
	pb00_ = pi00_;
    }
    pb0_ = pb00_;

    this->operator==(this->patchInternalField()
	  + this->gradient()/this->patch().deltaCoeffs());
}


template<class Type>
convectiveOutletFvPatchField<Type>::convectiveOutletFvPatchField
(
    const convectiveOutletFvPatchField<Type>& coptf,
#if FIELD_IS_DIMENSIONED
    const DimensionedField<Type, volMesh>& iF
#else
    const Field<Type>& iF
#endif
)
:
    fixedGradientFvPatchField<Type>(coptf, iF),
    convectiveVelocity_(coptf.convectiveVelocity_),
    curTimeIndex_(coptf.curTimeIndex_),
    snGradScheme_(coptf.snGradScheme_),
    ddtScheme_(coptf.ddtScheme_),
    updateValue_(coptf.updateValue_),
    writeValue_(coptf.writeValue_),
    fieldPatchName_(coptf.fieldPatchName_),
    gi0_(coptf.gi0_),
    gb0_(coptf.gb0_),
    pi0_(coptf.pi0_),
    pi00_(coptf.pi00_),
    pb0_(coptf.pb0_),
    pb00_(coptf.pb00_)
{
    Info << "4wahaha" << nl << endl;
}

#if NEEDS_COPY_CTOR
template<class Type>
convectiveOutletFvPatchField<Type>::convectiveOutletFvPatchField
(
    const convectiveOutletFvPatchField<Type>& coptf
)
    :
    fixedGradientFvPatchField<Type>(coptf),
    convectiveVelocity_(coptf.convectiveVelocity_),
    curTimeIndex_(coptf.curTimeIndex_),
    snGradScheme_(coptf.snGradScheme_),
    ddtScheme_(coptf.ddtScheme_),
    updateValue_(coptf.updateValue_),
    writeValue_(coptf.writeValue_),
    fieldPatchName_(coptf.fieldPatchName_),
    gi0_(coptf.gi0_),
    gb0_(coptf.gb0_),
    pi0_(coptf.pi0_),
    pi00_(coptf.pi00_),
    pb0_(coptf.pb0_),
    pb00_(coptf.pb00_)
{
    Info << "5wahaha" << nl << endl;
}
#endif

// * * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * //


// Update the coefficients associated with the patch field
template<class Type>
void convectiveOutletFvPatchField<Type>::updateCoeffs()
{
    Info << "6wahaha" << nl << endl;
    if (this->updated())
    {
        Info << "Boundary gradient of " << fieldPatchName_
	     << " has already been updated" << nl;
        return;
    }

    const Time& ldb = this->db().time();

    Info << "Updating boundary gradient of " << fieldPatchName_ << " by ";

    if (ddtScheme_ == "Euler")
    {
        Info << "Euler and ";

	if(snGradScheme_ == "normal")
        {
            Info << "normal schemes" << nl;
	    this->gradient() = (pb0_ - this->patchInternalField())
	      /(convectiveVelocity_*ldb.deltaT().value()
		+ 1.0/this->patch().deltaCoeffs());
        }
	else if (snGradScheme_ == "upwind")
        {
	  Info << "upwind schemes" << nl;
	  this->gradient() = -(this->patchInternalField() - pi0_)
	    /(convectiveVelocity_*ldb.deltaT().value());
	}
	else if(snGradScheme_ == "predictorCorrector")
	{
	  Info << "predictor-corrector schemes" << nl;
	  // predictor
	  this->gradient() = -(this->patchInternalField() - pi0_)
	    /(convectiveVelocity_*ldb.deltaT().value());
	  // corrector
	  this->gradient() = 0.5 * (this->gradient()
	  - ((this->patchInternalField()
	      + this->gradient()/this->patch().deltaCoeffs()) - pb0_)
			 /(convectiveVelocity_*ldb.deltaT().value()));
	}
    }
    else if (ddtScheme_ == "backward")
    {
        Info << "backward differencing and ";

	if(snGradScheme_ == "normal")
        {
            Info << "normal schemes" << nl;
	    this->gradient() = (4.0*pb0_ - pb00_
				- 3.0*this->patchInternalField())
	      /(2.0*convectiveVelocity_*ldb.deltaT().value()
		+ 3.0/this->patch().deltaCoeffs());
        }
	else if (snGradScheme_ == "upwind")
        {
            Info << "upwind schemes" << nl;
	    this->gradient() = -(3.0*this->patchInternalField() - 4.0*pi0_
		       + pi00_)/(2.0*convectiveVelocity_*ldb.deltaT().value());
	}
	else if(snGradScheme_ == "predictorCorrector")
	{
            Info << "predictor-corrector schemes" << nl;
	    // predictor
	    this->gradient() = -(3.0*this->patchInternalField() - 4.0*pi0_
		       + pi00_)/(2.0*convectiveVelocity_*ldb.deltaT().value());
	    // corrector
	    this->gradient() = 0.5*(this->gradient()
	    - (3.0*(this->patchInternalField()
		  + this->gradient()/this->patch().deltaCoeffs())
	   - 4.0*pb0_ + pb00_)/(2.0*convectiveVelocity_*ldb.deltaT().value()));
	}
    }
    else if (ddtScheme_ == "CrankNicholson")
    {
        Info << "CrankNicholson and ";

	if(snGradScheme_ == "normal")
        {
            Info << "normal schemes" << nl;
	    this->gradient() = (2.0*(pb0_ - this->patchInternalField())
	      - convectiveVelocity_*ldb.deltaT().value()*gb0_)
	      /(convectiveVelocity_*ldb.deltaT().value()
		+ 2.0/this->patch().deltaCoeffs());
        }
	else if (snGradScheme_ == "upwind")
        {
	  Info << "upwind schemes" << nl;
	  this->gradient() = -2.0*(this->patchInternalField() - pi0_)
	    /(convectiveVelocity_*ldb.deltaT().value()) - gb0_;
	}
	else if(snGradScheme_ == "predictorCorrector")
	{
	  Info << "predictor-corrector schemes" << nl;
	  // predictor
	  this->gradient() = -2.0*(this->patchInternalField() - pi0_)
	    /(convectiveVelocity_*ldb.deltaT().value()) - gi0_;
          // Use time index to save oldTime gradient values
          if (curTimeIndex_ != ldb.timeIndex())
          {
              gi0_ = this->gradient();
          }
	  // corrector
	  this->gradient() = 0.5*(this->gradient()
	  - 2.0*((this->patchInternalField()
	      + this->gradient()/this->patch().deltaCoeffs()) - pb0_)
			 /(convectiveVelocity_*ldb.deltaT().value()) - gb0_);
	}
    }

    if(updateValue_)
    {
        this->operator==(this->patchInternalField()
			 + this->gradient()/this->patch().deltaCoeffs());
    }

    // Use time index to save oldTime patchField values
    if (curTimeIndex_ != ldb.timeIndex())
    {
        Info << "Saving old-time boundary values of " << fieldPatchName_ << nl;
	gb0_ = this->gradient();

	if (ddtScheme_ == "backward")
	{
	    pi00_ = pi0_;
	    pb00_ = pb0_;
	}
        pi0_ = this->patchInternalField();
        pb0_ = this->patchInternalField()
	  + this->gradient()/this->patch().deltaCoeffs();
        curTimeIndex_ = ldb.timeIndex();
    }

    // Sets Updated to true
    fixedGradientFvPatchField<Type>::updateCoeffs();
}


// Write
template<class Type>
void convectiveOutletFvPatchField<Type>::write(Ostream& os) const
{
    Info << "7wahaha" << nl << endl;
    fixedGradientFvPatchField<Type>::write(os);
    os.writeKeyword("snGradScheme") << snGradScheme_ << token::END_STATEMENT
				    << endl;
    os.writeKeyword("ddtScheme") << ddtScheme_ << token::END_STATEMENT << endl;
    os.writeKeyword("updateValue") << updateValue_ << token::END_STATEMENT
				    << endl;
    os.writeKeyword("writeValue") << writeValue_ << token::END_STATEMENT
				    << endl;
    convectiveVelocity_.writeEntry("convectiveVelocity", os);

    if (writeValue_)
    {
        this->writeEntry("value", os);
    }
}

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace Foam

// ************************************************************************* //

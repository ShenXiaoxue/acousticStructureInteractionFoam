#include "vibrationInducedPressureFvPatchScalarField.H"

#include "addToRunTimeSelectionTable.H"
#include "fvPatchFieldMapper.H"
#include "volFields.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{

// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //


  vibrationInducedPressureFvPatchScalarField::vibrationInducedPressureFvPatchScalarField
  (
   const fvPatch& p,
//#if FIELD_IS_DIMENSIONED
   const DimensionedField<scalar, volMesh>& iF
//#else
//   const Field<scalar>& iF
//#endif
   )
    :
    fixedGradientFvPatchScalarField(p,iF),
    curTimeIndex_(0),
    snGradScheme_("normal"),
    ddtScheme_("CrankNicholson"),
    updateValue_(false),
    writeValue_(false),
    fieldPatchName_("::" + p.name()),
    gb0_(p.size(), pTraits<scalar>::zero),
    pi0_(p.size(), pTraits<scalar>::zero),
    pi00_(p.size(), pTraits<scalar>::zero),
    pb0_(p.size(), pTraits<scalar>::zero),
    pb00_(p.size(), pTraits<scalar>::zero),
    displacement0_(p.size(), vector::zero),
    displacement00_(p.size(), vector::zero)
  {
    Info << "1chushihua？"<< nl << endl;
    this->gradient() = 0.0;
  }


  vibrationInducedPressureFvPatchScalarField::vibrationInducedPressureFvPatchScalarField
  (
   const vibrationInducedPressureFvPatchScalarField& ptf,
   const fvPatch& p,
//#if FIELD_IS_DIMENSIONED
   const DimensionedField<scalar, volMesh>& iF,
//#else
//   const Field<scalar>& iF
//#endif
   const fvPatchFieldMapper& mapper
   )
    :
   fixedGradientFvPatchScalarField(ptf, p, iF, mapper),
    curTimeIndex_(ptf.curTimeIndex_),
    snGradScheme_(ptf.snGradScheme_),
    ddtScheme_(ptf.ddtScheme_),
    updateValue_(ptf.updateValue_),
    writeValue_(ptf.writeValue_),
    fieldPatchName_(ptf.fieldPatchName_),
   gb0_(ptf.gb0_, mapper),
   pi0_(ptf.pi0_, mapper),
   pi00_(ptf.pi00_, mapper),
   pb0_(ptf.pb0_, mapper),
   pb00_(ptf.pb00_, mapper),
   displacement0_(ptf.displacement0_, mapper),
   displacement00_(ptf.displacement00_, mapper)
  {
    Info << "2chushihua？"<< nl << endl;
  }

vibrationInducedPressureFvPatchScalarField::vibrationInducedPressureFvPatchScalarField
(
    const fvPatch& p,
//#if FIELD_IS_DIMENSIONED
    const DimensionedField<scalar, volMesh>& iF,
//#else
//    const Field<scalar>& iF,
//#endif
    const dictionary& dict
)
:
    fixedGradientFvPatchScalarField(p, iF),
    curTimeIndex_(0),
    snGradScheme_("normal"),
    ddtScheme_("CrankNicholson"),
    updateValue_(false),
    writeValue_(false),
    fieldPatchName_(dict.name().name()),
    gb0_(p.size(), pTraits<scalar>::zero),
    pi0_(p.size(), pTraits<scalar>::zero),
    pi00_(p.size(), pTraits<scalar>::zero),
    pb0_(p.size(), pTraits<scalar>::zero),
    pb00_(p.size(), pTraits<scalar>::zero),
    displacement0_("displacement", dict, p.size()),
    displacement00_(p.size(), vector::zero)
{
   Info << "3chushihua？"<< nl << endl;

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
        this->gradient() = Field<scalar>("gradient", dict, p.size());
   	pb00_ = pi00_ + this->gradient()/this->patch().deltaCoeffs();
//	Info << this->patch().deltaCoeffs() << nl << endl;
    }
    else
    {
        Info << "gradient r!" << nl <<endl;
	this->gradient() = Field<scalar>(this->size(), pTraits<scalar>::zero);
    	pb00_ = pi00_;
    }

    pb0_ = pb00_;

    this->operator==(this->patchInternalField()
	  + this->gradient()/this->patch().deltaCoeffs());
}   
   

vibrationInducedPressureFvPatchScalarField::vibrationInducedPressureFvPatchScalarField
(
    const vibrationInducedPressureFvPatchScalarField& coptf,
//#if FIELD_IS_DIMENSIONED
    const DimensionedField<scalar, volMesh>& iF
//#else
//    const Field<scalar>& iF
//#endif
)
:
    fixedGradientFvPatchScalarField(coptf, iF),
    curTimeIndex_(coptf.curTimeIndex_),
    snGradScheme_(coptf.snGradScheme_),
    ddtScheme_(coptf.ddtScheme_),
    updateValue_(coptf.updateValue_),
    writeValue_(coptf.writeValue_),
    fieldPatchName_(coptf.fieldPatchName_),
    gb0_(coptf.gb0_),
    pi0_(coptf.pi0_),
    pi00_(coptf.pi00_),
    pb0_(coptf.pb0_),
    pb00_(coptf.pb00_),
   displacement0_(coptf.displacement0_),
   displacement00_(coptf.displacement00_)
  {
    Info << "4chushihua？"<< nl << endl;
  }

#if NEEDS_COPY_CTOR
vibrationInducedPressureFvPatchScalarField::vibrationInducedPressureFvPatchScalarField
(
    const vibrationInducedPressureFvPatchScalarField& coptf
)
    :
    fixedGradientFvPatchScalarField(coptf),
    curTimeIndex_(coptf.curTimeIndex_),
    snGradScheme_(coptf.snGradScheme_),
    ddtScheme_(coptf.ddtScheme_),
    updateValue_(coptf.updateValue_),
    writeValue_(coptf.writeValue_),
    fieldPatchName_(coptf.fieldPatchName_),
    gb0_(coptf.gb0_),
    pi0_(coptf.pi0_),
    pi00_(coptf.pi00_),
    pb0_(coptf.pb0_),
    pb00_(coptf.pb00_),
    displacement0_(coptf.displacement0_),
    displacement00_(coptf.displacement00_)
{
    Info << "5chushihua？" << nl << endl;
}
#endif

// * * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * //

// Update the coefficients associated with the patch field
   void vibrationInducedPressureFvPatchScalarField::updateCoeffs()
   {
     Info << "6chushihua？" << nl <<endl;
	
     if (this->updated())
       {
	 Info << "Boundary gradient of " << fieldPatchName_
	      << " has already been updated" << nl;
	 return;
       }
     const Time& ldb = this->db().time();

     Info << "Updating boundary gradient of " << fieldPatchName_ << endl;
     
//     const int rho(7854);
//     const float rhoAcoustic(1.205);
     vectorField m = patch().nf();

     //     this->gradient() = - rhoAcoustic/rho *((displacement0_-displacement00_)&m)/ldb.deltaT().value() + gb0_;
     //this->gradient() = -((displacement0_-displacement00_)&m)/ldb.deltaT().value() + gb0_;
     this->gradient() = -displacement0_&m/ldb.deltaT().value();
   if(updateValue_)
     {
   this->operator==(this->patchInternalField()+ this->gradient()/this->patch().deltaCoeffs());
   }
   // Use time index to save oldTime patchField values
   if(curTimeIndex_ !=ldb.timeIndex())
     {
   Info << "Saving old-time boundary values of " << fieldPatchName_ << nl;
   gb0_ = this->gradient();
   pi00_ = pi0_;
   pb00_ = pb0_;
   pi0_ = this->patchInternalField();
   pb0_ = this->patchInternalField() + this->gradient()/this->patch().deltaCoeffs();

   curTimeIndex_ = ldb.timeIndex();

   displacement00_ = this->displacement0_;
   }

   // Sets Updated to true
   fixedGradientFvPatchScalarField::updateCoeffs();
   }
     
// Write
void vibrationInducedPressureFvPatchScalarField::write(Ostream& os) const
{
    Info << "7chushihua？" << nl << endl;
    fixedGradientFvPatchScalarField::write(os);
    os.writeKeyword("snGradScheme") << snGradScheme_ << token::END_STATEMENT
				    << endl;
    os.writeKeyword("ddtScheme") << ddtScheme_ << token::END_STATEMENT << endl;
    os.writeKeyword("updateValue") << updateValue_ << token::END_STATEMENT
				    << endl;
    os.writeKeyword("writeValue") << writeValue_ << token::END_STATEMENT
				    << endl;

    if (writeValue_)
    {
        this->writeEntry("value", os);
    }
}

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

makePatchTypeField(fvPatchScalarField, vibrationInducedPressureFvPatchScalarField);

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //
// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace Foam

// ************************************************************************* //     

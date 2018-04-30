/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     |
    \\  /    A nd           | Copyright (C) 2013-2016 OpenFOAM Foundation
     \\/     M anipulation  |
-------------------------------------------------------------------------------
License
    This file is part of OpenFOAM.

    OpenFOAM is free software: you can redistribute it and/or modify it
    under the terms of the GNU General Public License as published by
    the Free Software Foundation, either version 3 of the License, or
    (at your option) any later version.

    OpenFOAM is distributed in the hope that it will be useful, but WITHOUT
    ANY WARRANTY; without even the implied warranty of MERCHANTABILITY or
    FITNESS FOR A PARTICULAR PURPOSE.  See the GNU General Public License
    for more details.

    You should have received a copy of the GNU General Public License
    along with OpenFOAM.  If not, see <http://www.gnu.org/licenses/>.

\*---------------------------------------------------------------------------*/

#include "fieldBasedGradientFvPatchField.H"

// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

template<class Type>
Foam::fieldBasedGradientFvPatchField<Type>::fieldBasedGradientFvPatchField
(
    const fvPatch& p,
    const DimensionedField<Type, volMesh>& iF
)
:
    fixedGradientFvPatchField<Type>(p, iF),
    fieldName_("field"),
    gradCoeff_(0.0)
{}


template<class Type>
Foam::fieldBasedGradientFvPatchField<Type>::fieldBasedGradientFvPatchField
(
    const fvPatch& p,
    const DimensionedField<Type, volMesh>& iF,
    const Field<Type>& fld
)
:
    fixedGradientFvPatchField<Type>(p, iF, fld),
    fieldName_("field"),
    gradCoeff_(0.0)
{}


template<class Type>
Foam::fieldBasedGradientFvPatchField<Type>::fieldBasedGradientFvPatchField
(
    const fvPatch& p,
    const DimensionedField<Type, volMesh>& iF,
    const dictionary& dict
)
:
    fixedGradientFvPatchField<Type>(p, iF, dict),
    fieldName_(dict.lookupOrDefault<word>("field", "field")),
    gradCoeff_(readScalar(dict.lookup("gradCoeff")))
{}


template<class Type>
Foam::fieldBasedGradientFvPatchField<Type>::fieldBasedGradientFvPatchField
(
    const fieldBasedGradientFvPatchField<Type>& ptf,
    const fvPatch& p,
    const DimensionedField<Type, volMesh>& iF,
    const fvPatchFieldMapper& mapper
)
:
    fixedGradientFvPatchField<Type>(ptf, p, iF, mapper),
    fieldName_(ptf.fieldName_),
    gradCoeff_(ptf.gradCoeff_)
{}


template<class Type>
Foam::fieldBasedGradientFvPatchField<Type>::fieldBasedGradientFvPatchField
(
    const fieldBasedGradientFvPatchField<Type>& ptf
)
:
    fixedGradientFvPatchField<Type>(ptf),
    fieldName_(ptf.fieldName_),
    gradCoeff_(ptf.gradCoeff_)
{}


template<class Type>
Foam::fieldBasedGradientFvPatchField<Type>::fieldBasedGradientFvPatchField
(
    const fieldBasedGradientFvPatchField<Type>& ptf,
    const DimensionedField<Type, volMesh>& iF
)
:
    fixedGradientFvPatchField<Type>(ptf, iF),
    fieldName_(ptf.fieldName_),
    gradCoeff_(ptf.gradCoeff_)
{}


// * * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * //

template<class Type>
void Foam::fieldBasedGradientFvPatchField<Type>::updateCoeffs()
{
    if (this->updated())
    {
        return;
    }

    const Field<Type>& fieldp = this->patch().template lookupPatchField
        <
            GeometricField<Type, fvPatchField, volMesh>,
            Type
        >
    (
        fieldName_
    );
    
    this->gradient() = gradCoeff_ * fieldp;
    
    fixedGradientFvPatchField<Type>::updateCoeffs();
}


template<class Type>
void Foam::fieldBasedGradientFvPatchField<Type>::write(Ostream& os) const
{
    fixedGradientFvPatchField<Type>::write(os);
    os.writeKeyword("field") << fieldName_ << token::END_STATEMENT << nl;
    os.writeKeyword("gradCoeff") << gradCoeff_ << token::END_STATEMENT << nl;
    this->writeEntry("value", os);
}



// ************************************************************************* //

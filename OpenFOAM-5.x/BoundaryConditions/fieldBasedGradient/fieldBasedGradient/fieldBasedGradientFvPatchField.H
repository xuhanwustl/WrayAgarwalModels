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

Class
    Foam::fieldBasedGradientFvPatchField

Group
    grpGenericBoundaryConditions

Description
    This boundary condition provides a gradient condition based on a field.
    It is originally developed for applying rough wall boundary condition 
    required by roughness corrected Spalart-Allmaras and Wray-Agarwal model.

Usage
    Example of the boundary condition specification:
    \verbatim
    <patchName>
    {
        type            fieldBasedGradient;
        gradient        uniform 0;
        field           U;
        gradCoeff       2.0;
    }
    \endverbatim

Note
    The example above gives the usage for supplying:
        initial gradient value
        field name
        gradient coefficient
        
    The gradient is calculated as gradient = gradCoeff * field

See also
    Foam::fixedGradientFvPatchField

SourceFiles
    fieldBasedGradientFvPatchField.C

\*---------------------------------------------------------------------------*/

#ifndef fieldBasedGradientFvPatchField_H
#define fieldBasedGradientFvPatchField_H

#include "fixedGradientFvPatchFields.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{

/*---------------------------------------------------------------------------*\
              Class fieldBasedGradientFvPatchField Declaration
\*---------------------------------------------------------------------------*/

template<class Type>
class fieldBasedGradientFvPatchField
:
    public fixedGradientFvPatchField<Type>
{
    // Private data
    
        //- Name of the field field
        word fieldName_;
        
        //- Gradient coefficient
        scalar gradCoeff_;


public:

    //- Runtime type information
    TypeName("fieldBasedGradient");


    // Constructors

        //- Construct from patch and internal field
        fieldBasedGradientFvPatchField
        (
            const fvPatch&,
            const DimensionedField<Type, volMesh>&
        );

        //- Construct from patch and internal field and patch field
        fieldBasedGradientFvPatchField
        (
            const fvPatch&,
            const DimensionedField<Type, volMesh>&,
            const Field<Type>& fld
        );

        //- Construct from patch, internal field and dictionary
        fieldBasedGradientFvPatchField
        (
            const fvPatch&,
            const DimensionedField<Type, volMesh>&,
            const dictionary&
        );

        //- Construct by mapping given fieldBasedGradientFvPatchField
        //  onto a new patch
        fieldBasedGradientFvPatchField
        (
            const fieldBasedGradientFvPatchField<Type>&,
            const fvPatch&,
            const DimensionedField<Type, volMesh>&,
            const fvPatchFieldMapper&
        );

        //- Construct as copy
        fieldBasedGradientFvPatchField
        (
            const fieldBasedGradientFvPatchField<Type>&
        );

        //- Construct and return a clone
        virtual tmp<fvPatchField<Type>> clone() const
        {
            return tmp<fvPatchField<Type>>
            (
                new fieldBasedGradientFvPatchField<Type>(*this)
            );
        }

        //- Construct as copy setting internal field reference
        fieldBasedGradientFvPatchField
        (
            const fieldBasedGradientFvPatchField<Type>&,
            const DimensionedField<Type, volMesh>&
        );

        //- Construct and return a clone setting internal field reference
        virtual tmp<fvPatchField<Type>> clone
        (
            const DimensionedField<Type, volMesh>& iF
        ) const
        {
            return tmp<fvPatchField<Type>>
            (
                new fieldBasedGradientFvPatchField<Type>(*this, iF)
            );
        }


    // Member functions

        //- Update the coefficients associated with the patch field
        virtual void updateCoeffs();

        //- Write
        virtual void write(Ostream&) const;
};


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace Foam

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#ifdef NoRepository
    #include "fieldBasedGradientFvPatchField.C"
#endif

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#endif

// ************************************************************************* //

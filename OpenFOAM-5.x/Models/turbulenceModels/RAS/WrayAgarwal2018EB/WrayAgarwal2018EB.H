/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     |
    \\  /    A nd           | Copyright (C) 2011-2015 OpenFOAM Foundation
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
    Foam::RASModels::WrayAgarwal2018EB

Group
    grpRASTurbulence

Description
    Wray-Agarwal one equation Turbulence Model
    An elliptic blended wall distance free version of the WA-2017 model on 
    NASA Turbulence Modeling Resource (TMR) website 
    
    Coupled the Rnu equation with the elliptic blending equation (PR equation) 
    to improved the near wall velocity profile
    
    Require the coefficient Lref_ input as reference length scale in the model 
    coefficient dictionary
    
    Re-defined the f1 function in the WA-2017 model without using the wall 
    distance
    
    Added a bound on the second destruction term in the Rnu transport 
    equation to fix the significant eddy viscosity drop in the zero strain 
    rate region (e.g. channel center)
    
    Reference:
    \verbatim
        https://turbmodels.larc.nasa.gov/wray_agarwal.html
        Han, X., Rahman, M. M., and Agarwal, R. K., “Wall-Distance Free 
        Wray-Agarwal Turbulence Model with Elliptic Blending,” AIAA Aviation Meeting, 
        Atlanta, GA, June 2018.
    \endverbatim

    The default model coefficients are
        WrayAgarwal2018EBCoeffs
        {
            kappa       0.41
            Cw          5.97
            C1ke        0.094
            C1kw        0.2
            sigmake     1.0
            sigmakw     0.72
            C2ke        1.24
            C2kw        2.63
            C3kw        0.171
            Cmu         0.09
            Cm          8.0
            Lref        <reference length scale>
        }
    
SourceFiles
    WrayAgarwal2018EB.C

\*---------------------------------------------------------------------------*/

#ifndef WrayAgarwal2018EB_H
#define WrayAgarwal2018EB_H

#include "WrayAgarwal2018.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{
namespace RASModels
{

/*---------------------------------------------------------------------------*\
                      Class WrayAgarwal2018EB Declaration
\*---------------------------------------------------------------------------*/

template<class BasicTurbulenceModel>
class WrayAgarwal2018EB
:
    public WrayAgarwal2018<BasicTurbulenceModel>
{
    // Private Member Functions


protected:

    // Protected data

        // Model coefficients
        
            dimensionedScalar C3kw_;
            dimensionedScalar Lref_;

        // Fields

            volScalarField PR_;


    // Protected Member Functions
    
        tmp<volScalarField> LR2() const;

public:

    using typename BasicTurbulenceModel::transportModel;
    using typename BasicTurbulenceModel::rhoField;
    using typename BasicTurbulenceModel::alphaField;

    using WrayAgarwalRASModel<BasicTurbulenceModel>::C2kw_;
    using WrayAgarwalRASModel<BasicTurbulenceModel>::C2ke_;
    using WrayAgarwalRASModel<BasicTurbulenceModel>::Rnu_;
    using WrayAgarwalRASModel<BasicTurbulenceModel>::f1_;
    using WrayAgarwalRASModel<BasicTurbulenceModel>::S_;
    
    using WrayAgarwal2018<BasicTurbulenceModel>::Cm_;
    using WrayAgarwal2018<BasicTurbulenceModel>::W_;
    
    
    //- Runtime type information
    TypeName("WrayAgarwal2018EB");


    // Constructors

        //- Construct from components
        WrayAgarwal2018EB
        (
            const alphaField& alpha,
            const rhoField& rho,
            const volVectorField& U,
            const surfaceScalarField& alphaRhoPhi,
            const surfaceScalarField& phi,
            const transportModel& transport,
            const word& propertiesName = turbulenceModel::propertiesName,
            const word& type = typeName
        );

        // Disallow default bitwise copy construct and assignment
        WrayAgarwal2018EB(const WrayAgarwal2018EB&) = delete;
        WrayAgarwal2018EB& operator=(const WrayAgarwal2018EB&) = delete;


    //- Destructor
    virtual ~WrayAgarwal2018EB()
    {}


    // Member Functions

        //- Re-read model coefficients if they have changed
        virtual bool read();

        //- Solve the turbulence equations and correct the turbulence viscosity
        virtual void correct();
};


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace RASModels
} // End namespace Foam

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#ifdef NoRepository
#   include "WrayAgarwal2018EB.C"
#endif

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#endif

// ************************************************************************* //

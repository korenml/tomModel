/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     |
    \\  /    A nd           | Copyright (C) 2015-2016 OpenFOAM Foundation
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
    Foam::RASModels::zetaf0Base

Group
    grpRASTurbulence

Description
    Abstract base-class for zeta-f models to provide BCs access to the zeta and f
    fields.

See also
    Foam::RASModels::zetaf0

SourceFiles
    zetaf0Base.C

\*---------------------------------------------------------------------------*/

#ifndef zetaf0Base_H
#define zetaf0Base_H

#include "RASModel.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{
namespace RASModels
{

/*---------------------------------------------------------------------------*\
                            Class zetaf0Base Declaration
\*---------------------------------------------------------------------------*/

class zetaf0Base
{
public:

    //- Runtime type information
    TypeName("zetaf0Base");


    // Constructors

        zetaf0Base()
        {}


    //- Destructor
    virtual ~zetaf0Base()
    {}


    // Member Functions

        //- Return turbulence stress normal to streamlines
        virtual tmp<volScalarField> zeta() const = 0;

        //- Return the damping function
        virtual tmp<volScalarField> f0() const = 0;
};


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace RASModels
} // End namespace Foam

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //


#endif

// ************************************************************************* //

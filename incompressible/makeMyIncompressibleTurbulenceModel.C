/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     |
    \\  /    A nd           | Copyright (C) 2013-2015 OpenFOAM Foundation
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

#include "IncompressibleTurbulenceModel.H"
#include "incompressible/transportModel/transportModel.H"
#include "addToRunTimeSelectionTable.H"
#include "makeTurbulenceModel.H"

#include "laminarModel.H"
#include "RASModel.H"
#include "LESModel.H"


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

/*
Header based on: $FOAM_SRC/TurbulenceModels/incompressible/
				turbulentTransportModels/turbulentTransportModel.H
*/

namespace Foam
{
    typedef IncompressibleTurbulenceModel<transportModel> 
        transportModelIncompressibleTurbulenceModel;
    typedef RASModel<transportModelIncompressibleTurbulenceModel>
        RAStransportModelIncompressibleTurbulenceModel;
}

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#define makeRASModel(Type)                                                     \
    makeTemplatedTurbulenceModel                                               \
    (transportModelIncompressibleTurbulenceModel, RAS, Type)

#define makeLESModel(Type)                                                     \
    makeTemplatedTurbulenceModel                                               \
    (transportModelIncompressibleTurbulenceModel, LES, Type)

// -------------------------------------------------------------------------- //
// RAS models
// -------------------------------------------------------------------------- //

#include "kEpsilonPANSHR.H"
makeRASModel(kEpsilonPANSHR);

#include "kEpsilonPANSLR.H"
makeRASModel(kEpsilonPANSLR);

#include "kEpsilonPANSLRN.H"
makeRASModel(kEpsilonPANSLRN);

#include "kEpsilonLRN.H"
makeRASModel(kEpsilonLRN);

#include "kEpsilonPANSdyn.H"
makeRASModel(kEpsilonPANSdyn);

#include "kEpsilonPANSvar1.H"
makeRASModel(kEpsilonPANSvar1);

#include "kOmegaPANS.H"
makeRASModel(kOmegaPANS);

#include "kOmegaPANSLR.H"
makeRASModel(kOmegaPANSLR);

#include "kOmegaPANSdyn.H"
makeRASModel(kOmegaPANSdyn);

#include "kOmegaPANSdynFY.H"
makeRASModel(kOmegaPANSdynFY);

#include "kOmegaPANSvar1.H"
makeRASModel(kOmegaPANSvar1);


#include "kOmegaTNTPANS.H"
makeRASModel(kOmegaTNTPANS);

#include "kOmegaTNTPANSdyn.H"
makeRASModel(kOmegaTNTPANSdyn);

#include "kOmegaTNTPANSSSV.H"
makeRASModel(kOmegaTNTPANSSSV);

#include "kEpsilonPANSSSV.H"
makeRASModel(kEpsilonPANSSSV);

#include "zetaf0VarTL.H"
makeRASModel(zetaf0VarTL);

// -------------------------------------------------------------------------- //
// LES models
// -------------------------------------------------------------------------- //


// ----------------------------------------------------------------- end-of-file

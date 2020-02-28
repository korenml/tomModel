/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     |
    \\  /    A nd           | Copyright (C) 2012-2016 OpenFOAM Foundation
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
    Foam::RASModels::zetaf0

Group
    grpRASTurbulence

Description
    


    Wall boundary conditions are:

        k       = kLowReWallFunction
        epsilon = epsilonLowReWallFunction
        zeta      = zetaWallFunction
        f0       = fWallFunction

    These are applicable to both low- and high-Reynolds number flows.

    Inlet values can be approximated by:

        zeta      = 2/3
        f0       = zero-gradient

    References:
    \verbatim

		Hanjalić, K. & Popovac, M. & Hadžiabdić (2004).
		A robust near-wall elliptic-relaxation eddy-viscosity turbulence model for CFD.
		International Journal of Heat and Fluid Flow, 25, 1047-1051.

		Popovac, M (2010).
		Robust eddy viscosity turbulence modeling with elliptic relaxation for external building flow analysis.
		In: Fourth National Conference of IBPSA-USA, New York City, New Your, August 11-13, 2010.

    \endverbatim

    The default model coefficients are
    \verbatim
        zetaCoeffs
        {
            Cmu         0.22;
            CmuKEps     0.09;
            C1          0.4;
            C2          0.65;
            CL          0.36;
			CT			6.0;
            Ceta        50;
            Ceps2       1.9;
            Ceps3       -0.33;
            sigmaEps    1.3;
            sigmaK      1;
            sigmaZeta   1.2;
        }
    \endverbatim



Note
    Ceta was set on 50 based on reference from 2010. Ceta=85 did not predicted a correct near-wall flow dynamics for the channel flow. For futher information see the paper from 2010.

See also
    Foam::RASModels::zetafBase
    Foam::RASModels::kEpsilon
    Foam::kLowReWallFunctionFvPatchScalarField
    Foam::epsilonLowReWallFunctionFvPatchScalarField
    Foam::zetaWallFunctionFvPatchScalarField
    Foam::fWallFunctionFvPatchScalarField

SourceFiles
    zetaf0.C

\*---------------------------------------------------------------------------*/

#ifndef zetaf0_H
#define zetaf0_H

#include "zetaf0Base.H"
#include "RASModel.H"
#include "eddyViscosity.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{
namespace RASModels
{

/*---------------------------------------------------------------------------*\
                             Class zetaf0 Declaration
\*---------------------------------------------------------------------------*/

template<class BasicTurbulenceModel>
class zetaf0
:
    public eddyViscosity<RASModel<BasicTurbulenceModel>>,
    public zetaf0Base
{


protected:

    // Protected data

        // Model coefficients

            dimensionedScalar Cmu_;
            dimensionedScalar CmuKEps_;
            dimensionedScalar C1_;
            dimensionedScalar C2_;
            dimensionedScalar CL_;
            dimensionedScalar Ceta_;
            dimensionedScalar CT_;
            dimensionedScalar Ceps2_;
            dimensionedScalar Ceps3_;
            dimensionedScalar sigmaK_;
            dimensionedScalar sigmaEps_;
            dimensionedScalar sigmaZeta_;


        // Fields

            //- Turbulence kinetic energy
            volScalarField k_;

            //- Turbulence dissipation
            volScalarField epsilon_;

            //- Turbulence stress normal to streamlines
            volScalarField zeta_;

            //- Damping function
            volScalarField f0_;


        // Bounding values

            dimensionedScalar zetaMin_;
            dimensionedScalar f0Min_;

            dimensionedScalar TScMin_;
            dimensionedScalar TInMin_;
            dimensionedScalar LScMin_;

    // Protected Member Functions

        virtual void correctNut();

        //- Return time scale, Ts
        tmp<volScalarField> Ts() const;

        //- Return length scale, Ls
        tmp<volScalarField> Ls() const;


public:

    typedef typename BasicTurbulenceModel::alphaField alphaField;
    typedef typename BasicTurbulenceModel::rhoField rhoField;
    typedef typename BasicTurbulenceModel::transportModel transportModel;


    //- Runtime type information
    TypeName("zetaf0");


    // Constructors

        //- Construct from components
        zetaf0
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


    //- Destructor
    virtual ~zetaf0()
    {}


    // Member Functions

        //- Read RASProperties dictionary
        virtual bool read();

        //- Return the effective diffusivity for k
        tmp<volScalarField> DkEff() const
        {
            return tmp<volScalarField>
            (
                new volScalarField
                (
                    "DkEff",
                    this->nut_/sigmaK_ + this->nu()
                )
            );
        }

        //- Return the effective diffusivity for epsilon
        tmp<volScalarField> DepsilonEff() const
        {
            return tmp<volScalarField>
            (
                new volScalarField
                (
                    "DepsilonEff",
                    this->nut_/sigmaEps_ + this->nu()
                )
            );
        }

        //- Return the effective diffusivity for zeta
        tmp<volScalarField> DzetaEff() const
        {
            return tmp<volScalarField>
            (
                new volScalarField
                (
                    "DzetaEff",
                    this->nut_/sigmaZeta_ + this->nu()
                )
            );
        }

        //- Return the turbulence kinetic energy
        virtual tmp<volScalarField> k() const
        {
            return k_;
        }

        //- Return the turbulence kinetic energy dissipation rate
        virtual tmp<volScalarField> epsilon() const
        {
            return epsilon_;
        }

        //- Return turbulence stress normal to streamlines
        virtual tmp<volScalarField> zeta() const
        {
            return zeta_;
        }

        //- Return the damping function
        virtual tmp<volScalarField> f0() const
        {
            return f0_;
        }

        //- Solve the turbulence equations and correct the turbulence viscosity
        virtual void correct();
};


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace RASModels
} // End namespace Foam

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#ifdef NoRepository
    #include "zetaf0.C"
#endif

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#endif

// ************************************************************************* //

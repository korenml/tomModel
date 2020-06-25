/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     | Website:  https://openfoam.org
    \\  /    A nd           | Copyright (C) 2011-2018 OpenFOAM Foundation
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
    Foam::RASModels::kEpsilonLRN

Description
    PANS k-epsilon Low Reynolds Number turbulence model for incompressible and compressible
    flows.

    Reference:
    \verbatim
        Standard model:
            Abe, K., & Kondoh, T. (1994).
            A new turbulence model for predicting fluid flow and heat transfer
			 in separating and reattaching flows-I. Flow field calculations.
            International Journal of Heat and Mass Transfer,Vol 37 139-151.

			
    \endverbatim

    The default model coefficients are
    \verbatim
        kEpsilonLRNCoeffs
        {
            Cmu         0.09;
            C1          1.5;
            C2          1.9;
            sigmak      1.4;
            sigmaEps    1.4;
        }
    \endverbatim

SourceFiles
    kEpsilonLRN.C

\*---------------------------------------------------------------------------*/

#ifndef kEpsilonLRN_H
#define kEpsilonLRN_H

#include "RASModel.H"
#include "eddyViscosity.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{
namespace RASModels
{

/*---------------------------------------------------------------------------*\
                           Class kEpsilonLRN Declaration
\*---------------------------------------------------------------------------*/

template<class BasicTurbulenceModel>
class kEpsilonLRN
:
    public eddyViscosity<RASModel<BasicTurbulenceModel>>
{
    // Private Member Functions

        // Disallow default bitwise copy construct and assignment
        kEpsilonLRN(const kEpsilonLRN&);
        void operator=(const kEpsilonLRN&);


protected:

    // Protected data

        // Model coefficients

            dimensionedScalar Cmu_;
            dimensionedScalar C1_;
            dimensionedScalar C2_;
            dimensionedScalar C3_;
            dimensionedScalar sigmak_;
            dimensionedScalar sigmaEps_;

        // Fields

			const volScalarField& y_;

            volScalarField k_;
            volScalarField epsilon_;





    // Protected Member Functions

        tmp<volScalarField> fMu();
        tmp<volScalarField> f2();


        virtual void correctNut();
        virtual tmp<fvScalarMatrix> kSource() const;
        virtual tmp<fvScalarMatrix> epsilonSource() const;


public:

    typedef typename BasicTurbulenceModel::alphaField alphaField;
    typedef typename BasicTurbulenceModel::rhoField rhoField;
    typedef typename BasicTurbulenceModel::transportModel transportModel;


    //- Runtime type information
    TypeName("kEpsilonLRN");


    // Constructors

        //- Construct from components
        kEpsilonLRN
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
    virtual ~kEpsilonLRN()
    {}


    // Member Functions

        //- Re-read model coefficients if they have changed
        virtual bool read();

        //- Return the effective diffusivity for k
        tmp<volScalarField> DkEff() const
        {
            return tmp<volScalarField>
            (
                new volScalarField
                (
                    "DkEff",
                    (this->nut_/(sigmak_) + this->nu())
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
                    (this->nut_/(sigmaEps_) + this->nu())
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

        //- Solve the turbulence equations and correct the turbulence viscosity
        virtual void correct();
};


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace RASModels
} // End namespace Foam

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#ifdef NoRepository
    #include "kEpsilonLRN.C"
#endif

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#endif

// ************************************************************************* //

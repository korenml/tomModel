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
    Foam::RASModels::kEpsilonPANSSSV

Description
    PANS k-epsilon turbulence model for incompressible and compressible
    flows including rapid distortion theory (RDT) based compression term.

    Reference:
    \verbatim
        Standard model:
            Launder, B. E., & Spalding, D. B. (1972).
            Lectures in mathematical models of turbulence.

            Launder, B. E., & Spalding, D. B. (1974).
            The numerical computation of turbulent flows.
            Computer methods in applied mechanics and engineering,
            3(2), 269-289.

			Girimaji, S. S. (2006).
            Partially-averaged Navier-Stokes method for turbulence: A Reynolds
            averaged Navier-Stokes to direct numerical simulation bridging
            method.
            Journal of Applied Mechanics, Vol 73, 413-421.
    \endverbatim

    The default model coefficients are
    \verbatim
        kEpsilonPANSSSVCoeffs
        {
            Cmu         0.09;
            C1          1.44;
            C2          1.92;
            C3          -0.33;
            sigmak      1.0;
            sigmaEps    1.3;
        }
    \endverbatim

SourceFiles
    kEpsilonPANSSSV.C

\*---------------------------------------------------------------------------*/

#ifndef kEpsilonPANSSSV_H
#define kEpsilonPANSSSV_H

#include "RASModel.H"
#include "eddyViscosity.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{
namespace RASModels
{

/*---------------------------------------------------------------------------*\
                           Class kEpsilonPANSSSV Declaration
\*---------------------------------------------------------------------------*/

template<class BasicTurbulenceModel>
class kEpsilonPANSSSV
:
    public eddyViscosity<RASModel<BasicTurbulenceModel>>
{
    // Private Member Functions

        // Disallow default bitwise copy construct and assignment
        kEpsilonPANSSSV(const kEpsilonPANSSSV&);
        void operator=(const kEpsilonPANSSSV&);


protected:

    // Protected data

        // Model coefficients

            dimensionedScalar Cmu_;
            dimensionedScalar C1_;
            dimensionedScalar C2_;
            dimensionedScalar C3_;
            dimensionedScalar sigmak_;
            dimensionedScalar sigmaEps_;
			dimensionedScalar fEpsilon_;
			dimensionedScalar fkLow_;

        // Fields

            autoPtr<Foam::LESdelta> delta_;
			volScalarField fk_;

            volScalarField k_;
            volScalarField kU_;
            volScalarField kSSV_;

            volScalarField C2U_;
            volScalarField epsilon_;
            volScalarField epsilonU_;


    // Protected Member Functions

        virtual void correctNut();
        virtual tmp<fvScalarMatrix> kSource() const;
        virtual tmp<fvScalarMatrix> epsilonSource() const;


public:

    typedef typename BasicTurbulenceModel::alphaField alphaField;
    typedef typename BasicTurbulenceModel::rhoField rhoField;
    typedef typename BasicTurbulenceModel::transportModel transportModel;


    //- Runtime type information
    TypeName("kEpsilonPANSSSV");


    // Constructors

        //- Construct from components
        kEpsilonPANSSSV
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
    virtual ~kEpsilonPANSSSV()
    {}


    // Member Functions

        //- Re-read model coefficients if they have changed
        virtual bool read();

        //- The filter width
        inline const volScalarField& delta() const
        {
            return delta_();
        }

        //- Return the effective diffusivity for k
        tmp<volScalarField> DkEff() const
        {
            return tmp<volScalarField>
            (
                new volScalarField
                (
                    "DkEff",
                    (this->nut_/(sigmak_*fk_*fk_/fEpsilon_) + this->nu())
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
                    (this->nut_/(sigmaEps_*fk_*fk_/fEpsilon_) + this->nu())
                )
            );
        }

        //- Return the turbulence kinetic energy
        virtual tmp<volScalarField> k() const
        {
            return k_;
        }

        //- Return the unresolved turbulence kinetic energy
        virtual tmp<volScalarField> kU() const
        {
            return kU_;
        }

        //- Return the resolved turbulence kinetic energy
        virtual tmp<volScalarField> kSSV() const
        {
            return kSSV_;
        }



        //- Return the turbulence kinetic energy dissipation rate
        virtual tmp<volScalarField> epsilon() const
        {
            return epsilon_;
        }

        //- Return the unresolved turbulence kinetic energy dissipation rate
        virtual tmp<volScalarField> epsilonU() const
        {
            return epsilonU_;
        }

        //- Solve the turbulence equations and correct the turbulence viscosity
        virtual void correct();
};


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace RASModels
} // End namespace Foam

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#ifdef NoRepository
    #include "kEpsilonPANSSSV.C"
#endif

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#endif

// ************************************************************************* //

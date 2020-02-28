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
    Foam::RASModels::kOmegaTNTPANS

Description
    PANS k-omega TNT turbulence model for
    incompressible and compressible flows.

    References:
    \verbatim
		Kok, J. C. (1999).
		Resolving the dependence on free-stream values for the k-omega turbulence model.
		National Aerospace Laboratory NLR Report.

		Basara, B. & Krajnović, S. & Grimaji, S. & Pavlovic, Z. (2011).
		Near-Wall Formulation of the Partially Averaged Navier-Stokes Turbulence Model.
		AIAA Journal Vol. 49, 2627-2636.

    \endverbatim

    The default model coefficients are
    \verbatim
        kOmegaTNTPANSdynCoeffs
        {
            Cmu         0.09;  // Equivalent to betaStar
            alpha       0.52;
            beta        0.072;
            alphak      0.5;
            alphaOmega  0.5;
        }
    \endverbatim

SourceFiles
    kOmegaTNTPANSdyn.C

\*---------------------------------------------------------------------------*/

#ifndef kOmegaTNTPANSdyn_H
#define kOmegaTNTPANSdyn_H

#include "RASModel.H"
#include "eddyViscosity.H"
#include "LESdelta.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{
namespace RASModels
{

/*---------------------------------------------------------------------------*\
                           Class kOmegaTNTPANSdyn Declaration
\*---------------------------------------------------------------------------*/

template<class BasicTurbulenceModel>
class kOmegaTNTPANSdyn
:
    public eddyViscosity<RASModel<BasicTurbulenceModel>>
{

protected:

    // Protected data

        // Model coefficients

            dimensionedScalar Cmu_;
            dimensionedScalar beta_;
            dimensionedScalar gamma_;
            dimensionedScalar alphaK_;
            dimensionedScalar alphaOmega_;
            dimensionedScalar alphaD_;

			dimensionedScalar fEpsilon_;

            autoPtr<Foam::LESdelta> delta_;

        // Fields
			volScalarField fk_;
			volScalarField fOmega_;
            volScalarField k_;
            volScalarField kU_;
            volScalarField omega_;
            volScalarField omegaU_;



    // Protected Member Functions

        virtual void correctNut();


public:

    typedef typename BasicTurbulenceModel::alphaField alphaField;
    typedef typename BasicTurbulenceModel::rhoField rhoField;
    typedef typename BasicTurbulenceModel::transportModel transportModel;


    //- Runtime type information
    TypeName("kOmegaTNTPANSdyn");


    // Constructors

        //- Construct from components
        kOmegaTNTPANSdyn
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
    virtual ~kOmegaTNTPANSdyn()
    {}


    // Member Functions

        //- Read RASProperties dictionary
        virtual bool read();

        //- Access function to filter width
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
                    (alphaK_*fOmega_/fk_)*this->nut_ + this->nu()
                )
            );
        }

        //- Return the effective diffusivity for omega
        tmp<volScalarField> DomegaEff() const
        {
            return tmp<volScalarField>
            (
                new volScalarField
                (
                    "DomegaEff",
                    (alphaOmega_*fOmega_/fk_)*this->nut_ + this->nu()
                )
            );
        }

        //- Return the ratio of unresolved to total turbulence kinetic energy
        virtual tmp<volScalarField> fk() const
        {
            return fk_;
        }
        //- Return the ratio of unresolved to total specific dissipation rate
        virtual tmp<volScalarField> fOmega() const
        {
            return fOmega_;
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

        //- Return the turbulence specific dissipation rate
        virtual tmp<volScalarField> omega() const
        {
            return omega_;
        }

        //- Return the unresolved turbulence specific dissipation rate
        virtual tmp<volScalarField> omegaU() const
        {
            return omegaU_;
        }


        //- Return the turbulence kinetic energy dissipation rate
        virtual tmp<volScalarField> epsilon() const
        {
            return tmp<volScalarField>
            (
                new volScalarField
                (
                    IOobject
                    (
                        "epsilon",
                        this->mesh_.time().timeName(),
                        this->mesh_
                    ),
                    Cmu_*k_*omega_,
                    omega_.boundaryField().types()
                )
            );
        }

        //- Return the unresolved turbulence kinetic energy dissipation rate
        virtual tmp<volScalarField> epsilonU() const
        {
            return tmp<volScalarField>
            (
                new volScalarField
                (
                    IOobject
                    (
                        "epsilonU",
                        this->mesh_.time().timeName(),
                        this->mesh_
                    ),
                    Cmu_*kU_*omegaU_,
                    omega_.boundaryField().types()
                )
            );
        }

        //- Solve the turbulence equations and correct the turbulence viscosity
        virtual void correct();
};


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace RASModels
} // End namespace Foam

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //
#ifdef NoRepository
    #include "kOmegaTNTPANSdyn.C"
#endif

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#endif

// ************************************************************************* //

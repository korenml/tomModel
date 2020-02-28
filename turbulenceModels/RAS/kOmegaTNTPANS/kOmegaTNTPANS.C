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

\*---------------------------------------------------------------------------*/

#include "kOmegaTNTPANS.H"
#include "fvOptions.H"
#include "bound.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{
namespace RASModels
{

// * * * * * * * * * * * * Protected Member Functions  * * * * * * * * * * * //

template<class BasicTurbulenceModel>
void kOmegaTNTPANS<BasicTurbulenceModel>::correctNut()
{
    this->nut_ = kU_/omegaU_;
    this->nut_.correctBoundaryConditions();
    fv::options::New(this->mesh_).correct(this->nut_);

    BasicTurbulenceModel::correctNut();
}


// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

template<class BasicTurbulenceModel>
kOmegaTNTPANS<BasicTurbulenceModel>::kOmegaTNTPANS
(
    const alphaField& alpha,
    const rhoField& rho,
    const volVectorField& U,
    const surfaceScalarField& alphaRhoPhi,
    const surfaceScalarField& phi,
    const transportModel& transport,
    const word& propertiesName,
    const word& type
)
:
    eddyViscosity<RASModel<BasicTurbulenceModel>>
    (
        type,
        alpha,
        rho,
        U,
        alphaRhoPhi,
        phi,
        transport,
        propertiesName
    ),

    Cmu_
    (
        dimensioned<scalar>::lookupOrAddToDict
        (
            "betaStar",
            this->coeffDict_,
            0.09
        )
    ),
    beta_
    (
        dimensioned<scalar>::lookupOrAddToDict
        (
            "beta",
            this->coeffDict_,
            0.075
        )
    ),
    gamma_
    (
        dimensioned<scalar>::lookupOrAddToDict
        (
            "gamma",
            this->coeffDict_,
            0.553
        )
    ),
    alphaK_
    (
        dimensioned<scalar>::lookupOrAddToDict
        (
            "alphaK",
            this->coeffDict_,
            0.667
        )
    ),
    alphaOmega_
    (
        dimensioned<scalar>::lookupOrAddToDict
        (
            "alphaOmega",
            this->coeffDict_,
            0.5
        )
    ),
	alphaD_
    (
        dimensioned<scalar>::lookupOrAddToDict
        (
            "alphaD",
            this->coeffDict_,
            0.5
        )
    ),
	fk_
    (
        dimensioned<scalar>::lookupOrAddToDict
        (
            "fk",
            this->coeffDict_,
            0.25
        )
    ),
	fEpsilon_
    (
        dimensioned<scalar>::lookupOrAddToDict
        (
            "fEpsilon",
            this->coeffDict_,
            1.0
        )
    ),
	fOmega_
    (
        dimensioned<scalar>::lookupOrAddToDict
        (
            "fOmega",
            this->coeffDict_,
            4.0
        )
    ),

    k_
    (
        IOobject
        (
            IOobject::groupName("k", alphaRhoPhi.group()),
            this->runTime_.timeName(),
            this->mesh_,
            IOobject::MUST_READ,
            IOobject::AUTO_WRITE
        ),
        this->mesh_
    ),
    kU_
    (
        IOobject
        (
            IOobject::groupName("kU", alphaRhoPhi.group()),
            this->runTime_.timeName(),
            this->mesh_,
            IOobject::NO_READ,
            IOobject::AUTO_WRITE
        ),
        k_*fk_,
		k_.boundaryField().types()
    ),
    omega_
    (
        IOobject
        (
            IOobject::groupName("omega", alphaRhoPhi.group()),
            this->runTime_.timeName(),
            this->mesh_,
            IOobject::MUST_READ,
            IOobject::AUTO_WRITE
        ),
        this->mesh_
    ),
    omegaU_
    (
        IOobject
        (
            IOobject::groupName("omegaU", alphaRhoPhi.group()),
            this->runTime_.timeName(),
            this->mesh_,
            IOobject::NO_READ,
            IOobject::AUTO_WRITE
        ),
		omega_*fOmega_,
		omega_.boundaryField().types()
    )

{
	
	fOmega_ = fEpsilon_/fk_;
    bound(k_, this->kMin_);
    bound(omega_, this->omegaMin_);

    bound(kU_, fk_*this->kMin_);
    bound(omegaU_, fOmega_*this->omegaMin_);

    if (type == typeName)
    {
        this->printCoeffs(type);
    }
	
}


// * * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * //

template<class BasicTurbulenceModel>
bool kOmegaTNTPANS<BasicTurbulenceModel>::read()
{
    if (eddyViscosity<RASModel<BasicTurbulenceModel>>::read())
    {
        Cmu_.readIfPresent(this->coeffDict());
        beta_.readIfPresent(this->coeffDict());
        gamma_.readIfPresent(this->coeffDict());
        alphaK_.readIfPresent(this->coeffDict());
        alphaOmega_.readIfPresent(this->coeffDict());
        alphaD_.readIfPresent(this->coeffDict());
		fk_.readIfPresent(this->coeffDict());
		fEpsilon_.readIfPresent(this->coeffDict());
		fOmega_.readIfPresent(this->coeffDict());

        return true;
    }
    else
    {
        return false;
    }
}


template<class BasicTurbulenceModel>
void kOmegaTNTPANS<BasicTurbulenceModel>::correct()
{
    if (!this->turbulence_)
    {
        return;
    }

    // Local references
    const alphaField& alpha = this->alpha_;
    const rhoField& rho = this->rho_;
    const surfaceScalarField& alphaRhoPhi = this->alphaRhoPhi_;
    const volVectorField& U = this->U_;
    volScalarField& nut = this->nut_;
    fv::options& fvOptions(fv::options::New(this->mesh_));

    eddyViscosity<RASModel<BasicTurbulenceModel>>::correct();

    volScalarField divU(fvc::div(fvc::absolute(this->phi(), U)));

    tmp<volTensorField> tgradU = fvc::grad(U);
    volScalarField G
    (
        this->GName(),
        nut*(tgradU() && dev(twoSymm(tgradU())))
    );
    tgradU.clear();

    // Update omega and G at the wall
    omegaU_.boundaryFieldRef().updateCoeffs();

    volScalarField CDkOmega
    (
        max(alphaD_*(fOmega_/fk_)*(fvc::grad(kU_) & fvc::grad(omegaU_))/omegaU_, dimensionedScalar("1.0e-10", dimless/sqr(dimTime), 1.0e-10))
    );

    // Turbulence specific dissipation rate equation
    tmp<fvScalarMatrix> omegaUEqn
    (
        fvm::ddt(alpha, rho, omegaU_)
      + fvm::div(alphaRhoPhi, omegaU_)
      - fvm::laplacian(alpha*rho*DomegaEff(), omegaU_)
     ==
        gamma_*alpha*rho*G*omegaU_/kU_
      - fvm::SuSp(((2.0/3.0)*gamma_)*alpha*rho*divU, omegaU_)
      - fvm::Sp(
				(gamma_*beta_ - gamma_*beta_/fOmega_ + beta_/fOmega_)*alpha*rho*omegaU_, omegaU_
				)
      + fvm::SuSp
        (
            alpha()*rho()*CDkOmega()/omegaU_(),
            omegaU_
        )
      + fvOptions(alpha, rho, omegaU_)
    );

    omegaUEqn.ref().relax();
    fvOptions.constrain(omegaUEqn.ref());
    omegaUEqn.ref().boundaryManipulate(omegaU_.boundaryFieldRef());
    solve(omegaUEqn);
    fvOptions.correct(omegaU_);
    bound(omegaU_, fOmega_*this->omegaMin_);


    // Turbulent kinetic energy equation
    tmp<fvScalarMatrix> kUEqn
    (
        fvm::ddt(alpha, rho, kU_)
      + fvm::div(alphaRhoPhi, kU_)
      - fvm::laplacian(alpha*rho*DkEff(), kU_)
     ==
        alpha*rho*G
      - fvm::SuSp((2.0/3.0)*alpha*rho*divU, kU_)
      - fvm::Sp(Cmu_*alpha*rho*omegaU_, kU_)
      + fvOptions(alpha, rho, kU_)
    );

    kUEqn.ref().relax();
    fvOptions.constrain(kUEqn.ref());
    solve(kUEqn);
    fvOptions.correct(kU_);
    bound(kU_, fk_*this->kMin_);

	this->k_ = kU_/fk_;
	this->k_.correctBoundaryConditions();

	this->omega_ = omegaU_/fOmega_;
	this->omega_.correctBoundaryConditions();

    correctNut();
}


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace RASModels
} // End namespace Foam

// ************************************************************************* //


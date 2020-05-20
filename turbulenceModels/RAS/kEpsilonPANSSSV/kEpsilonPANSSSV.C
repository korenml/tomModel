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

#include "kEpsilonPANSSSV.H"
#include "fvOptions.H"
#include "bound.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{
namespace RASModels
{

// * * * * * * * * * * * * Protected Member Functions  * * * * * * * * * * * //

template<class BasicTurbulenceModel>
void kEpsilonPANSSSV<BasicTurbulenceModel>::correctNut()
{
    this->nut_ = Cmu_*sqr(kU_)/epsilonU_;
    this->nut_.correctBoundaryConditions();
    fv::options::New(this->mesh_).correct(this->nut_);

    BasicTurbulenceModel::correctNut();
}


template<class BasicTurbulenceModel>
tmp<fvScalarMatrix> kEpsilonPANSSSV<BasicTurbulenceModel>::kSource() const
{
    return tmp<fvScalarMatrix>
    (
        new fvScalarMatrix
        (
            kU_,
            dimVolume*this->rho_.dimensions()*kU_.dimensions()
            /dimTime
        )
    );
}


template<class BasicTurbulenceModel>
tmp<fvScalarMatrix> kEpsilonPANSSSV<BasicTurbulenceModel>::epsilonSource() const
{
    return tmp<fvScalarMatrix>
    (
        new fvScalarMatrix
        (
            epsilonU_,
            dimVolume*this->rho_.dimensions()*epsilonU_.dimensions()
            /dimTime
        )
    );
}


// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

template<class BasicTurbulenceModel>
kEpsilonPANSSSV<BasicTurbulenceModel>::kEpsilonPANSSSV
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
            "Cmu",
            this->coeffDict_,
            0.09
        )
    ),
    C1_
    (
        dimensioned<scalar>::lookupOrAddToDict
        (
            "C1",
            this->coeffDict_,
            1.44
        )
    ),
    C2_
    (
        dimensioned<scalar>::lookupOrAddToDict
        (
            "C2",
            this->coeffDict_,
            1.92
        )
    ),
    C3_
    (
        dimensioned<scalar>::lookupOrAddToDict
        (
            "C3",
            this->coeffDict_,
            0
        )
    ),

    sigmak_
    (
        dimensioned<scalar>::lookupOrAddToDict
        (
            "sigmak",
            this->coeffDict_,
            1.0
        )
    ),
    sigmaEps_
    (
        dimensioned<scalar>::lookupOrAddToDict
        (
            "sigmaEps",
            this->coeffDict_,
            1.3
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

    fkLow_
    (
        dimensioned<scalar>::lookupOrAddToDict
        (
            "fkLow",
            this->coeffDict_,
            0.1
        )
    ),

    delta_
    (
        LESdelta::New
        (
            IOobject::groupName("delta", U.group()),
            *this,
            this->coeffDict_
        )
    ),

    fk_
    (
        IOobject
        (
            IOobject::groupName("fk", alphaRhoPhi.group()),
            this->runTime_.timeName(),
            this->mesh_,
            IOobject::MUST_READ,
            IOobject::AUTO_WRITE
        ),
        this->mesh_,
        dimensionedScalar("low", fkLow_)
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

    kSSV_
    (
        IOobject
        (
            IOobject::groupName("kSSV", alphaRhoPhi.group()),
            this->runTime_.timeName(),
            this->mesh_,
            IOobject::NO_READ,
            IOobject::AUTO_WRITE
        ),
        k_*(1.0-fk_),
		k_.boundaryField().types()
    ),

    C2U_
    (
        IOobject
        (
            "C2U",
            this->runTime_.timeName(),
            this->mesh_
        ),
		C1_ + (fk_/fEpsilon_)*(C2_ - C1_)
    ),

    epsilon_
    (
        IOobject
        (
            IOobject::groupName("epsilon", alphaRhoPhi.group()),
            this->runTime_.timeName(),
            this->mesh_,
            IOobject::MUST_READ,
            IOobject::AUTO_WRITE
        ),
        this->mesh_
    ),

    epsilonU_
    (
        IOobject
        (
            IOobject::groupName("epsilonU", alphaRhoPhi.group()),
            this->runTime_.timeName(),
            this->mesh_,
            IOobject::NO_READ,
            IOobject::AUTO_WRITE
        ),
        epsilon_*fEpsilon_,
		epsilon_.boundaryField().types()
    )

{

    bound(k_, this->kMin_);
    bound(epsilon_, this->epsilonMin_);


    bound(kU_, min(fk_)*this->kMin_);
    bound(kSSV_, (1.0-min(fk_))*this->kMin_);
    bound(epsilonU_, fEpsilon_*this->epsilonMin_);

    if (type == typeName)
    {
        this->printCoeffs(type);
    }
}


// * * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * //

template<class BasicTurbulenceModel>
bool kEpsilonPANSSSV<BasicTurbulenceModel>::read()
{
    if (eddyViscosity<RASModel<BasicTurbulenceModel>>::read())
    {
        Cmu_.readIfPresent(this->coeffDict());
        C1_.readIfPresent(this->coeffDict());
        C2_.readIfPresent(this->coeffDict());
        C3_.readIfPresent(this->coeffDict());
        sigmak_.readIfPresent(this->coeffDict());
        sigmaEps_.readIfPresent(this->coeffDict());
		fEpsilon_.readIfPresent(this->coeffDict());
		fkLow_.readIfPresent(this->coeffDict());

        return true;
    }
    else
    {
        return false;
    }
}


template<class BasicTurbulenceModel>
void kEpsilonPANSSSV<BasicTurbulenceModel>::correct()
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

    volScalarField::Internal divU
    (
        fvc::div(fvc::absolute(this->phi(), U))().v()
    );
	


    tmp<volTensorField> tgradU = fvc::grad(U);
    volScalarField::Internal G
    (
        this->GName(),
        nut.v()*(dev(twoSymm(tgradU().v())) && tgradU().v())
    );
    tgradU.clear();

	volScalarField::Internal fp_
	(
		((1.0-epsilonU_/G)/fk_+epsilonU_/(fEpsilon_*G))
	);

    // Update epsilonU and G at the wall
    epsilonU_.boundaryFieldRef().updateCoeffs();

    // Dissipation equation
    tmp<fvScalarMatrix> epsUEqn
    (
        fvm::ddt(alpha, rho, epsilonU_)
      + fvm::div(alphaRhoPhi, epsilonU_)
      - fvm::laplacian(alpha*rho*DepsilonEff(), epsilonU_)
     ==
        C1_*alpha()*rho()*G*epsilonU_()/kU_()
      - fvm::SuSp(((2.0/3.0)*C1_ - C3_)*alpha()*rho()*divU, epsilonU_)
//      - fvm::Sp(alpha()*rho()*epsilonU_()/kU_(), epsilonU_)
      - fvm::Sp(C2U_*alpha()*rho()*epsilonU_()/kU_(), epsilonU_)
      + epsilonSource()
      + fvOptions(alpha, rho, epsilonU_)
    );

    epsUEqn.ref().relax();
    fvOptions.constrain(epsUEqn.ref());
    epsUEqn.ref().boundaryManipulate(epsilonU_.boundaryFieldRef());
    solve(epsUEqn);
    fvOptions.correct(epsilonU_);
    bound(epsilonU_, fEpsilon_*this->epsilonMin_);

    // Unresolved turbulent kinetic energy equation
    tmp<fvScalarMatrix> kUEqn
    (
        fvm::ddt(alpha, rho, kU_)
      + fvm::div(alphaRhoPhi, kU_)
      - fvm::laplacian(alpha*rho*DkEff(), kU_)
     ==
        alpha()*rho()*G
      - fvm::SuSp((2.0/3.0)*alpha()*rho()*divU, kU_)
      - fvm::Sp(alpha()*rho()*epsilonU_()/kU_(), kU_)
      + kSource()
      + fvOptions(alpha, rho, kU_)
    );

    kUEqn.ref().relax();
    fvOptions.constrain(kUEqn.ref());
    solve(kUEqn);
    fvOptions.correct(kU_);
    bound(kU_, min(fk_)*this->kMin_);

    // scale supplying turbulent kinetic energy equation
    tmp<fvScalarMatrix> kSSVEqn
    (
        fvm::ddt(alpha, rho, kSSV_)
      + fvm::div(alphaRhoPhi, kSSV_)
      - fvm::laplacian(alpha*rho*DkEff(), kSSV_)
     ==
        fp_*alpha()*rho()*G
	  -	fk_*fp_*alpha()*rho()*G
      - fvm::SuSp((2.0/3.0)*alpha()*rho()*divU, kSSV_)
      - fvm::Sp(1.0*alpha()*rho()*epsilonU_()/(fEpsilon_*kSSV_()), kSSV_)
	  +	fvm::Sp(fk_*alpha()*rho()*epsilonU_()/(fEpsilon_*kSSV_()), kSSV_)
      + fvOptions(alpha, rho, kSSV_)
    );

    kSSVEqn.ref().relax();
    fvOptions.constrain(kSSVEqn.ref());
    solve(kSSVEqn);
    fvOptions.correct(kSSV_);
    bound(kSSV_, (1.0-min(fk_))*this->kMin_);

	this->k_ = kU_ + kSSV_;
	this->k_.correctBoundaryConditions();

	this->epsilon_ = epsilonU_/fEpsilon_;
	this->epsilon_.correctBoundaryConditions();

    correctNut();


	volScalarField Lt(pow(k_,1.5)/(epsilon_));

	//dynamic fk_
	fk_.primitiveFieldRef() = min
	(
							max((1.0/(sqrt(Cmu_))*pow(delta()/Lt,2.0/3.0)), fkLow_),
							1.0
	);

	C2U_ = C1_ + (fk_/fEpsilon_)*(C2_ - C1_);
}


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace RASModels
} // End namespace Foam

// ************************************************************************* //

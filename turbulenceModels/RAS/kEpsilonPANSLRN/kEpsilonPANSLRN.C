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

#include "kEpsilonPANSLRN.H"
#include "fvOptions.H"
#include "bound.H"
#include "wallDist.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{
namespace RASModels
{

// * * * * * * * * * * * * Protected Member Functions  * * * * * * * * * * * //

template<class BasicTurbulenceModel>
tmp<volScalarField> kEpsilonPANSLRN<BasicTurbulenceModel>::fMu() const
{
//yStar = (pow(epsilonU_*this->nu(), 0.25)*y_/this->nu())
//Ret = (sqr(kU_)/(this->nu()*epsilonU_))
//	tmp<volScalarField> yStar1 = pow(epsilonU_*this->nu(), 0.25)*y_/this->nu();
//	tmp<volScalarField> Ret1 = sqr(kU_)/(this->nu()*epsilonU_);

    return
		sqr(scalar(1) - exp(-(pow(epsilonU_*this->nu(), 0.25)*y_/this->nu())/14.0))*(scalar(1) + (5.0/pow((sqr(kU_)/(this->nu()*epsilonU_)), 3.0/4.0))*exp(-sqr((sqr(kU_)/(this->nu()*epsilonU_))/200.0)));


		//sqr(scalar(1) - exp(-yStar1/14.0))*(scalar(1) + (5.0/pow(Ret1, 3.0/4.0))*exp(-sqr(Ret1/200.0)));
}


template<class BasicTurbulenceModel>
tmp<volScalarField> kEpsilonPANSLRN<BasicTurbulenceModel>::f2() const
{
	tmp<volScalarField> yStar = pow(epsilonU_*this->nu(), 0.25)*y_/this->nu();
	tmp<volScalarField> Ret = sqr(kU_)/(this->nu()*epsilonU_);

    return 
		sqr(scalar(1) - exp(-yStar/14.0))*(scalar(1) - 0.3*exp(-sqr(Ret/6.5)));
}

template<class BasicTurbulenceModel>
void kEpsilonPANSLRN<BasicTurbulenceModel>::correctNut()
{
    this->nut_ = Cmu_*fMu()*sqr(kU_)/epsilonU_;
    this->nut_.correctBoundaryConditions();
    fv::options::New(this->mesh_).correct(this->nut_);

    BasicTurbulenceModel::correctNut();
}


template<class BasicTurbulenceModel>
tmp<fvScalarMatrix> kEpsilonPANSLRN<BasicTurbulenceModel>::kSource() const
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
tmp<fvScalarMatrix> kEpsilonPANSLRN<BasicTurbulenceModel>::epsilonSource() const
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
kEpsilonPANSLRN<BasicTurbulenceModel>::kEpsilonPANSLRN
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
            1.5
        )
    ),
    C2_
    (
        dimensioned<scalar>::lookupOrAddToDict
        (
            "C2",
            this->coeffDict_,
            1.9
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
            1.4
        )
    ),
    sigmaEps_
    (
        dimensioned<scalar>::lookupOrAddToDict
        (
            "sigmaEps",
            this->coeffDict_,
            1.4
        )
    ),

    fk_
    (
        dimensioned<scalar>::lookupOrAddToDict
        (
            "fk",
            this->coeffDict_,
            0.5
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

	y_(wallDist::New(this->mesh_).y()),

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


    bound(kU_, fk_*this->kMin_);
    bound(epsilonU_, fEpsilon_*this->epsilonMin_);

    if (type == typeName)
    {
        this->printCoeffs(type);
    }
}


// * * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * //

template<class BasicTurbulenceModel>
bool kEpsilonPANSLRN<BasicTurbulenceModel>::read()
{
    if (eddyViscosity<RASModel<BasicTurbulenceModel>>::read())
    {
        Cmu_.readIfPresent(this->coeffDict());
        C1_.readIfPresent(this->coeffDict());
        C2_.readIfPresent(this->coeffDict());
        C3_.readIfPresent(this->coeffDict());
        sigmak_.readIfPresent(this->coeffDict());
        sigmaEps_.readIfPresent(this->coeffDict());
		fk_.readIfPresent(this->coeffDict());
		fEpsilon_.readIfPresent(this->coeffDict());

        return true;
    }
    else
    {
        return false;
    }
}


template<class BasicTurbulenceModel>
void kEpsilonPANSLRN<BasicTurbulenceModel>::correct()
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
	const volScalarField fMu(this->fMu());
	const volScalarField f2(this->f2());


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
      - fvm::Sp((C1_+(fk_/fEpsilon_)*(f2()*C2_-C1_))*alpha()*rho()*epsilonU_()/kU_(), epsilonU_)
      + epsilonSource()
      + fvOptions(alpha, rho, epsilonU_)
    );

    epsUEqn.ref().relax();
    fvOptions.constrain(epsUEqn.ref());
    epsUEqn.ref().boundaryManipulate(epsilonU_.boundaryFieldRef());
    solve(epsUEqn);
    fvOptions.correct(epsilonU_);
    bound(epsilonU_, fEpsilon_*this->epsilonMin_);

    // Turbulent kinetic energy equation
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
    bound(kU_, fk_*this->kMin_);

	k_ = kU_/fk_;
	k_.correctBoundaryConditions();

	epsilon_ = epsilonU_/fEpsilon_;
	epsilon_.correctBoundaryConditions();


    correctNut();
}


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace RASModels
} // End namespace Foam

// ************************************************************************* //

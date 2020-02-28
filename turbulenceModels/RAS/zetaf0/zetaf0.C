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

\*---------------------------------------------------------------------------*/

#include "zetaf0.H"
#include "fvOptions.H"
#include "bound.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{
namespace RASModels
{

// * * * * * * * * * * * * * Private Member Functions  * * * * * * * * * * * //

template<class BasicTurbulenceModel>
tmp<volScalarField> zetaf0<BasicTurbulenceModel>::Ts() const
{
    return max(min(k_/(epsilon_ + this->epsilonMin_),
			CT_/(TInMin_ + zeta_*Cmu_*sqrt(6*2*magSqr(dev(symm(fvc::grad(this->U_))))))),
			6.0*sqrt(this->nu()/(epsilon_ + this->epsilonMin_)));
	//return max(k_/epsilon_, 6.0*sqrt(this->nu()/epsilon_));
}


template<class BasicTurbulenceModel>
tmp<volScalarField> zetaf0<BasicTurbulenceModel>::Ls() const
{
    return CL_*max(min(pow(k_, 1.5)/(epsilon_ + this->epsilonMin_),
			sqrt(k_)/(TInMin_ + zeta_*Cmu_*sqrt(6*2*magSqr(dev(symm(fvc::grad(this->U_))))))),
			Ceta_*pow025(pow3(this->nu())/(epsilon_ + this->epsilonMin_)));
	//return CL_*max(pow(k_, 1.5)/epsilon_, Ceta_*pow025(pow3(this->nu())/epsilon_));
}


template<class BasicTurbulenceModel>
void zetaf0<BasicTurbulenceModel>::correctNut()
{
	this->nut_ =  Cmu_*zeta_*k_*Ts();
    this->nut_.correctBoundaryConditions();
    fv::options::New(this->mesh_).correct(this->nut_);

    BasicTurbulenceModel::correctNut();
}


// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

template<class BasicTurbulenceModel>
zetaf0<BasicTurbulenceModel>::zetaf0
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
	zetaf0Base(),

    Cmu_
    (
        dimensioned<scalar>::lookupOrAddToDict
        (
            "Cmu",
            this->coeffDict_,
            0.22
        )
    ),
    CmuKEps_
    (
        dimensioned<scalar>::lookupOrAddToDict
        (
            "CmuKEps",
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
            1.4
        )
    ),
    C2_
    (
        dimensioned<scalar>::lookupOrAddToDict
        (
            "C2",
            this->coeffDict_,
            0.65
        )
    ),
    CL_
    (
        dimensioned<scalar>::lookupOrAddToDict
        (
            "CL",
            this->coeffDict_,
            0.36
        )
    ),
    Ceta_
    (
        dimensioned<scalar>::lookupOrAddToDict
        (
            "Ceta",
            this->coeffDict_,
            50.0
        )
    ),
    CT_
    (
        dimensioned<scalar>::lookupOrAddToDict
        (
            "CT",
            this->coeffDict_,
            0.6
        )
    ),
    Ceps2_
    (
        dimensioned<scalar>::lookupOrAddToDict
        (
            "Ceps2",
            this->coeffDict_,
            1.9
        )
    ),
    Ceps3_
    (
        dimensioned<scalar>::lookupOrAddToDict
        (
            "Ceps3",
            this->coeffDict_,
            -0.33
        )
    ),

    sigmaK_
    (
        dimensioned<scalar>::lookupOrAddToDict
        (
            "sigmaK",
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

    sigmaZeta_
    (
        dimensioned<scalar>::lookupOrAddToDict
        (
            "sigmaZeta",
            this->coeffDict_,
            1.2
        )
    ),

    k_
    (
        IOobject
        (
            IOobject::groupName("k", U.group()),
            this->runTime_.timeName(),
            this->mesh_,
            IOobject::MUST_READ,
            IOobject::AUTO_WRITE
        ),
        this->mesh_
    ),

    epsilon_
    (
        IOobject
        (
            IOobject::groupName("epsilon", U.group()),
            this->runTime_.timeName(),
            this->mesh_,
            IOobject::MUST_READ,
            IOobject::AUTO_WRITE
        ),
        this->mesh_
    ),

    zeta_
    (
        IOobject
        (
            IOobject::groupName("zeta", U.group()),
            this->runTime_.timeName(),
            this->mesh_,
            IOobject::MUST_READ,
            IOobject::AUTO_WRITE
        ),
        this->mesh_
    ),

    f0_
    (
        IOobject
        (
            IOobject::groupName("f0", U.group()),
            this->runTime_.timeName(),
            this->mesh_,
            IOobject::MUST_READ,
            IOobject::AUTO_WRITE
        ),
        this->mesh_
    ),

    zetaMin_(dimensionedScalar("zetaMin", dimless, SMALL)),
    f0Min_(dimensionedScalar("f0Min", dimless/dimTime, SMALL)),
    TScMin_(dimensionedScalar("TScMin", dimTime, SMALL)),
    TInMin_(dimensionedScalar("TInMin", dimless/dimTime, SMALL)),
    LScMin_(dimensionedScalar("LScMin", dimLength, SMALL))
{
    bound(k_, this->kMin_);
    bound(epsilon_, this->epsilonMin_);
    bound(zeta_, zetaMin_);
    bound(f0_, f0Min_);

    if (type == typeName)
    {
        this->printCoeffs(type);
    }
}


// * * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * //

template<class BasicTurbulenceModel>
bool zetaf0<BasicTurbulenceModel>::read()
{
    if (eddyViscosity<RASModel<BasicTurbulenceModel>>::read())
    {
        Cmu_.readIfPresent(this->coeffDict());
        CmuKEps_.readIfPresent(this->coeffDict());
        C1_.readIfPresent(this->coeffDict());
        C2_.readIfPresent(this->coeffDict());
        CL_.readIfPresent(this->coeffDict());
		Ceta_.readIfPresent(this->coeffDict());
		CT_.readIfPresent(this->coeffDict());
        Ceps2_.readIfPresent(this->coeffDict());
        Ceps3_.readIfPresent(this->coeffDict());
        sigmaK_.readIfPresent(this->coeffDict());
        sigmaEps_.readIfPresent(this->coeffDict());
        sigmaZeta_.readIfPresent(this->coeffDict());

        return true;
    }
    else
    {
        return false;
    }
}


template<class BasicTurbulenceModel>
void zetaf0<BasicTurbulenceModel>::correct()
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

    const volTensorField gradU(fvc::grad(U));
    const volScalarField S2(2*magSqr(dev(symm(gradU))));

    const volScalarField G(this->GName(), nut*S2);
    const volScalarField Ts(this->Ts() + TScMin_);
    const volScalarField L2(type() + ":L2", sqr(this->Ls()+LScMin_));

    const volScalarField Ceps1
    (
        "Ceps1",
        1.4*(1.0 + 0.012/(zeta_ + zetaMin_))
    );

    // Update epsilon (and possibly G) at the wall
    epsilon_.boundaryFieldRef().updateCoeffs();

    // Dissipation equation
    tmp<fvScalarMatrix> epsEqn
    (
        fvm::ddt(alpha, rho, epsilon_)
      + fvm::div(alphaRhoPhi, epsilon_)
      - fvm::laplacian(alpha*rho*DepsilonEff(), epsilon_)
     ==
        Ceps1*alpha*rho*G/Ts
      - fvm::Sp(Ceps2_*alpha*rho/Ts, epsilon_)
      + fvOptions(alpha, rho, epsilon_)
    );

    epsEqn.ref().relax();
    fvOptions.constrain(epsEqn.ref());
    epsEqn.ref().boundaryManipulate(epsilon_.boundaryFieldRef());
    solve(epsEqn);
    fvOptions.correct(epsilon_);
    bound(epsilon_, this->epsilonMin_);


    // Turbulent kinetic energy equation
    tmp<fvScalarMatrix> kEqn
    (
        fvm::ddt(alpha, rho, k_)
      + fvm::div(alphaRhoPhi, k_)
      - fvm::laplacian(alpha*rho*DkEff(), k_)
     ==
        alpha*rho*G
//      - fvm::Sp(alpha*rho*epsilon_/k_, k_)
	  - fvm::Sp(alpha*rho/Ts, k_)
      + fvOptions(alpha, rho, k_)
    );

    kEqn.ref().relax();
    fvOptions.constrain(kEqn.ref());
    solve(kEqn);
    fvOptions.correct(k_);
    bound(k_, this->kMin_);


    // Relaxation function equation
    tmp<fvScalarMatrix> f0Eqn
    (
      - fvm::laplacian(f0_)
     ==
      - fvm::Sp(1.0/L2, f0_)
      - (1.0/(Ts*L2))*(C1_ - 1.0 + C2_*G/(epsilon_ + this->epsilonMin_))*(zeta_ - (2.0/3.0))
	  - zeta_/(Ts*L2)
    );

    f0Eqn.ref().relax();
    fvOptions.constrain(f0Eqn.ref());
    solve(f0Eqn);
    fvOptions.correct(f0_);
    bound(f0_, f0Min_);


    // Turbulence stress normal to streamlines equation
    tmp<fvScalarMatrix> zetaEqn
    (
        fvm::ddt(alpha, rho, zeta_)
      + fvm::div(alphaRhoPhi, zeta_)
      - fvm::laplacian(alpha*rho*DzetaEff(), zeta_)
      ==
      - fvm::Sp(alpha*rho*(G+epsilon_)/k_,zeta_)
	  +	alpha*rho*f0_
      + fvOptions(alpha, rho, zeta_)
    );

    zetaEqn.ref().relax();
    fvOptions.constrain(zetaEqn.ref());
    solve(zetaEqn);
    fvOptions.correct(zeta_);
    bound(zeta_, zetaMin_);

    correctNut();
}


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace RASModels
} // End namespace Foam

// ************************************************************************* //

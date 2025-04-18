/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     |
    \\  /    A nd           | Copyright (C) 2011-2016 OpenFOAM Foundation
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


    Org Author:
	Ajit Kumar
	Shiv Nadar University

	Current Author:
	TAY WB
	TL@NUS 

\*---------------------------------------------------------------------------*/

#include "moving_sixDoFDynamicMotion.H"
#include "addToRunTimeSelectionTable.H"

#include "polyMesh.H"
#include "pointPatchDist.H"
#include "pointConstraints.H"
#include "uniformDimensionedFields.H"
#include "forces.H"
#include "mathematicalConstants.H"

// * * * * * * * * * * * * * * Static Data Members * * * * * * * * * * * * * //

namespace Foam
{
namespace solidBodyMotionFunctions
{
    defineTypeNameAndDebug(moving_sixDoFDynamicMotion, 0);
    addToRunTimeSelectionTable
    (
        solidBodyMotionFunction,
        moving_sixDoFDynamicMotion,
        dictionary
    );
}
}


// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

Foam::solidBodyMotionFunctions::moving_sixDoFDynamicMotion::moving_sixDoFDynamicMotion
(
    const dictionary& SBMFCoeffs,
    const Time& runTime
)
:
    solidBodyMotionFunction(SBMFCoeffs, runTime),
	dict_(SBMFCoeffs),
	//motion_(SBMFCoeffs,SBMFCoeffs),
	//motion_(time_.time()),
	motion_(SBMFCoeffs,SBMFCoeffs,time_.time()),
	curTimeIndex_(-1),
	patches_(wordRes(SBMFCoeffs.lookup("patches"))),
	rhoInf_(SBMFCoeffs_.getOrDefault<scalar>("rhoInf", 1.0)),
	rhoName_(SBMFCoeffs.lookupOrDefault<word>("rho", "rho")),
	
	firstRun_(true),
	fileWriteInterval_(runTime.controlDict().get<scalar>("writeInterval")),
	mass_(SBMFCoeffs_.get<scalar>("mass")),
	centreOfMass_(SBMFCoeffs_.get<vector>("centreOfMass")),
    momentOfInertia_(SBMFCoeffs_.get<vector>("momentOfInertia")),
    accelerationRelaxation_(SBMFCoeffs_.getOrDefault<scalar>("accelerationRelaxation", 1.0)),
	accelerationDamping_(SBMFCoeffs_.getOrDefault<scalar>("accelerationDamping", 1.0)),
	g_(SBMFCoeffs_.get<vector>("g")),
	gamma_new_(SBMFCoeffs_.getOrDefault<scalar>("gamma_new", 0.5)),
	beta_new_(SBMFCoeffs_.getOrDefault<scalar>("beta_new", 0.25)),
	x_displ_old(0),
	vel_x_old(0),
	acc_x_old(0),
	y_displ_old(0),
	vel_y_old(0),
	acc_y_old(0),
	//y_displ_ramp(0),
	//vel_y_ramp(0),
	//acc_y_ramp(0),
	saved_time_(-1),
	saved_no(1)
	
{
    read(SBMFCoeffs);

//    g_ = SBMFCoeffs.lookupObject<uniformDimensionedVectorField>("g");

}


// * * * * * * * * * * * * * * * * Destructor  * * * * * * * * * * * * * * * //

Foam::solidBodyMotionFunctions::moving_sixDoFDynamicMotion::~moving_sixDoFDynamicMotion()
{}


// * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * * //

Foam::septernion
Foam::solidBodyMotionFunctions::moving_sixDoFDynamicMotion::transformation() const
{
    
    /*
    Ajit: Apr 9
    This class is called by motion solver through this function.
    The aim should be to 
    - integrate force on the solid body patch,
    - calculate the six dof motion 
    - and return the motion as a septernion    
	- added moving grid functionality
    */
	
	scalar tmp_value;
	scalar delta_time;
	scalar acc_x;
	scalar vel_x;
	scalar x_displ;
	scalar delta_x;
	scalar acc_y;
	scalar vel_y;
	scalar y_displ;
	scalar delta_y;
    scalar t = time_.value();
	//scalar write_interval = time_.writeInterval();
	vector displacement;

	//Store the motion state at the beginning of the time-step
    bool firstIter = false;
    if (curTimeIndex_ != time_.time().timeIndex())
    {
        motion_.newTime();
        curTimeIndex_ = time_.time().timeIndex();
        firstIter = true;
    }

    // const scalar ramp = min(max((this->db().time().value() - 5)/10, 0), 1);
    const scalar ramp = 1.0;


    // add force
    dictionary forcesDict;
    forcesDict.add("type", functionObjects::forces::typeName);
    forcesDict.add("patches", patches_);

    forcesDict.add("rhoInf", rhoInf_);
    forcesDict.add("rho", rhoName_);
    forcesDict.add("CofR", motion_.centreOfRotation());


    functionObjects::forces f("forces", time_, forcesDict);

    f.calcForcesMoments();

    dimensionedVector g("g", dimAcceleration, Zero);
    if (dict_.found("g"))
    {
        dict_.lookup("g") >> g;
    }

    Info << "g = " << g << endl;

    // TODO: save the motion state before update
    
    // current center of mass
    vector pi = motion_.state().centreOfRotation ();

    // current orientation
    tensor Qi = motion_.state().Q();


    motion_.update(
        firstIter,
        ramp * (f.forceEff() + motion_.mass() * g.value()),
        ramp * (f.momentEff() + motion_.mass() * (motion_.momentArm() ^ g.value())),
        time_.deltaTValue(),
        time_.deltaT0Value());


    // TODO: extract the motion state after  the update
    // new center of mass
    vector pn = motion_.state().centreOfRotation ();

    // new orientation
    tensor Qn = motion_.state().Q();

    // TODO: create a septernion by comparting the before and after state
    // return that value

    // check
    // R  = new * inv(old)
    // call the constructor of rotation quaternion using the rotation tensor
    //R Qi vl = Qo vl

    //quaternion R(Qn & inv(Qi));

    //septernion TR(pn-pi,R );


    scalar t = time_.value();

    // // // Translation of centre of gravity with constant velocity
    const vector displacement = velocity_*t;
    // //const vector displacement = pn-pi;

    quaternion R(Qn);
 // TODO: add rotation later
    septernion TR(septernion(-pn)*R);


    // septernion TR(septernion(-displacement)*R);

    DebugInFunction << "Time = " << t << " transformation: " << TR << endl;

    
    // ofstream myfile;
    // myfile.open ("path.csv",std::ios::app);
    // myfile << t << '\t' << f.forceEff()[0] << '\t' << f.forceEff()[1] << '\n';
    // myfile.close();

    return TR;
}


bool Foam::solidBodyMotionFunctions::moving_sixDoFDynamicMotion::read
(
    const dictionary& SBMFCoeffs
)
{
    solidBodyMotionFunction::read(SBMFCoeffs);

    SBMFCoeffs_.lookup("velocity") >> velocity_;

    return true;
}


// ************************************************************************* //

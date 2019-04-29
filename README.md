## Coupling solidBodyMotionFunction with sixDOF solver in OpenFOAM-v1812

### Simple-pendulum simulation in OpenFoam

<center>
<table>
    <tr>
        <td>
          <img src="https://github.com/krajit/sixDoFDynamicMotion/blob/master/images/pendulum_U.gif?raw=true" alt="simple pendulum openfoam" height="300rem" width="300rem">
        </td>        
        <td>
          <img src="https://github.com/krajit/sixDoFDynamicMotion/blob/master/images/pendulum_meshMotion.gif?raw=true" alt="simple pendulum openfoam" height="300rem" width="300rem">
        </td>        
</tr>
</table>

### A projectile motion simulation in OpenFoam

<table>
    <tr>
        <td>
          <img src="https://github.com/krajit/sixDoFDynamicMotion/blob/master/images/projectileVelocity.gif" alt="projectile simulation in openfoam" height="300rem" width="500rem">
        </td>        
        <td>
          <img src="https://github.com/krajit/sixDoFDynamicMotion/blob/master/images/projectileMesh.gif?raw=true" alt="projectile simulation in openfoam" height="300rem" width="500rem">
        </td>        
</tr>
</table>


OpenFOAM has this class called `solidBodyMotionFunction` which is very convenient for solving flow around moving bodies. However, all the option under this class like, linearMotion, rotatingMotion, etc all seem to solve flow around presribed motion. I needed to solve for floating body whose motion needs to be solved and the motion is coupled with flow around it. OpenFoam also has this `sixDoFRigidBodyMotion` class, which allows intgration for pressure and shear force around a patch, and solve for the resulting motion of the rigid body associated with this patch. 

![](https://github.com/krajit/sixDoFDynamicMotion/blob/master/images/inheritence.png?raw=true)


This repository is is the source code for  `sixDoFDynamicMotion`, is a child of `solidBodyMotionFunction`, which as private member `motion_` of the type `sixDoFRigidBodyMotion`. This member, `motion_`, is responsible for integrating pressure and shear force around the given patch, and solve for the resulting  translation and rotations. The resulting motion is returned as a `septernion`, through the `transform()` function. This is tested on two examples. Trajectory of a projectile motion, and a simple pendulum. Results look alright to me. 

----------------

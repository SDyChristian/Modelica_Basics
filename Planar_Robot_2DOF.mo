model Planar_Robot_2DOF
  // Robot parameters
  parameter Real L1 = 1.0 "Length of the first link";
  parameter Real L2 = 1.0 "Length of the second link";
  parameter Real m1 = 1.0 "Mass of the first link";
  parameter Real m2 = 1.0 "Mass of the second link";
  parameter Real I1 = 0.1 "Moment of inertia of the first link";
  parameter Real I2 = 0.1 "Moment of inertia of the second link";
  parameter Real d1 = 0.1 "Friction coeficient of the first link";
  parameter Real d2 = 0.1 "Friction coeficient of the second link";
  parameter Real g = 9.81 "Acceleration due to gravity";
  
  // State variables (joint angles and velocities)
  Real theta1(start = 0.0) "Angle of the first joint (in radians)";
  Real theta2(start = 0.0) "Angle of the second joint (in radians)";
  Real dtheta1(start = 0.0) "Angular velocity of the first joint";
  Real dtheta2(start = 0.0) "Angular velocity of the second joint";
  
  // Input torques applied to the joints
  input Real tau1 = 0 "Torque applied to the first joint";
  input Real tau2 = 0 "Torque applied to the second joint";
  
  // Matrices for system dynamics
  Real M[2,2] "Mass (inertia) matrix";
  Real C[2,2] "Coriolis matrix";
  Real D[2,2] "Friction matrix";
  Real G[2] "Gravitational forces vector";
  Real tau[2] "Torques vector";
  Real ddtheta[2] "Angular accelerations";
  
  // Equations for system dynamics
  equation
    // Calculate the elements of the mass (inertia) matrix
    M[1,1] = I1 + m1*L1^2 + m2*(L1^2 + L2^2 + 2*L1*L2*cos(theta2));
    M[1,2] = m2*(L2^2 + L1*L2*cos(theta2));
    M[2,1] = M[1,2];
    M[2,2] = I2 + m2*L2^2;
    
    // Calculate the Coriolis matrix (terms involving velocity)
    C[1,1] = -m2*L1*L2*sin(theta2)*dtheta2;
    C[1,2] = -m2*L1*L2*sin(theta2)*(dtheta1 + dtheta2);
    C[2,1] = m2*L1*L2*sin(theta2)*dtheta1;
    C[2,2] = 0;
    
    // Set Friction Matrix (Constant Matrix)
    D[1,1] = d1;
    D[1,2] = 0;
    D[2,1] = 0;
    D[2,2] = d2;
    
    // Calculate the gravitational forces vector
    G[1] = (m1*L1 + m2*L1)*g*cos(theta1) + m2*L2*g*cos(theta1 + theta2);
    G[2] = m2*L2*g*cos(theta1 + theta2);
    
    // Calculate the torques vector
    tau = {tau1, tau2};
   
    // System of equations: M * ddtheta + C * dtheta + D * dtheta + G = tau
    // Solve for angular accelerations ddtheta = M^-1 * (tau - C * dtheta - G)
    
    // Joint velocities update
    der(dtheta1) = ddtheta[1];
    der(dtheta2) = ddtheta[2];

    // Joint angles update
    der(theta1) = dtheta1;
    der(theta2) = dtheta2;

    // Solve the linear system M * ddtheta = tau - C * dtheta - D * dtheta - G
    M*ddtheta = tau - C * {dtheta1, dtheta2} - D * {dtheta1, dtheta2} - G;
    
    
end Planar_Robot_2DOF;

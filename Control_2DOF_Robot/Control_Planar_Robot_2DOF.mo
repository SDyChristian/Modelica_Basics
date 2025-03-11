/* 
Description: This model represents the dynamics 2DOF Planar Robot. 
- The parameters used were based on a real small industrial robots.
- The code was hand-written, not auto-generated
Author: Christian Trejo
Date: 02/March/2025 
Updated: 10/March/2025 
Simulation Parameters:
- Step size: 0.001 seconds
- Solver: Runge-Kutta
*/
model Control_Planar_Robot_2DOF
  /* Robot parameters */
  parameter Real L1 = 1.0  "[M] Length of the first link";
  parameter Real L2 = 0.8  "[M] Length of the second link";
  parameter Real m1 = 5.0  "[Kg] Mass of the first link";
  parameter Real m2 = 4.0  "[Kg] Mass of the second link";
  parameter Real I1 = 0.15 "[Kg.m^2] Moment of inertia of the first link";
  parameter Real I2 = 0.12 "[Kg.m^2] Moment of inertia of the second link";
  parameter Real d1 = 0.05 "[N.Rad/sec] Friction coefficient of the first link";
  parameter Real d2 = 0.05 "[N.Rad/sec]Friction coefficient of the second link";
  parameter Real g = 9.81  "[m/sec^s]Acceleration due to gravity";
  parameter Real kp1 = 2  "Proportional gain for th1";
  parameter Real kp2 = 2  "Proportional gain for th2";
  parameter Real kd1 = 1  "Proportional gain for dth1";
  parameter Real kd2 = 1  "Proportional gain for dth2";
  parameter Real Kp[2,2] = [kp1, 0; 0, kp2] "Proportional gain matrix";
  parameter Real Kd[2,2] = [kd1, 0; 0, kd2]"Derivative gain matrix";
  
  
  /* State variables (joint angles and velocities) */
  Real th1(start = 0.0)  "[Rad] Angle of the first joint";
  Real th2(start = 0.0)  "[Rad] Angle of the second joint";
  Real dth1(start = 0.0) "[Rad/sec] Angular velocity of the first joint";
  Real dth2(start = 0.0) "[Rad/sec] Angular velocity of the second joint";
  Real th_d[2]           "Desired angular position";
  Real dth_d[2]          "Desired angular velocity";
  Real ddth_d[2]         "Desired angular acceleration";
  Real e[2]              "Error in position";
  Real de[2]             "Error in velocity";
  
  /* Input torques applied to the joints */
  input Real tau1 = 0 "Torque applied to the first joint";
  input Real tau2 = 0 "Torque applied to the second joint";
  
  /* Matrices for system dynamics */
  Real M[2,2]  "Mass (inertia) matrix";
  Real C[2,2]  "Coriolis matrix";
  Real D[2,2]  "Friction matrix";
  Real G[2]    "Gravitational forces vector";
  Real ddth[2] "Angular accelerations";
  Real tau[2]  "Torques vector";
  Real h[2]    "Non-linear compensation terms (Control purpose)";
  Real u_k[2]  "Second Order Kinematic Control";
  
  
  /* Equations for system dynamics */
  equation
    /* Calculate the elements of the mass (inertia) matrix */
    M[1,1] = I1 + m1*L1^2 + m2*(L1^2 + L2^2 + 2*L1*L2*cos(th2));
    M[1,2] = m2*(L2^2 + L1*L2*cos(th2));
    M[2,1] = M[1,2];
    M[2,2] = I2 + m2*L2^2;
    
    /* Calculate the Coriolis matrix (terms involving velocity) */
    C[1,1] = -m2*L1*L2*sin(th2)*dth2;
    C[1,2] = -m2*L1*L2*sin(th2)*(dth1 + dth2);
    C[2,1] = m2*L1*L2*sin(th2)*dth1;
    C[2,2] = 0;
    
    /* Set Friction Matrix (Constant Matrix) */
    D[1,1] = d1;
    D[1,2] = 0;
    D[2,1] = 0;
    D[2,2] = d2;
    
    /* Calculate the gravitational forces vector */
    G[1] = (m1*L1 + m2*L1)*g*cos(th1) + m2*L2*g*cos(th1 + th2);
    G[2] = m2*L2*g*cos(th1 + th2);
    
    /* ---- Controller ---- */
    
    /* Desired trayectory */
    th_d[1] = 2.0*sin(time);
    th_d[2] = 1.5*cos(time);
    dth_d = der(th_d);
    ddth_d = der(dth_d);
    
    /* Errors */
    e = {th1,th2}-th_d;
    de = der(e);
    
    u_k = ddth_d - Kd * de - Kp * e; 
    
    h = C * {dth1, dth2} + D * {dth1, dth2} + G;
    
    tau = M*u_k + h;
    
    /* Joint velocities update */
    der(dth1) = ddth[1];
    der(dth2) = ddth[2];

    /* Joint angles update */
    der(th1) = dth1;
    der(th2) = dth2;

    /* Solve ODE equation: M * ddth = tau - C * dth - D * dth - G */
    M*ddth = tau - C * {dth1, dth2} - D * {dth1, dth2} - G;
    
end Control_Planar_Robot_2DOF;

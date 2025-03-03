model DC_Motor_Model
  // Model Components
  Modelica.Electrical.Analog.Basic.Ground GND annotation(
    Placement(transformation(origin = {-44, -30}, extent = {{-10, -10}, {10, 10}})));
  Modelica.Electrical.Analog.Basic.Resistor R_a(R = 2) annotation(
    Placement(transformation(origin = {-24, 40}, extent = {{-10, -10}, {10, 10}})));
  Modelica.Electrical.Analog.Basic.Inductor L_a(L = 0.23) annotation(
    Placement(transformation(origin = {16, 40}, extent = {{-10, -10}, {10, 10}})));
  Modelica.Electrical.Analog.Basic.RotationalEMF EMF(k = 0.235) annotation(
    Placement(transformation(origin = {36, 12}, extent = {{-10, -10}, {10, 10}})));
  Modelica.Electrical.Analog.Sources.SignalVoltage Voltage annotation(
    Placement(transformation(origin = {-44, 6}, extent = {{10, -10}, {-10, 10}}, rotation = 90)));
  Modelica.Blocks.Sources.Step Signal(height = 12) annotation(
    Placement(transformation(origin = {-80, 6}, extent = {{-10, -10}, {10, 10}})));
  Modelica.Mechanics.Rotational.Components.Inertia Shaft_Inertia(J = 0.000052) annotation(
    Placement(transformation(origin = {74, 12}, extent = {{-10, -10}, {10, 10}})));
    
  // Variables
    
equation
  // Connections bewteen components
  connect(R_a.n, L_a.p) annotation(
    Line(points = {{-14, 40}, {6, 40}}, color = {0, 0, 255}));
  connect(L_a.n, EMF.p) annotation(
    Line(points = {{26, 40}, {36, 40}, {36, 22}}, color = {0, 0, 255}));
  connect(EMF.n, GND.p) annotation(
    Line(points = {{36, 2}, {36, -20}, {-44, -20}}, color = {0, 0, 255}));
  connect(Voltage.p, R_a.p) annotation(
    Line(points = {{-44, 16}, {-44, 40}, {-34, 40}}, color = {0, 0, 255}));
  connect(Voltage.n, GND.p) annotation(
    Line(points = {{-44, -4}, {-44, -20}}, color = {0, 0, 255}));
  connect(Voltage.v, Signal.y) annotation(
    Line(points = {{-56, 6}, {-69, 6}}, color = {0, 0, 127}));
  connect(Shaft_Inertia.flange_a, EMF.flange) annotation(
    Line(points = {{64, 12}, {46, 12}}));
  annotation(
    uses(Modelica(version = "4.0.0")));
end DC_Motor_Model;

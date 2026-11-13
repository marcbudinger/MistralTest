within ;
package SubmarineSTI
  model FullModel
    Modelica.Electrical.Analog.Basic.Ground ground
      annotation (Placement(transformation(extent={{-30,-8},{-10,12}})));
    Modelica.Electrical.Analog.Basic.Resistor resistor(R=8e-3)
      annotation (Placement(transformation(extent={{-74,44},{-54,64}})));
    Modelica.Electrical.Analog.Basic.Inductor inductor(L=20e-3)
      annotation (Placement(transformation(extent={{-46,44},{-26,64}})));
    Modelica.Electrical.Analog.Sources.TableVoltage tableVoltage(table=[0,0; 50,
          700; 100,700; 100,0; 150,0]) annotation (Placement(transformation(
          extent={{-10,-10},{10,10}},
          rotation=270,
          origin={-80,32})));
    Modelica.Electrical.Analog.Basic.Ground ground1
      annotation (Placement(transformation(extent={{-90,-8},{-70,12}})));
    Modelica.Mechanics.Rotational.Components.Inertia inertia_motor(J=3e3)
      annotation (Placement(transformation(extent={{2,20},{22,40}})));
    Modelica.Mechanics.Rotational.Components.Damper friction_motor(d=10)
      annotation (Placement(transformation(
          extent={{-10,-10},{10,10}},
          rotation=90,
          origin={0,12})));
    Modelica.Mechanics.Rotational.Components.Fixed fixed
      annotation (Placement(transformation(extent={{-10,-12},{10,8}})));
    Modelica.Mechanics.Rotational.Components.SpringDamper shaft_stiffness_loss(
        c = 6.5e6, d=25000, phi_rel(displayUnit = "rad"))
      annotation (Placement(transformation(extent={{32,20},{52,40}})));
    Modelica.Mechanics.Rotational.Components.Inertia Propeler_Inertia(J=150e3)
      annotation (Placement(transformation(extent={{92,16},{120,44}})));
    Modelica.Mechanics.Rotational.Components.BearingFriction Propeler_friction(
        tau_pos=[0,9000])
      annotation (Placement(transformation(extent={{62,20},{82,40}})));
    Modelica.Mechanics.Rotational.Sources.QuadraticSpeedDependentTorque
      Propeler_quadratic_torque(
      w_nominal=120*2*3.14/60,
      TorqueDirection=false,
      tau_nominal=-1.75e5)
      annotation (Placement(transformation(extent={{156,19},{134,41}})));
  Modelica.Electrical.Analog.Basic.EMF           emf(k = 55.7)  annotation (
      Placement(visible = true, transformation(origin = {-20, 32}, extent = {{-10, -10}, {10, 10}}, rotation = 0)));
  equation
    connect(Propeler_Inertia.flange_b, Propeler_quadratic_torque.flange) annotation (
      Line(points = {{120, 30}, {134, 30}, {134, 30}, {134, 30}}));
    connect(inductor.p, resistor.n) annotation (
      Line(points = {{-46, 54}, {-54, 54}}, color = {0, 0, 255}, smooth = Smooth.None));
    connect(tableVoltage.p, resistor.p) annotation (
      Line(points = {{-80, 42}, {-80, 54}, {-74, 54}}, color = {0, 0, 255}, smooth = Smooth.None));
    connect(ground1.p, tableVoltage.n) annotation (
      Line(points = {{-80, 12}, {-80, 22}}, color = {0, 0, 255}, smooth = Smooth.None));
    connect(friction_motor.flange_b, inertia_motor.flange_a) annotation (
      Line(points = {{0, 22}, {0, 30}, {2, 30}}, color = {0, 0, 0}, smooth = Smooth.None));
    connect(fixed.flange, friction_motor.flange_a) annotation (
      Line(points = {{0, -2}, {0, 2}}, color = {0, 0, 0}, smooth = Smooth.None));
    connect(inertia_motor.flange_b, shaft_stiffness_loss.flange_a) annotation (
      Line(points = {{22, 30}, {32, 30}}, color = {0, 0, 0}, smooth = Smooth.None));
    connect(shaft_stiffness_loss.flange_b, Propeler_friction.flange_a) annotation (
      Line(points = {{52, 30}, {62, 30}}, color = {0, 0, 0}, smooth = Smooth.None));
    connect(Propeler_friction.flange_b, Propeler_Inertia.flange_a) annotation (
      Line(points = {{82, 30}, {92, 30}}, color = {0, 0, 0}, smooth = Smooth.None));
  connect(inductor.n, emf.p) annotation (
      Line(points = {{-26, 54}, {-20, 54}, {-20, 42}}, color = {0, 0, 255}));
  connect(ground.p, emf.n) annotation (
      Line(points = {{-20, 12}, {-20, 22}}, color = {0, 0, 255}));
  connect(emf.flange, inertia_motor.flange_a) annotation (
      Line(points = {{-10, 32}, {2, 32}, {2, 30}}));
    annotation (Diagram(coordinateSystem(preserveAspectRatio=false, extent={{
              -100,-20},{180,100}}), graphics), Icon(coordinateSystem(extent={{
              -100,-20},{180,100}})));
  end FullModel;

  model ReducedModel1
    Modelica.Electrical.Analog.Basic.Ground ground
      annotation (Placement(transformation(extent={{-30,-8},{-10,12}})));
    Modelica.Electrical.Analog.Basic.Resistor resistor(R=8e-3)
      annotation (Placement(transformation(extent={{-74,44},{-54,64}})));
    Modelica.Electrical.Analog.Basic.Inductor inductor(L=20e-3)
      annotation (Placement(transformation(extent={{-46,44},{-26,64}})));
    Modelica.Electrical.Analog.Sources.TableVoltage tableVoltage(table=[0,0; 50,
          700; 100,700; 100,0; 150,0]) annotation (Placement(transformation(
          extent={{-10,-10},{10,10}},
          rotation=270,
          origin={-80,32})));
    Modelica.Electrical.Analog.Basic.Ground ground1
      annotation (Placement(transformation(extent={{-90,-8},{-70,12}})));
    Modelica.Mechanics.Rotational.Components.Damper friction_motor(d=10)
      annotation (Placement(transformation(
          extent={{-10,-10},{10,10}},
          rotation=90,
          origin={0,12})));
    Modelica.Mechanics.Rotational.Components.Fixed fixed
      annotation (Placement(transformation(extent={{-10,-12},{10,8}})));
    Modelica.Mechanics.Rotational.Components.Inertia Propeler_Inertia(J=150e3)
      annotation (Placement(visible = true, transformation(extent = {{14, 16}, {42, 44}}, rotation = 0)));
    Modelica.Mechanics.Rotational.Sources.QuadraticSpeedDependentTorque
      Propeler_quadratic_torque(
      w_nominal=120*2*3.14/60,
      tau_nominal=-1.75e5,
      TorqueDirection=true)
      annotation (Placement(transformation(extent={{82,19},{60,41}})));
  Modelica.Thermal.HeatTransfer.Interfaces.HeatPort_a heatPort(Q_flow = -LossPower, T = T_heatPort) annotation (
      Placement(visible = true, transformation(extent = {{-10, -110}, {10, -90}}, rotation = 0), iconTransformation(extent = {{-10, -110}, {10, -90}}, rotation = 0)));
  Modelica.Electrical.Analog.Basic.RotationalEMF emf(k = 55.7) annotation (
      Placement(visible = true, transformation(origin = {-20, 32}, extent = {{-10, -10}, {10, 10}}, rotation = 0)));
  equation
    connect(inductor.p, resistor.n) annotation (
      Line(points = {{-46, 54}, {-54, 54}}, color = {0, 0, 255}, smooth = Smooth.None));
    connect(tableVoltage.p, resistor.p) annotation (
      Line(points = {{-80, 42}, {-80, 54}, {-74, 54}}, color = {0, 0, 255}, smooth = Smooth.None));
    connect(ground1.p, tableVoltage.n) annotation (
      Line(points = {{-80, 12}, {-80, 22}}, color = {0, 0, 255}, smooth = Smooth.None));
    connect(fixed.flange, friction_motor.flange_a) annotation (
      Line(points = {{0, -2}, {0, 2}}, color = {0, 0, 0}, smooth = Smooth.None));
  connect(Propeler_Inertia.flange_b, Propeler_quadratic_torque.flange) annotation (
      Line(points = {{42, 30}, {60, 30}}));
  connect(friction_motor.flange_b, Propeler_Inertia.flange_a) annotation (
      Line(points = {{0, 22}, {0, 26}, {14, 26}, {14, 30}}));
  connect(ground.p, emf.n) annotation (
      Line(points = {{-20, 12}, {-20, 22}}, color = {0, 0, 255}));
  connect(emf.p, inductor.n) annotation (
      Line(points = {{-20, 42}, {-26, 42}, {-26, 54}}, color = {0, 0, 255}));
  connect(emf.flange, Propeler_Inertia.flange_a) annotation (
      Line(points = {{-10, 32}, {14, 32}, {14, 30}}));
    annotation (Diagram(coordinateSystem(preserveAspectRatio=false, extent={{-100,
              -20},{180,100}}),      graphics), Icon(coordinateSystem(extent={{
              -100,-20},{180,100}})));
  end ReducedModel1;


  model ReducedModel2
    Modelica.Electrical.Analog.Basic.EMF emf(k=700/(120*2*3.14/60))
      annotation (Placement(transformation(extent={{-30,20},{-10,40}})));
    Modelica.Electrical.Analog.Basic.Ground ground
      annotation (Placement(transformation(extent={{-30,-8},{-10,12}})));
    Modelica.Electrical.Analog.Basic.Resistor resistor(R=8e-3)
      annotation (Placement(transformation(extent={{-74,44},{-54,64}})));
    Modelica.Electrical.Analog.Basic.Inductor inductor(L=20e-3)
      annotation (Placement(transformation(extent={{-46,44},{-26,64}})));
    Modelica.Electrical.Analog.Sources.TableVoltage tableVoltage(table=[0,0; 50,
          700; 100,700; 100,0; 150,0]) annotation (Placement(transformation(
          extent={{-10,-10},{10,10}},
          rotation=270,
          origin={-80,32})));
    Modelica.Electrical.Analog.Basic.Ground ground1
      annotation (Placement(transformation(extent={{-90,-8},{-70,12}})));
    Modelica.Mechanics.Rotational.Components.Inertia Propeler_Inertia(J=150e3)
      annotation (Placement(transformation(extent={{16,16},{44,44}})));
    Modelica.Blocks.Continuous.TransferFunction transferFunction(a={20e-3*150e3
          /55.73^2,150e3/55.73^2*8e-3,1}, b={0,150e3/55.73,0})
      annotation (Placement(transformation(extent={{112,20},{132,40}})));
    Modelica.Blocks.Sources.TimeTable timeTable(table=[0,0; 50,700; 100,700;
          100,0; 150,0])
      annotation (Placement(transformation(extent={{80,20},{100,40}})));
  equation
    connect(emf.p, inductor.n) annotation (Line(
        points={{-20,40},{-20,54},{-26,54}},
        color={0,0,255},
        smooth=Smooth.None));
    connect(emf.n, ground.p) annotation (Line(
        points={{-20,20},{-20,12}},
        color={0,0,255},
        smooth=Smooth.None));
    connect(inductor.p, resistor.n) annotation (Line(
        points={{-46,54},{-54,54}},
        color={0,0,255},
        smooth=Smooth.None));
    connect(tableVoltage.p, resistor.p) annotation (Line(
        points={{-80,42},{-80,54},{-74,54}},
        color={0,0,255},
        smooth=Smooth.None));
    connect(ground1.p, tableVoltage.n) annotation (Line(
        points={{-80,12},{-80,22}},
        color={0,0,255},
        smooth=Smooth.None));
    connect(emf.flange, Propeler_Inertia.flange_a) annotation (Line(
        points={{-10,30},{16,30}},
        color={0,0,0},
        smooth=Smooth.None));
    connect(timeTable.y, transferFunction.u) annotation (Line(
        points={{101,30},{110,30}},
        color={0,0,127},
        smooth=Smooth.None));
    annotation (Diagram(coordinateSystem(preserveAspectRatio=false, extent={{-100,
              -20},{180,100}}),      graphics), Icon(coordinateSystem(extent={{
              -100,-20},{180,100}})));
  end ReducedModel2;
  annotation (uses(Modelica(version="3.2.3")));
end SubmarineSTI;

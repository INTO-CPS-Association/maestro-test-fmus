within ;
model massive

  Modelica.Blocks.Interfaces.RealInput signal_1_in
    annotation (Placement(transformation(extent={{-98,46},{-58,86}})));
  Modelica.Blocks.Interfaces.RealInput signal_2_in
    annotation (Placement(transformation(extent={{-96,6},{-56,46}})));
  Modelica.Blocks.Interfaces.RealInput signal_ext_in
    annotation (Placement(transformation(extent={{-98,-28},{-58,12}})));
  Modelica.Blocks.Interfaces.RealOutput signal_1_out
    annotation (Placement(transformation(extent={{30,48},{50,68}})));
  Modelica.Blocks.Interfaces.RealOutput signal_2_out
    annotation (Placement(transformation(extent={{30,16},{50,36}})));
  Modelica.Blocks.Interfaces.RealOutput signal_ext_out
    annotation (Placement(transformation(extent={{30,-14},{50,6}})));
  Modelica.Blocks.Continuous.FirstOrder in1(
    T=1e-8,
    initType=Modelica.Blocks.Types.Init.InitialState,
    y_start=sig1Init)
    annotation (Placement(transformation(extent={{-50,58},{-30,78}})));
  Modelica.Blocks.Continuous.FirstOrder in2(
    T=1e-8,
    initType=Modelica.Blocks.Types.Init.InitialState,
    y_start=sig2Init)
    annotation (Placement(transformation(extent={{-48,16},{-28,36}})));
  Modelica.Blocks.Continuous.FirstOrder ext(
    T=1e-8,
    initType=Modelica.Blocks.Types.Init.InitialState,
    y_start=extInit)
    annotation (Placement(transformation(extent={{-48,-18},{-28,2}})));

    parameter Real sig1Init = 0;
    parameter Real sig2Init = 0;
    parameter Real extInit = 0;
equation
  connect(signal_1_in, in1.u) annotation (Line(
      points={{-78,66},{-66,66},{-66,68},{-52,68}},
      color={0,0,127},
      smooth=Smooth.None));
  connect(in1.y, signal_ext_out) annotation (Line(
      points={{-29,68},{-10,68},{-10,-4},{40,-4}},
      color={0,0,127},
      smooth=Smooth.None));
  connect(signal_2_in, in2.u) annotation (Line(
      points={{-76,26},{-50,26}},
      color={0,0,127},
      smooth=Smooth.None));
  connect(in2.y, signal_2_out) annotation (Line(
      points={{-27,26},{40,26}},
      color={0,0,127},
      smooth=Smooth.None));
  connect(signal_ext_in, ext.u) annotation (Line(
      points={{-78,-8},{-50,-8}},
      color={0,0,127},
      smooth=Smooth.None));
  connect(ext.y, signal_1_out) annotation (Line(
      points={{-27,-8},{2,-8},{2,58},{40,58}},
      color={0,0,127},
      smooth=Smooth.None));
  annotation (uses(Modelica(version="3.2.1")), Diagram(coordinateSystem(
          preserveAspectRatio=false, extent={{-100,-100},{100,100}}), graphics));
end massive;

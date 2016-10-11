within ;
model sink
  Modelica.Blocks.Interfaces.RealInput sig1
    annotation (Placement(transformation(extent={{-72,30},{-32,70}})));
  Modelica.Blocks.Interfaces.RealInput sig2
    annotation (Placement(transformation(extent={{-72,-12},{-32,28}})));
  Modelica.Blocks.Interfaces.RealInput sig3
    annotation (Placement(transformation(extent={{-74,-50},{-34,-10}})));
  Modelica.Blocks.Interaction.Show.RealValue realValue
    annotation (Placement(transformation(extent={{-4,40},{16,60}})));
  Modelica.Blocks.Interaction.Show.RealValue realValue1
    annotation (Placement(transformation(extent={{-4,-2},{16,18}})));
  Modelica.Blocks.Interaction.Show.RealValue realValue2
    annotation (Placement(transformation(extent={{-4,-40},{16,-20}})));
equation
  connect(sig1, realValue.numberPort)
    annotation (Line(points={{-52,50},{-5.5,50},{-5.5,50}}, color={0,0,127}));
  connect(sig2, realValue1.numberPort)
    annotation (Line(points={{-52,8},{-5.5,8},{-5.5,8}}, color={0,0,127}));
  connect(sig3, realValue2.numberPort)
    annotation (Line(points={{-54,-30},{-5.5,-30}}, color={0,0,127}));
  annotation (uses(Modelica(version="3.2.1")));
end sink;

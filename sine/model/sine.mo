within ;
model sine

  Modelica.Blocks.Sources.Sine sine(
    amplitude=amp,
    freqHz=sineFreq,
    phase=phase,
    offset=sineOffset,
    startTime=start)
    annotation (Placement(transformation(extent={{-56,22},{-36,42}})));
  Modelica.Blocks.Interfaces.RealOutput y
    annotation (Placement(transformation(extent={{6,22},{26,42}})));
  parameter Real pi = Modelica.Constants.pi;
  parameter Real amp = 1 "Sinus Amplite";
  parameter Real sineFreq = 0.4 "Sinus frequency";
  parameter Real phase = 0 "Sinus phase";
  parameter Real sineOffset = 0 "Sinus offset";
  parameter Real start = 0 "Sinus startTime";
equation
  connect(sine.y, y)
    annotation (Line(points={{-35,32},{-8,32},{16,32}}, color={0,0,127}));
  annotation (uses(Modelica(version="3.2.1")));
end sine;

within ;
model pump

  Modelica.Blocks.Interfaces.RealOutput mdot_out
    annotation (Placement(transformation(extent={{-74,0},{-54,20}})));
    parameter Real mdotInit = 1;
    parameter Real tBegin = 2;
    parameter Real tLength = 1e-6;
    parameter Real h = 0;
    Real tStart;
    Real tEnd = tBegin + tLength;
    Real T;
    Real mdot;
    Real mdotFake;
  Modelica.Blocks.Sources.Clock clock
    annotation (Placement(transformation(extent={{-82,-42},{-62,-22}})));
initial equation
  mdot = mdotInit;
  tStart = tBegin;
  der(tStart)=0;
  mdotFake = mdotInit;
equation
  der(tStart)=0;
  der(mdotFake)=mdotInit*h;
  T = clock.y;
  if ((T>=tStart) and (T<=tEnd)) then
    der(mdot) = h/tLength;
  else
    der(mdot) = 0;
  end if;
  mdot = mdot_out;
  annotation (uses(Modelica(version="3.2.1")));
end pump;

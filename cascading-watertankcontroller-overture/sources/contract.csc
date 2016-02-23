
//shared design parameters
sdp real maxlevel1;
sdp real minlevel1; 
sdp real maxlevel2;
sdp real minlevel2;

-- Monitored variables (seen from the DE controller)
monitored real level1;
monitored real level2;

-- Controlled variables (seen from the DE controller)
controlled bool valveState1;
controlled bool valveState2;


package HydroPowerPlant
  model Forebay
    //parameter Modelica.SIunits.Volume absolute_volume = 1;
    //Modelica.SIunits.Volume available_volume;
    //WaterInput current_income_water;
  equation
//watch for filling forebay
  end Forebay;

  model Tailrace
    //parameter Modelica.Fluid.Sources.Boundary_pT bounds = 1;???
    //parameter Modelica.Fluid.Pipes.DynamicPipe pipe = 1;??
    //WaterInput  from_river;
    //WaterOutput for_forebay;
  equation

  end Tailrace;

  model IntakeStructure
    // parameter Modelica.SIunits.Volume basin = 1;
    // parameter Modelica.SIunits.Volume raw_collector = 1;
    // parameter Modelica.SIunits.Time delta_for_purification = 1; darsi
    // parameter Modelica.SIunits.Volume tank = 1;
    // parameter Modelica.Fluid.Machines.PrescribedPump pump = 1;
    // parameter Modelica.Fluid.Vessels.OpenTank tank_for_pureficated_water = 1;
    // WaterInput from_river;
    // WaterOutput to_surge (or to purification station);
  equation

  end IntakeStructure;

  model Penstock
    // parameter Modelica.SIunits.Area tube_cross-sectional_area;
    // WaterInput from_surge_chamber;
    // WaterOutput to_turbine;
  equation

  end Penstock;

  model SurgeChamber
    // parameter Modelica.SIunits.Volume chamber_volume;
    // WaterInput from_intake_structure;
    // WaterOutput to_penstock;
  equation

  end SurgeChamber;

  model HydraulicTurbines
    // parameter Modelica.Mechanics.Rotational.Interfaces.
    // WaterInput from_penstock;
    // WaterOutput to_draft_tube;
  equation

  end HydraulicTurbines;

  model Powerhouse
  equation

  end Powerhouse;

  model DraftTube
    // parameter Modelica.SIunits.Area dtube_cross-sectional_area;
    // WaterInput from_turbine;
    // WaterOutput to_drop_point;
  equation

  end DraftTube;

  connector DynamicInput
  end DynamicInput;

  connector DynamicOutput
  end DynamicOutput;

  connector WaterInput
    input Modelica.SIunits.Pressure p "Contact pressure";
    input Modelica.SIunits.Temperature T "Contact temperature";
    input Modelica.SIunits.MassFlowRate m_dot "Mass flow rate through the contact";
  end WaterInput;

  connector WaterOutput
    output Modelica.SIunits.Pressure p "Contact pressure";
    output Modelica.SIunits.Temperature T "Contact temperature";
    output Modelica.SIunits.MassFlowRate m_dot "Mass flow rate through the contact";
  end WaterOutput;

  connector ElectricityInput
    //input Modelica.Electrical.Analog.Sources.ConstantCurrent in_current;
    //input ... in_res;
    //input ... in_amp;
    //input ... in_pow;
  end ElectricityInput;

  connector ElectricityOutput
    //output Modelica.Electrical.Analog.Sources.ConstantCurrent out_current;
    //output ... out_res;
    //output ... out_amp;
    //output ... out_pow;
  end ElectricityOutput;

  model River
    parameter Modelica.SIunits.Height H_r = 50; //"Initial water level above intake"
    parameter Modelica.SIunits.Length L = 500; //"Length of the reservoir"
    parameter Modelica.SIunits.Length w = 100; //"Bed width of the reservoir"
    parameter Modelica.SIunits.Conversions.NonSIunits.Angle_deg alpha = 30; //"The angle of the reservoir walls (zero angle corresponds to vertical walls)"
    parameter Real f = 0.0008; //"Friction factor of the reservoir"
    parameter Modelica.SIunits.Temperature T_i = Const.T_i; //"Initial temperature of the water"
    
  
    Modelica.SIunits.Area A; //"Vertical cros section";
    Modelica.SIunits.Mass m; //"Water mass";
    Modelica.SIunits.MassFlowRate m_dot; //"Water mass flow rate";
    Modelica.SIunits.VolumeFlowRate V_o_dot "Outlet flow rate", V_i_dot "Inlet flow rate", V_dot "Vertical flow rate";
    Modelica.SIunits.Velocity v; //"Water velocity";
    Modelica.SIunits.Momentum M; //"Water momentum";
    Modelica.SIunits.Force F_f; //"Friction force";
    Modelica.SIunits.Height H; //"Water height";
    Modelica.SIunits.Pressure p_2; //"Outside pressure";
    
    OpenHPL.Interfaces.Contact n(p=p_2); //"Outflow from reservoir"
  equation
    // Define vertiacal cross section of the reservoir
    A = H * (w + 2 * H * Modelica.Math.tan(Modelica.SIunits.Conversions.from_deg(alpha)));
    // Define water mass
    m = Const.rho * A * L;
    // Define volumetric water flow rate
    V_dot = V_i_dot - V_o_dot;
    // Define mass water flow rate
    m_dot = Const.rho * V_dot;
    // Define water velocity
    v = m_dot / Const.rho / A;
    // Define momentrumn
    M = L * m_dot;
    // Define friction term
    F_f = 1 / 8 * Const.rho * f * L * (w + 2 * H / Modelica.Math.cos(alpha)) * v * abs(v);
     
    // define derivatives of momentum and mass
    der(M) = A * (Const.p_a - p_2) + Const.g * Const.rho * A * H - F_f + Const.rho / A * (V_i_dot ^ 2 - V_o_dot ^ 2);
    der(m) = m_dot;
    // define output pressure
    p_2 = Const.p_a + Const.g * Const.rho * H;
    // output flow conector
    n.m_dot = -Const.rho * V_o_dot;
    // output temperature conector
    n.T = T_i;
  end River;

record ListOfConstants
parameter Modelica.SIunits.Acceleration g = Modelica.Constants.g_n "gravity";
  parameter Modelica.SIunits.Density rho = 997.0 "density";
  parameter Modelica.SIunits.DynamicViscosity mu = 0.89e-3 "dynamic viscosity of water";
  parameter Modelica.SIunits.Height eps = 5e-2 "pipe roughness height";
  parameter Modelica.SIunits.Pressure p_a = 1.013e5 "Atmospheric pressure";
  parameter Modelica.SIunits.Compressibility beta = 4.5e-10 "water compressibility";
  parameter Modelica.SIunits.Compressibility beta_total = 1 / rho / 1000 ^ 2 "total compressibility";
  parameter Modelica.SIunits.VolumeFlowRate V_0 = 21.000 "Initial flow rate through the system";
  

end ListOfConstants;
  annotation(
    uses(Modelica(version = "3.2.3")));
end HydroPowerPlant;

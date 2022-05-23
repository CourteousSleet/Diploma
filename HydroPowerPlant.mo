package HydroPowerPlant
  extends Icons.MainPicture;
  //within HydroPowerPlant;

  package Icons "Пакет с иконками для элементов системы"
    extends Modelica.Icons.IconsPackage;

    package MainPicture "Логотип проекта"
      annotation(
        Icon(graphics = {Ellipse(origin = {-62, 61}, fillPattern = FillPattern.Solid, extent = {{-8, 7}, {8, -7}}), Line(origin = {0, 54}, points = {{-62, 0}, {62, 0}, {62, 0}}), Ellipse(origin = {63, 61}, fillPattern = FillPattern.Solid, extent = {{-7, 7}, {7, -7}}), Line(origin = {1, 68}, points = {{-63, 0}, {63, 0}, {63, 0}}), Line(origin = {-70, 18}, points = {{0, 44}, {0, -44}, {0, -44}}), Line(origin = {70, 17}, points = {{0, 43}, {0, -43}, {0, -43}}), Ellipse(origin = {-62, -23}, fillPattern = FillPattern.Solid, extent = {{-8, -7}, {8, 7}}), Ellipse(origin = {63, -22}, fillPattern = FillPattern.Solid, extent = {{-7, 8}, {7, -8}}), Line(origin = {-54, 15}, points = {{0, -39}, {0, 39}, {0, 39}}), Line(origin = {56, 16}, points = {{0, -38}, {0, 38}, {0, 38}}), Ellipse(origin = {-12, 24}, fillPattern = FillPattern.Solid, extent = {{-4, 4}, {4, -4}}), Ellipse(origin = {12, 24}, fillPattern = FillPattern.Solid, extent = {{-4, 4}, {4, -4}}), Ellipse(origin = {-46, 24}, fillPattern = FillPattern.Solid, extent = {{-4, 4}, {4, -4}}), Ellipse(origin = {-26, 24}, fillPattern = FillPattern.Solid, extent = {{-4, 4}, {4, -4}}), Ellipse(origin = {26, 24}, fillPattern = FillPattern.Solid, extent = {{-4, 4}, {4, -4}}), Ellipse(origin = {46, 24}, fillPattern = FillPattern.Solid, extent = {{-4, 4}, {4, -4}}), Line(origin = {-36, 28}, points = {{-10, 0}, {10, 0}, {10, 0}}), Line(origin = {0, 28}, points = {{-12, 0}, {12, 0}, {12, 0}}), Line(origin = {36, 28}, points = {{-10, 0}, {10, 0}, {10, 0}}), Line(origin = {-36, 20}, points = {{-10, 0}, {10, 0}, {10, 0}}), Line(origin = {0, 20}, points = {{-12, 0}, {12, 0}, {12, 0}}), Line(origin = {36, 20}, points = {{-10, 0}, {10, 0}, {10, 0}}), Line(origin = {-50, 3}, points = {{0, 21}, {0, -21}, {0, -21}}), Line(origin = {-22, 3}, points = {{0, 21}, {0, -21}, {0, -21}}), Line(origin = {-16, 3}, points = {{0, 21}, {0, -21}, {0, -21}}), Line(origin = {16, 3}, points = {{0, 21}, {0, -21}, {0, -21}}), Line(origin = {-42, 3}, points = {{0, 21}, {0, -21}, {0, -21}}), Line(origin = {-30, 3}, points = {{0, 21}, {0, -21}, {0, -21}}), Line(origin = {-8, 3}, points = {{0, 21}, {0, -21}, {0, -21}}), Line(origin = {8, 3}, points = {{0, 21}, {0, -21}, {0, -21}}), Line(origin = {23.2071, 3.14645}, points = {{-1.20711, 20.8536}, {-1.20711, -21.1464}, {0.792893, -19.1464}}), Line(origin = {30, 3}, points = {{0, 21}, {0, -21}, {0, -21}}), Line(origin = {42, 3}, points = {{0, 21}, {0, -21}, {0, -21}}), Line(origin = {50, 3}, points = {{0, 21}, {0, -21}, {0, -21}}), Ellipse(origin = {-46, -18}, fillPattern = FillPattern.Solid, extent = {{-4, 4}, {4, -4}}), Ellipse(origin = {-26, -18}, fillPattern = FillPattern.Solid, extent = {{-4, 4}, {4, -4}}), Ellipse(origin = {-12, -18}, fillPattern = FillPattern.Solid, extent = {{-4, 4}, {4, -4}}), Ellipse(origin = {12, -18}, fillPattern = FillPattern.Solid, extent = {{-4, 4}, {4, -4}}), Ellipse(origin = {26, -18}, fillPattern = FillPattern.Solid, extent = {{-4, 4}, {4, -4}}), Ellipse(origin = {46, -18}, fillPattern = FillPattern.Solid, extent = {{-4, 4}, {4, -4}}), Line(origin = {0.44, -40}, points = {{-64, -4}, {-55, -20}, {-44, 0}, {-33, -20}, {-20, 0}, {-11, -20}, {1.5, 2}, {12.75, -20}, {24.375, 0}, {35.1875, -18}, {45.5938, -2}, {56.7969, -18}, {66, -6}}, thickness = 3, smooth = Smooth.Bezier), Line(origin = {-1.52875, -58.5938}, points = {{-64, -4}, {-55, -20}, {-44, 0}, {-33, -20}, {-20, 0}, {-11, -20}, {1.5, 2}, {12.75, -20}, {24.375, 0}, {35.1875, -18}, {45.5938, -2}, {56.7969, -18}, {66, -6}}, thickness = 3, smooth = Smooth.Bezier), Line(origin = {-2.84125, -79.3751}, points = {{-64, -4}, {-55, -20}, {-44, 0}, {-33, -20}, {-20, 0}, {-11, -20}, {1.5, 2}, {12.75, -20}, {24.375, 0}, {35.1875, -18}, {45.5938, -2}, {56.7969, -18}, {66, -6}}, thickness = 3, smooth = Smooth.Bezier), Rectangle(origin = {1, 61}, fillPattern = FillPattern.Solid, extent = {{-63, 7}, {63, -7}}), Rectangle(origin = {-62, 19}, fillPattern = FillPattern.Solid, extent = {{-8, 43}, {8, -43}}), Rectangle(origin = {63, 20}, fillPattern = FillPattern.Solid, extent = {{-7, 42}, {7, -42}}), Rectangle(origin = {-46, 3}, fillPattern = FillPattern.Solid, extent = {{-4, 21}, {4, -21}}), Rectangle(origin = {-36, 24}, fillPattern = FillPattern.Solid, extent = {{-10, 4}, {10, -4}}), Rectangle(origin = {-26, 3}, fillPattern = FillPattern.Solid, extent = {{-4, -21}, {4, 21}}), Rectangle(origin = {-12, 3}, fillPattern = FillPattern.Solid, extent = {{-4, 21}, {4, -21}}), Rectangle(origin = {12, 3}, fillPattern = FillPattern.Solid, extent = {{-4, 21}, {4, -21}}), Rectangle(origin = {-1, 24}, fillPattern = FillPattern.Solid, extent = {{-13, 4}, {13, -4}}), Rectangle(origin = {26, 3}, fillPattern = FillPattern.Solid, extent = {{-4, 21}, {4, -21}}), Rectangle(origin = {36, 24}, fillPattern = FillPattern.Solid, extent = {{-10, 4}, {10, -4}}), Rectangle(origin = {46, 3}, fillPattern = FillPattern.Solid, extent = {{4, 21}, {-4, -21}})}));
    end MainPicture;
  end Icons;

  package Connectors "Данный пакет содержит все элементы соединений, необходимые для симуляций в модели"
    extends Modelica.Icons.InterfacesPackage;

    package IWater
    
      connector WaterConnector
      
        Modelica.SIunits.Pressure p "Давление на контакте";
        flow Modelica.SIunits.MassFlowRate mdot "Массовый расход на контакте";
      
      end WaterConnector;
    
      connector WaterInput "Входной элемент соединения"
        extends WaterConnector;
      end WaterInput;
      
      connector WaterOutput "Выходной элемент соединения"
        extends WaterConnector;
      end WaterOutput;

      partial model TwoPortWaterConnector "Модель двойного элемента соединения"
      
        WaterInput i;
        
        WaterOutput o;
      
      end TwoPortWaterConnector;

      partial model MeticulousTwoPortWaterConnector
      
        Modelica.SIunits.MassFlowRate mdot "Массовый расход с условием постоянности";
        extends TwoPortWaterConnector;
      
      equation
        
        0 = i.mdot + o.mdot;
        mdot = i.mdot;
      
      end MeticulousTwoPortWaterConnector;

      partial model ExuberantTwoPortWaterConnector
      
        Modelica.SIunits.Pressure p "Давление на контакте";
        Modelica.SIunits.MassFlowRate mdot "Массовый расход без условия постоянности";
        extends TwoPortWaterConnector;
      
      equation
      
        p = i.p;
        i.p = o.p;
        mdot = i.mdot + o.mdot;
      
      end ExuberantTwoPortWaterConnector;
    
    end IWater;

    package ITurbine 
partial model TurbineConnector
  
  extends IWater.MeticulousTwoPortWaterConnector;
  
  input Modelica.Blocks.Interfaces.RealInput u_t "Сигнал разгона лопаток турбины";
  
  Modelica.Blocks.Interfaces.RealOutput P_out(unit = "W") "Выходящая мощность";

end TurbineConnector;

      partial model SmartTurbineConnector
      
        extends IWater.MeticulousTwoPortWaterConnector;
        
        parameter Boolean include_P_out = false "Условие вывода мощности на генератор"
        annotation (choices(checkBox = true), Dialog(group="Outputs"));
        
        input Modelica.Blocks.Interfaces.RealInput u_t "Сигнал разгона лопаток турбины";
        
        Modelica.Blocks.Interfaces.RealOutput P_out(unit = "W") if include_P_out "Выходящая мощность";
        
      end SmartTurbineConnector;
        
    end ITurbine;
  end Connectors;

  package Utilities
    extends Modelica.Icons.UtilitiesPackage;

    package Darcy_Weisbach
      extends Modelica.Icons.FunctionsPackage;

      function DarcyFrictionFactor "Функция для определения силы трения на основе коэффициента трения Дарси"
      
        extends Modelica.Icons.Function;
        import Modelica.Constants.pi;
        
        input Modelica.SIunits.Velocity v           "Скорость потока";
        input Modelica.SIunits.Diameter D           "Диаметр трубы";
        input Modelica.SIunits.Length L             "Длина трубы";
        input Modelica.SIunits.Density rho          "Плотность среды";
        input Modelica.SIunits.DynamicViscosity mu  "Динамическая вязкость воды";
        input Modelica.SIunits.Height p_eps         "Высота шероховатости трубы";
        
        output Modelica.SIunits.Force F_f           "Сила трения";
        
        protected
          Modelica.SIunits.ReynoldsNumber N_Re      "Число Рейнольдса";
          Real f                                    "Коэффициент трения Дарси";
      
      algorithm
        
        N_Re := rho * abs(v) * D / mu;
        f := CalculateByReynoldsNumber(N_Re, D, p_eps);
        F_f := 0.5 * pi * f * rho * L * v * abs(v) * D / 4;

      end DarcyFrictionFactor;

      function CalculateByReynoldsNumber
      
        extends Modelica.Icons.Function;
        
        input Modelica.SIunits.ReynoldsNumber N_Re "Число Рейнольдса";
        input Modelica.SIunits.Diameter D          "Диаметр трубы";
        input Modelica.SIunits.Height p_eps        "Высота шероховатости трубы";
        
        output Real fD                             "Коэффициент трения Дарси";
        
        protected
          Real arg;
        
          // Вычисляем удельную энтальпию
          Modelica.SIunits.ReynoldsNumber N_Re_lam = 2100, N_Re_tur = 2300;
          Real X[4, 4], Y[4], K[4];
      
      algorithm
        
        X := [N_Re_lam ^ 3, N_Re_lam ^ 2, N_Re_lam, 1; N_Re_tur ^ 3, N_Re_tur ^ 2, N_Re_tur, 1; 3 * N_Re_lam ^ 2, 2 * N_Re_lam, 1, 0; 3 * N_Re_tur ^ 2, 2 * N_Re_tur, 1, 0];
        
        Y := {64 / N_Re_lam, 1 / (2 * log10(p_eps / 3.7 / D + 5.74 / N_Re_tur ^ 0.9)) ^ 2, -64 / N_Re_lam ^ 2, -0.25 * 0.316 / N_Re_tur ^ 1.25};
        
        K := Modelica.Math.Matrices.inv(X) * Y;
        
        arg := p_eps / 3.7 / D + 5.74 / (N_Re + Modelica.Constants.eps) ^ 0.9;
        
        if N_Re <= 0 then
          fD := 0;
        elseif N_Re <= 2100 then
          fD := 64 / N_Re;
        elseif N_Re < 2300 then
          fD := K[1] * N_Re ^ 3 + K[2] * N_Re ^ 2 + K[3] * N_Re + K[4];
        else
          fD := 1 / (2 * log10(arg)) ^ 2;
        end if;

      end CalculateByReynoldsNumber;
    
    end Darcy_Weisbach;
  
  end Utilities;

  record Presets "Запись с необходимыми константами, чтобы не импортировать каждый раз"
  
    extends Modelica.Icons.Record;
    
    parameter Modelica.SIunits.Acceleration g = Modelica.Constants.g_n          "Ускорение свободного падения";
    parameter Modelica.SIunits.Density rho = 997.0                              "Плотность воды";
    parameter Modelica.SIunits.DynamicViscosity mu = 0.89e-3                    "Динамическая вязкость воды";
    parameter Modelica.SIunits.Height p_eps = 5e-2                              "Высота шероховатости трубы";
    parameter Modelica.SIunits.Pressure p_a = 1.013e5                           "Атмосферное давление";
    parameter Modelica.SIunits.Compressibility beta = 4.5e-10                   "Сжимаемость воды";
    parameter Modelica.SIunits.Compressibility beta_total = 1 / rho / 1000 ^ 2  "Полная сжимаемость";
    parameter Modelica.SIunits.VolumeFlowRate V_0 = 15.000                      "Начальный объёмный расход в системе";
    parameter Modelica.SIunits.Frequency f_0 = 50                               "Частота тока для TorqueGenerator";
    parameter Modelica.SIunits.MolarMass M_a = 28.97e-3                         "Молярная масса воздуха при н.у.";
    
  end Presets;

  package WaterSources
    model Reservoir "Водохранилище"
        
        outer Presets data;
        
        parameter Modelica.SIunits.Height H_0 = 50 "Начальный уровень воды над водозаборным узлом";
    
        parameter Modelica.SIunits.Length l = 500 "Длина водохранилища";
        parameter Modelica.SIunits.Length w = 100 "Ширина водохранилища";
        parameter Real f = 0.0008                 "Коэффициент трения Дарси для водохранилища";
        parameter Modelica.SIunits.Mass m_0 = data.rho * H_0 * w * l;
        
        Modelica.SIunits.Area A                   "Площадь вертикального сечения объёма воды";
        Modelica.SIunits.Mass m(start = m_0)      "Масса воды";
        Modelica.SIunits.MassFlowRate mdot        "Массовый расход воды";
        Modelica.SIunits.VolumeFlowRate Vdot_i    "Входящй объёмный расход";
        Modelica.SIunits.VolumeFlowRate Vdot_o    "Исходящий объёмный расход";
        Modelica.SIunits.VolumeFlowRate Vdot      "Объёмный расход внутри водохранилища";
        Modelica.SIunits.Velocity v               "Скорость потока";
        Modelica.SIunits.Momentum M               "Количества движения потока воды";
        Modelica.SIunits.Force F_f                "Сила трения";
        Modelica.SIunits.Height h(start = H_0)    "Переменная уровня воды";
        Modelica.SIunits.Pressure p_o             "Исходящее давление";
        
        HydroPowerPlant.Connectors.IWater.WaterConnector w_o(p = p_o) "Элемент выходящего соединения";
    
    equation
        
        A = h * w;
        m = data.rho * A * l        "Масса воды в водохранилище";
        Vdot = Vdot_i - Vdot_o      "Объёмный расход";
        mdot = data.rho * Vdot      "Массовый расход";
        v = mdot / data.rho / A     "Скорость потока воды";
        M = l * mdot                "Количество движения потока";
        F_f = 1 / 8 * data.rho * f * l * w * v * abs(v) "Сила трения для движения воды вдоль водохраналища";
        
        Vdot_i = 0;
        p_o = data.p_a + data.g * data.rho * h;
        der(m) = mdot;
        
        w_o.mdot = -data.rho * Vdot_o "Элемент соеденинения для выходящего массового расхода";
    
    end Reservoir;
  end WaterSources;

  package Watercourse

    model SurgeChamber "Модель уравнительного резервуара"
    
      outer Presets data;
      import Modelica.Constants.pi;
      
      
      parameter Modelica.SIunits.Height H = 120                     "Высота стержня резервуара";
      parameter Modelica.SIunits.Length L = 140                     "Длина стержня резервуара";
      parameter Modelica.SIunits.Diameter D = 3.4                   "Диаметр резервуара";
      parameter Modelica.SIunits.Height p_eps = data.p_eps          "Высота шероховатости стенок";
      
      parameter Modelica.SIunits.VolumeFlowRate Vdot_0 = 0          "Начальный объёмный расход";
      parameter Modelica.SIunits.Height h_0 = 69.9                  "Начальный уровень воды";
      parameter Modelica.SIunits.Pressure p_ac = 4 * data.p_a       "Давление воздушной подушки внутри резервуара";
      parameter Modelica.SIunits.Temperature T_ac(displayUnit = "degC") = 298.15 "Температура воздушной подушки";
      
      Modelica.SIunits.Mass m                                       "Масса воды";
      Modelica.SIunits.Mass m_a = p_ac * A * (L - h_0 / cos_theta) * data.M_a / (Modelica.Constants.R * T_ac) "Масса воздуха внутри резервуара";
      Modelica.SIunits.Momentum M                                   "Количество движения потока";
      Modelica.SIunits.Force Mdot                                   "Разница количества движения притока и оттока";
      Modelica.SIunits.Force F                                      "Сумма сил внутри резервуара";
      Modelica.SIunits.Area A = pi * D ^ 2 / 4                      "Площадь горизонтального сечения";
      Real cos_theta = H / L                                        "Коэффициент наклона";
      Modelica.SIunits.Length l = h / cos_theta                     "Длина объёма воды";
      Modelica.SIunits.Velocity v                                   "Скорость потока";
      Modelica.SIunits.Force F_p                                    "Сила давления";
      Modelica.SIunits.Force F_f                                    "Сила трения";
      Modelica.SIunits.Force F_g                                    "Сила гравитации";
      Modelica.SIunits.Pressure p_t                                 "Давление в верхней точке резервуара";
      Modelica.SIunits.Pressure p_b                                 "Давление в нижней точке резервуара";
      
      Modelica.SIunits.Height h(start = h_0)                        "Переменная высота водного столба";
      Modelica.SIunits.VolumeFlowRate Vdot(start = Vdot_0)          "Объёмный расход";
    
      extends HydroPowerPlant.Connectors.IWater.ExuberantTwoPortWaterConnector;
    
    initial equation
      
      der(M) = 0;
      der(m) = 0;
    
    equation
    
      der(m) = mdot     "Сохранение массы в системе";
      der(M) = Mdot + F "Сохранение количества движения в системе";
    
      v = Vdot / A;
      m = data.rho * A * l;
      M = m * v;
      p_t = data.p_a;
      F_f = Utilities.Darcy_Weisbach.DarcyFrictionFactor(v, D, l, data.rho, data.mu, p_eps);
      F_p = (p_b - p_t) * A;
      
      mdot = data.rho * Vdot;
      Mdot = mdot * v;
      F = F_p - F_f - F_g;
      p_b = p           "Соединяем выходной элемент соединения";
      F_g = m * data.g * cos_theta;
    
    end SurgeChamber;

    model Pipeline
      
      outer Presets data;
      import Modelica.Constants.pi;
    
      parameter Modelica.SIunits.Length H = 25                    "Разница высоты между входным и выходным отверстиями";
      parameter Modelica.SIunits.Length L = 66                    "Длина трубы";
      parameter Modelica.SIunits.Diameter D_i = 5.8               "Диаметр входного отверстия";
      parameter Modelica.SIunits.Diameter D_o = D_i               "Диаметр выходного отверстия";
      parameter Modelica.SIunits.Height p_eps = data.p_eps        "Высота шероховатости трубы";
      
      parameter Modelica.SIunits.VolumeFlowRate Vdot_0 = data.V_0 "Начальный объёмный расход в трубе";
      
      Modelica.SIunits.Diameter D_ = 0.5 * (D_i + D_o)            "Средний диаметр трубы";
      Modelica.SIunits.Mass m                                     "Масса воды";
      Modelica.SIunits.Area A_i = D_i ^ 2 * pi / 4                "Площадь входного отверстия";
      Modelica.SIunits.Area A_o = D_o ^ 2 * pi / 4                "Площадь выходного отверстия";
      Modelica.SIunits.Area A_ = D_ ^ 2 * pi / 4                  "Средня площадь сечения трубы";
      Real cos_theta = H / L                                      "Коэффициент наклона";
      Modelica.SIunits.Velocity v                                 "Скорость воды";
      Modelica.SIunits.Force F_f                                  "Сила трения";
      Modelica.SIunits.Momentum M                                 "Количество движения воды";
      Modelica.SIunits.Pressure p_i                               "Входное давление";
      Modelica.SIunits.Pressure p_o                               "Выходное давление";
      Modelica.SIunits.VolumeFlowRate Vdot(start = Vdot_0)        "Объёмный расход";
    
      extends HydroPowerPlant.Connectors.IWater.MeticulousTwoPortWaterConnector;
    
    equation
      
      Vdot = mdot / data.rho;             // Объёмный расход в трубе

      v = Vdot / A_;                      // Скорость воды

      M = data.rho * L * Vdot;            // Количество движения и масса водного потока
      m = data.rho * A_ * L;

      F_f = Utilities.Darcy_Weisbach.DarcyFrictionFactor(v, D_, L, data.rho, data.mu, p_eps);                                 // Определяем коэффициент трения Дарси
    
      der(M) = data.rho * Vdot ^ 2 * (1 / A_i - 1 / A_o) + p_i * A_i - p_o * A_o - F_f + m * data.g * cos_theta;

      p_i = i.p;                          // Давление водного потока
      p_o = o.p;
      
    end Pipeline;
  end Watercourse;

  package ElectricalInstallations
    package HydraulicTurbines
      extends Modelica.Icons.GearIcon;

      model GeneralTurbine "Модель простой гидротурбины"
      
        outer Presets data;
        
        parameter Real C_v = 3.7                                "Пропускная способность клапана";
        parameter Modelica.SIunits.Height H_n = 460             "Номинальный напор";
        parameter Modelica.SIunits.VolumeFlowRate Vdot_n = 23.4 "Номинальный массовый расход";
        parameter Modelica.SIunits.PerUnit u_n = 0.95           "Номинальная пропускная способность клапана";
        parameter Modelica.SIunits.Efficiency eta_h = 0.9       "Гидравлическая эффективность турбины";
        
        extends PowerHouse.TorqueGenerator(power(y = Wdot_s));
        extends Connectors.ITurbine.SmartTurbineConnector;
        
        Modelica.SIunits.Pressure dp                            "Перепад давления";
        Modelica.SIunits.EnergyFlowRate Kdot_i_tr               "Поток кинетической энергии";
        Modelica.SIunits.VolumeFlowRate Vdot                    "Объёмный расход";
        
        output Modelica.SIunits.EnergyFlowRate Wdot_s           "Мощность вала турбины";
        Modelica.Blocks.Math.Feedback lossCorrection            "Блок для коррекции потерь от трения";
      
      equation
      
        Vdot = mdot / (data.rho * (1 + data.beta * (i.p - data.p_a))) "Проверка на сжимаемость воды";
        
        dp = Vdot ^ 2 * data.p_a / (C_v * u_t) ^ 2                    "Уравнение перепада давления для турбинного клапана";
        
        dp = i.p - o.p                                                "Передаём информацию о перепадах давления на элементы соединения";
        
        Kdot_i_tr = dp * Vdot                                         "Уравнение сохранения энергии";
       
        Wdot_s = eta_h * Kdot_i_tr;
        
        connect(P_out, lossCorrection.y);
        connect(lossCorrection.u1, power.y);
        connect(frictionLoss.power, lossCorrection.u2);
        
      end GeneralTurbine;
    end HydraulicTurbines;

    package PowerHouse
      package Generators
        model CommonGenerator "Модель простого генератора"
        
          extends TorqueGenerator(power(y = -Pload));
        
          Modelica.Blocks.Interfaces.RealInput Pload(unit = "W") "Нагрузка электросети";

        end CommonGenerator;
      end Generators;

      partial model TorqueGenerator
      
        outer Presets data;
        import Modelica.Constants.pi;
      
        parameter Modelica.SIunits.Power Pmax = 100e6                           "Максимальная номинальная мощность (для ограничения крутящего момента и расчета H)";
        parameter Modelica.SIunits.Time H = 2.75                                "Постоянная инерции H. Как правило, 2с (гидроусилители с высоким напором) или 6с (газовые или гидроусилители с низким напором)";
        parameter Modelica.SIunits.MomentOfInertia J = 2e5                      "Момент инерции блока";
        parameter Integer p(min = 2) = 12                                       "Количество полюсов";
        parameter Modelica.SIunits.Power Ploss = 0                              "Потери от трения блока при номинальной скорости";
        parameter Modelica.SIunits.AngularVelocity w_0 = data.f_0 * 4 * pi / p  "Начальная механическая угловая скорость";
      
        Modelica.Blocks.Math.Division power2torque;
        Modelica.Mechanics.Rotational.Sensors.SpeedSensor speedSensor;
        Modelica.Mechanics.Rotational.Components.Inertia inertia(J = J,  w(start = w_0));
        
        Modelica.Blocks.Interfaces.RealOutput f(unit = "Hz");
        Modelica.Electrical.Machines.Losses.Friction friction(frictionParameters(PRef = Ploss, wRef = data.f_0 * 4 * pi / p));
        Modelica.Mechanics.Rotational.Components.Fixed fixed;
        Modelica.Mechanics.Rotational.Sources.Torque torque;
        Modelica.Blocks.Nonlinear.Limiter div0protect(uMax = Modelica.Constants.inf, uMin = Modelica.Constants.small);
        Modelica.Blocks.Math.Gain toHz(k = Modelica.SIunits.Conversions.to_Hz(p / 2));
        Modelica.Blocks.Nonlinear.Limiter torqueLimit(uMax = Pmax / w_0);
        Modelica.Blocks.Interfaces.RealOutput w                                  "Угловая скорость генератора";
      
        Modelica.Mechanics.Rotational.Interfaces.Flange_b flange(phi(start = 0, fixed = true)) "Гребень правого вала";
        Modelica.Blocks.Sources.RealExpression power;
        Modelica.Mechanics.Rotational.Sensors.PowerSensor frictionLoss;
        Modelica.Mechanics.Rotational.Sources.ConstantSpeed nomSpeed(w_fixed = w_0 * p / 2);
        Modelica.Mechanics.Rotational.Components.IdealGear toSysSpeed(ratio = 2 / p) "Дифференциал с преобразованием скорости p = 2";
        
      equation
        
        connect(toHz.u, speedSensor.w);
        
        connect(w, speedSensor.w);
        
        connect(div0protect.y, power2torque.u2);
        
        connect(f , toHz.y);
        
        connect(power2torque.y, torqueLimit.u);
        
        connect(torqueLimit.y, torque.tau);
        
        connect(speedSensor.w, div0protect.u);
        
        connect(inertia.flange_b, speedSensor.flange);
        
        connect(friction.support, fixed.flange);
        
        connect(torque.flange, inertia.flange_a);
        
        connect(w, w);
        
        connect(power.y, power2torque.u1);
        
        connect(frictionLoss.flange_a, inertia.flange_b);
        
        connect(frictionLoss.flange_b, friction.flange);
        
        connect(nomSpeed.flange, flange);
        
        connect(flange, toSysSpeed.flange_b);
        
        connect(toSysSpeed.flange_a, inertia.flange_b);
      
      end TorqueGenerator;
    end PowerHouse;
  end ElectricalInstallations;

  package WorkingModels
    extends Modelica.Icons.ExamplesPackage;

    model BasicExample "Базовая модель плотинной ГЭС"
    
      extends Modelica.Icons.Example;
      inner HydroPowerPlant.Presets data;
    
      HydroPowerPlant.WaterSources.Reservoir Reservoir(H_0 = 60);
      HydroPowerPlant.Watercourse.Pipeline IntakeStructure(H = 30);
      HydroPowerPlant.Watercourse.SurgeChamber SurgeChamber(h_0 = 8);
      HydroPowerPlant.Watercourse.Pipeline Penstock(D_i = 3, D_o = 3, H = 4, L = 60);
      HydroPowerPlant.ElectricalInstallations.HydraulicTurbines.GeneralTurbine Turbine;
      Modelica.Blocks.Sources.Ramp Control(duration = 1, height = -0.05, offset = 1, startTime = 600);
      HydroPowerPlant.Watercourse.Pipeline Discharge(H = 0.5, L = 60);
      HydroPowerPlant.WaterSources.Reservoir Tailrace(H_0 = 5);
    
    equation
    
      connect(Reservoir.w_o, IntakeStructure.i);
      connect(IntakeStructure.o, SurgeChamber.i);
      connect(SurgeChamber.o, Penstock.i);
      connect(Penstock.o, Turbine.i);
      connect(Control.y, Turbine.u_t);
      connect(Turbine.o, Discharge.i);
      connect(Discharge.o, Tailrace.w_o);
    
    end BasicExample;
  end WorkingModels;
  annotation(
    version = "0.9.0",
    versionDate = "17/04/2022",
    uses(Modelica(version = "3.2.3")),
    preferredView = "info",
    Documentation(info = "<html>
      <p>HydroPowerPlant is a software package available for using in the OpenModelica Connection Editor development environment and used to simulate various mechanical and electrical systems and assemblies of hydroelectric power plants.</p>
      <p>The original <a href=\"https://github.com/CourteousSleet/Diploma\">repository</a>.</p>
</html>"));
end HydroPowerPlant;

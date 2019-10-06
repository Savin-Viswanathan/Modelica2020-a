package OceanEngineering
  extends Modelica.Icons.Package;

  package Components
    package Waves
      package RegularWave
        model Regular_Airy_Wave
          extends Modelica.Blocks.Icons.Block;
          OceanEngineering.Connectors.WaveDataConnector wdoc(n_omega_i = 1) annotation(
            Placement(transformation(extent = {{80, -20}, {120, 20}})));
          constant Real pi = Modelica.Constants.pi "Value of pi";
          constant Real g = Modelica.Constants.g_n "Acceleration due to gravity";
          parameter Modelica.SIunits.Length d = 100 "Water depth";
          parameter Modelica.SIunits.Density rho_w = 1025 "Density of Water";
          parameter Modelica.SIunits.Length Hr = 1 "Wave Height";
          parameter Real Tr = 7 "Wave Period";
          parameter Real Trmp = 20 "Wave ramp time";
          parameter Real Tdel = 0 "Wave delay";
          parameter Integer n_omega_i = 1 "Number of frequency components";
          Real T[size(omega, 1)] "Wave period";
          Modelica.SIunits.AngularFrequency omega[n_omega_i] "Frequency components of time series";
          Modelica.SIunits.Length zeta0i[n_omega_i] "Amplitude of wave component";
          Modelica.SIunits.Length SSE_X0 "Sea surface elevation at x=0 - to plot generated waves";
          Real epsilon[n_omega_i] "Random phase of wave components";
          Real k[n_omega_i](each unit = "m^-1") = OceanEngineering.Functions.waveNumberIterator(d, omega) "Wave number of component waves";
        equation
          for i in 1:n_omega_i loop
            omega[i] = 2 * pi / Tr;
            epsilon[i] = 0;
            if time < Tdel then
              zeta0i[i] = 0;
            elseif time < Tdel + Trmp then
              zeta0i[i] = sin(pi / 2 * (time - Tdel) / Trmp) * Hr / 2;
            else
              zeta0i[i] = Hr / 2;
            end if;
            T[i] = Tr;
          end for;
          SSE_X0 = sum(zeta0i .* cos(omega * time - 2 * pi * epsilon));
          wdoc.d = d;
          wdoc.rho_w = rho_w;
          wdoc.omega = omega;
          wdoc.T = T;
          wdoc.k = k;
          wdoc.epsilon = epsilon;
          wdoc.zeta0i = zeta0i;
          wdoc.SSE_X0 = SSE_X0;
          annotation(
            Icon(graphics = {Text(origin = {-54, 65}, extent = {{-40, 17}, {40, -17}}, textString = "H,T"), Rectangle(origin = {100, 0}, fillColor = {170, 255, 0}, fillPattern = FillPattern.Solid, extent = {{-20, -20}, {20, 20}}), Line(origin = {-15, 0}, points = {{-81, 0}, {81, 0}}), Line(origin = {-14.97, 1.01}, points = {{-81.0281, -1.00867}, {-33.0281, 42.9913}, {36.9719, -43.0087}, {80.9719, -1.00867}}, color = {0, 0, 255}, thickness = 0.5, smooth = Smooth.Bezier)}, coordinateSystem(initialScale = 0.1)),
            experiment(StartTime = 0, StopTime = 400, Tolerance = 1e-06, Interval = 0.5));
        end Regular_Airy_Wave;
      end RegularWave;

      package IrregularWave
        model IRW_PM_RDFCWI
          extends Modelica.Blocks.Icons.Block;
          OceanEngineering.Connectors.WaveDataConnector wdoc(n_omega_i = 100) annotation(
            Placement(transformation(extent = {{80, -20}, {120, 20}})));
          constant Real pi = Modelica.Constants.pi "Value of pi";
          constant Real g = Modelica.Constants.g_n "Acceleration due to gravity";
          parameter Modelica.SIunits.Length d = 100 "Water depth";
          parameter Modelica.SIunits.Density rho_w = 1025 "Density of Water";
          parameter Modelica.SIunits.Length Hs = 1 "Significant Wave Height";
          parameter Modelica.SIunits.AngularFrequency omega_min = 0.03141 "Lowest frequency component/frequency interval";
          parameter Modelica.SIunits.AngularFrequency omega_max = 3.141 "Highest frequency component";
          parameter Integer n_omega_i = 100 "Number of frequency components";
          parameter Integer localSeed = 614657 "Local seed for random number generator";
          parameter Integer globalSeed = 30020 "Global seed for random number generator";
          parameter Real rnd_shft[n_omega_i] = OceanEngineering.Functions.randomNumberGenerator(localSeed, globalSeed, n_omega_i);
          parameter Integer localSeed1 = 614757 "Local seed for random number generator";
          parameter Integer globalSeed1 = 40020 "Global seed for random number generator";
          parameter Real epsilon[n_omega_i] = OceanEngineering.Functions.randomNumberGenerator(localSeed1, globalSeed1, n_omega_i);
          parameter Integer SSE_X0on = 1 "Flag for generating SSE profile";
          parameter Integer Tdel = 10 "Time delay before waves are generated";
          parameter Real Trmp = 10 "Interval for ramping up of waves during start phase";
          parameter Real omega[n_omega_i] = OceanEngineering.Functions.frequencySelector(omega_min, omega_max, rnd_shft);
          parameter Real S[n_omega_i] = OceanEngineering.Functions.spectrumGenerator_PM(Hs, omega);
          parameter Modelica.SIunits.Length zeta0i[n_omega_i] = sqrt(2 * S * omega_min) "Amplitude of wave component";
          parameter Real T[n_omega_i] = 2 * pi ./ omega "Wave period";
          parameter Real k[n_omega_i] = OceanEngineering.Functions.waveNumberIterator(d, omega) "Wave number";
          Modelica.SIunits.Length SSE_X0 "Sea surface elevation at x=0 - to plot generated waves";
          Real zeta0i_rmp[n_omega_i] "Ramp value of zeta0i";
        equation
          for i in 1:n_omega_i loop
            if time < Tdel then
              zeta0i_rmp[i] = 0;
            elseif time < Tdel + Trmp then
              zeta0i_rmp[i] = sin(pi / 2 * (time - Tdel) / Trmp) * zeta0i[i];
            else
              zeta0i_rmp[i] = zeta0i[i];
            end if;
          end for;
          if SSE_X0on == 1 then
            SSE_X0 = sum(zeta0i_rmp .* cos(omega * time - 2 * pi * epsilon));
          else
            SSE_X0 = 0;
          end if;
          wdoc.d = d;
          wdoc.rho_w = rho_w;
          wdoc.omega = omega;
          wdoc.T = T;
          wdoc.k = k;
          wdoc.epsilon = epsilon;
          wdoc.zeta0i = zeta0i_rmp;
          wdoc.SSE_X0 = SSE_X0;
          annotation(
            Icon(graphics = {Line(origin = {-50.91, 48.08}, points = {{-33.2809, -22.5599}, {-21.2809, -20.5599}, {-13.2809, 27.4401}, {6.71907, -20.5599}, {24.7191, -24.5599}, {42.7191, -24.5599}, {44.7191, -24.5599}}, color = {255, 0, 0}, smooth = Smooth.Bezier), Line(origin = {-37, 51}, points = {{-51, 29}, {-51, -29}, {37, -29}}), Text(origin = {6, 55}, extent = {{-40, 17}, {40, -17}}, textString = "Hs"), Line(origin = {22, 4}, points = {{0, 22}, {0, -22}}, thickness = 1, arrow = {Arrow.None, Arrow.Filled}), Line(origin = {-7.57, -61.12}, points = {{-82.4341, -12.8774}, {-76.4341, -2.87735}, {-72.4341, -6.87735}, {-62.4341, 13.1226}, {-50.4341, -26.8774}, {-46.4341, -20.8774}, {-38.4341, -26.8774}, {-34.4341, -18.8774}, {-34.4341, 3.12265}, {-26.4341, 1.12265}, {-20.4341, 7.12265}, {-12.4341, 9.12265}, {-8.43408, 19.1226}, {1.56592, -4.87735}, {7.56592, -24.8774}, {19.5659, -6.87735}, {21.5659, 9.12265}, {31.5659, 13.1226}, {39.5659, -0.87735}, {43.5659, 11.1226}, {55.5659, 15.1226}, {63.5659, 27.1226}, {79.5659, -22.8774}}, color = {0, 0, 255}, smooth = Smooth.Bezier), Rectangle(origin = {100, 0}, fillColor = {85, 255, 127}, fillPattern = FillPattern.Solid, extent = {{-20, 20}, {20, -20}})}, coordinateSystem(initialScale = 0.1)),
            experiment(StartTime = 0, StopTime = 400, Tolerance = 1e-06, Interval = 0.2));
        end IRW_PM_RDFCWI;
      end IrregularWave;
    end Waves;

    package CurrentProfiles
      block CurrentProfile_4pt
        extends Modelica.Blocks.Icons.Block;
        OceanEngineering.Connectors.CurrentDataConnector cdc annotation(
          Placement(transformation(extent = {{80, -20}, {120, 20}})));
        constant Real pi = Modelica.Constants.pi "Value of pi";
        parameter Real zcg[:] = {-100, -20, -10, 0};
        parameter Real Uf[:] = {0, 0.5, 1, 2};
        parameter Real Trmp = 5;
        Real Ucg[size(Uf, 1)];
      equation
        for i in 1:size(Uf, 1) loop
          if time < Trmp then
            Ucg[i] = sin(pi / 2 * time / Trmp) * Uf[i];
          else
            Ucg[i] = Uf[i];
          end if;
        end for;
        cdc.zcg = zcg;
        cdc.Ucg = Ucg;
        annotation(
          Icon(graphics = {Rectangle(origin = {0, -30}, fillColor = {0, 170, 255}, fillPattern = FillPattern.Solid, extent = {{-100, 70}, {100, -70}}), Line(origin = {-80, -29}, points = {{0, 69}, {0, -69}}), Line(origin = {-50.0068, -29.0068}, points = {{62.0068, 69.0068}, {10.0068, 17.0068}, {-9.99322, -58.9932}, {-29.9932, -68.9932}}), Line(origin = {-53, 16.0755}, points = {{-27, 0}, {41, 0}}, arrow = {Arrow.None, Arrow.Filled}), Line(origin = {-60, -14}, points = {{-20, 0}, {20, 0}}, arrow = {Arrow.None, Arrow.Filled}), Line(origin = {-65, -48}, points = {{-15, 0}, {15, 0}}, arrow = {Arrow.None, Arrow.Filled}), Line(origin = {-70, -86}, points = {{-10, 0}, {10, 0}}, arrow = {Arrow.None, Arrow.Filled}), Line(origin = {-47.6886, 40.3302}, points = {{-32, 0}, {60, 0}}, arrow = {Arrow.None, Arrow.Filled}), Rectangle(origin = {100, 0}, fillColor = {255, 255, 0}, fillPattern = FillPattern.Solid, extent = {{-20, -20}, {20, 20}})}, coordinateSystem(initialScale = 0.1)),
          experiment(StartTime = 0, StopTime = 10, Tolerance = 1e-06, Interval = 0.5));
      end CurrentProfile_4pt;
    end CurrentProfiles;
  end Components;

  package Functions
    function waveNumberIterator
      input Real d "Water depth";
      input Real omega[:] "Wave frequency";
      output Real k[size(omega, 1)] "Wave number";
    protected
      constant Real g = Modelica.Constants.g_n;
      constant Real pi = Modelica.Constants.pi;
      parameter Integer n = size(omega, 1);
      Real T[size(omega, 1)] "Wave period";
      Real L0[size(omega, 1)] "Deepwater wavelength";
      Real Ll(start = 0, fixed = true) "Intermediate loop value";
      Real Llc(start = 0, fixed = true) "Intermediate loop comparator value";
      Real L[size(omega, 1)] "Iterated wave length";
    algorithm
      T := 2 * pi ./ omega;
      L0 := g * T .^ 2 / (2 * pi);
      for i in 1:size(omega, 1) loop
        Ll := L0[i];
        Llc := 0;
        while abs(Llc - Ll) > 0.001 loop
          Llc := Ll;
          L[i] := g * T[i] ^ 2 / (2 * pi) * tanh(2 * pi / Ll * d);
          Ll := L[i];
        end while;
      end for;
      k := 2 * pi ./ L;
    end waveNumberIterator;

    function randomNumberGenerator
      input Integer ls = 614657;
      input Integer gs = 30020;
      input Integer n = 100;
      output Real r64[n];
    protected
      Integer state64[2](each start = 0, each fixed = true);
    algorithm
      state64[1] := 0;
      state64[2] := 0;
      for i in 1:n loop
        if i == 1 then
          state64 := Modelica.Math.Random.Generators.Xorshift64star.initialState(ls, gs);
          r64[i] := 0;
        else
          (r64[i], state64) := Modelica.Math.Random.Generators.Xorshift64star.random(pre(state64));
        end if;
      end for;
    end randomNumberGenerator;

    function frequencySelector
      input Real omega_min;
      input Real omega_max;
      input Real epsilon[:];
      output Real omega[size(epsilon, 1)];
    protected
      parameter Real delta_omega = (omega_max - omega_min) / (size(epsilon, 1) - 1);
      parameter Real ref_omega[size(epsilon, 1)] = omega_min:delta_omega:omega_max;
    algorithm
      omega[1] := omega_min;
      for i in 2:size(epsilon, 1) - 1 loop
        omega[i] := ref_omega[i] + epsilon[i] * omega_min;
      end for;
      omega[size(epsilon, 1)] := omega_max;
    end frequencySelector;

    function spectrumGenerator_PM
      input Real Hs = 2 "Significant wave height";
      input Real omega[:] "Frequency components";
      output Real spec[size(omega, 1)] "Spectral values for input frequencies";
    protected
      constant Real pi = Modelica.Constants.pi;
      constant Real g = Modelica.Constants.g_n;
    algorithm
      for i in 1:size(omega, 1) loop
        spec[i] := 0.0081 * g ^ 2 / omega[i] ^ 5 * exp(-0.032 * (g / (Hs * omega[i] ^ 2)) ^ 2);
      end for;
    end spectrumGenerator_PM;
  end Functions;

  package Connectors
    connector WaveDataConnector
      parameter Integer n_omega_i = 100 "number of wave components";
      Modelica.Blocks.Interfaces.RealOutput d;
      Modelica.Blocks.Interfaces.RealOutput rho_w;
      Modelica.Blocks.Interfaces.RealOutput omega[n_omega_i];
      Modelica.Blocks.Interfaces.RealOutput T[n_omega_i];
      Modelica.Blocks.Interfaces.RealOutput k[n_omega_i];
      Modelica.Blocks.Interfaces.RealOutput epsilon[n_omega_i];
      Modelica.Blocks.Interfaces.RealOutput zeta0i[n_omega_i];
      Modelica.Blocks.Interfaces.RealOutput SSE_X0;
      annotation(
        Icon(graphics = {Rectangle(fillColor = {170, 255, 0}, fillPattern = FillPattern.Solid, extent = {{-100, 100}, {100, -100}}), Text(origin = {-1, 3}, extent = {{-67, 41}, {67, -41}}, textString = "%name")}, coordinateSystem(initialScale = 0.1)));
    end WaveDataConnector;

    expandable connector EnvironmentBus
      annotation(
        Icon(graphics = {Rectangle(fillColor = {170, 170, 255}, fillPattern = FillPattern.Solid, extent = {{-100, 100}, {100, -100}}), Line(origin = {0, 80}, points = {{-80, 0}, {80, 0}}), Line(origin = {-0.566026, -0.471698}, points = {{-80, 0}, {80, 0}}), Line(origin = {-0.471698, 39.9057}, points = {{-80, 0}, {80, 0}}), Line(origin = {0.51888, -40}, points = {{-80, 0}, {80, 0}}), Line(origin = {-0.707535, -79.8585}, points = {{-80, 0}, {80, 0}}), Rectangle(origin = {-2, -1}, fillColor = {255, 0, 127}, fillPattern = FillPattern.Solid, extent = {{-52, 91}, {52, -91}}), Rectangle(origin = {-100, 0}, fillColor = {85, 85, 127}, fillPattern = FillPattern.Solid, extent = {{-20, 20}, {20, -20}}), Line(origin = {-86, 50}, points = {{6, 30}, {-6, -30}}), Line(origin = {-83.1673, 29.8327}, points = {{2.66543, 10}, {-3.33457, -10}}), Line(origin = {-88.171, -50}, points = {{7, -30}, {-5, 30}}), Line(origin = {-84, -30}, points = {{4, -10}, {-2, 10}}), Rectangle(origin = {100, 0}, fillColor = {85, 85, 255}, fillPattern = FillPattern.Solid, extent = {{-20, -20}, {20, 20}}), Line(origin = {87.8284, 50}, points = {{-7, 30}, {7, -30}}), Line(origin = {83, 30}, points = {{-3, 10}, {3, -10}}), Line(origin = {84, -30}, points = {{-4, -10}, {4, 10}}), Line(origin = {87.8284, -50}, points = {{-9, -30}, {7, 30}}), Text(origin = {1, 113}, lineColor = {0, 0, 255}, extent = {{-150, 60}, {150, -10}}, textString = "%name")}, coordinateSystem(initialScale = 0.1)));
    end EnvironmentBus;

    connector CurrentDataConnector
      Modelica.Blocks.Interfaces.RealOutput zcg[4];
      Modelica.Blocks.Interfaces.RealOutput Ucg[4];
      annotation(
        Icon(graphics = {Rectangle(fillColor = {170, 255, 0}, fillPattern = FillPattern.Solid, extent = {{-100, 100}, {100, -100}}), Text(origin = {-1, 3}, extent = {{-67, 41}, {67, -41}}, textString = "%name")}, coordinateSystem(initialScale = 0.1)));
    end CurrentDataConnector;
  end Connectors;

  package SampleSimulations
    extends Modelica.Icons.ExamplesPackage;

    package PropertyChk_RegularWave
      model RegWavePropProbe
        extends Modelica.Blocks.Icons.Block;
        OceanEngineering.SampleSimulations.PropertyChk_RegularWave.GeneralDataConnector gdc(n_omega_i = 1) annotation(
          Placement(visible = true, transformation(extent = {{-80, -20}, {-120, 20}}, rotation = 0), iconTransformation(extent = {{-118, -20}, {-78, 20}}, rotation = 0)));
        // Value of pi
        constant Real pi = Modelica.Constants.pi;
        // Value of acceleration due to gravity
        constant Real g = Modelica.Constants.g_n;
        // x coordinates for plotting the progressive wave
        parameter Real x[:] = 0:0.1:30;
        //z coordinates for plotting water particle trajectories and velocities
        parameter Real z[:] = (-10):1:0;
        //x coordinates to be used in plot for water particle trajectories
        parameter Real x_del[:] = (-6):0.1:6;
        // x coordinates to be used in plot for  water particle velocities
        parameter Real x_lng[:] = (-20):1:40;
        // x coordinates to be used in plot showing waves of different wave lengths
        parameter Real x_lng_lng[:] = (-110):1:110;
        //z coordinates for plotting pressure
        parameter Real z_pr[:] = (-10):0.1:1;
        //Wave Length
        Real L;
        //Locations at which water particle velocities and displacements are to be calculated
        Real x0;
        Real x1;
        Real x2;
        Real x3;
        Real x4;
        //Water particle velocity components at above locations
        Real u0[size(z, 1)];
        Real w0[size(z, 1)];
        Real u1[size(z, 1)];
        Real w1[size(z, 1)];
        Real u2[size(z, 1)];
        Real w2[size(z, 1)];
        Real u3[size(z, 1)];
        Real w3[size(z, 1)];
        Real u4[size(z, 1)];
        Real w4[size(z, 1)];
        //Hydrodstatic and hydrodynamic pressures at x0,x1
        Real p_hydx0[size(z_pr, 1)];
        Real p_dynx0[size(z_pr, 1)];
        Real p_hydx1[size(z_pr, 1)];
        Real p_dynx1[size(z_pr, 1)];
        Real p_hydx2[size(z_pr, 1)];
        Real p_dynx2[size(z_pr, 1)];
        //Sea surface elevation based on x[:]
        Real eta[size(x, 1)];
        //Water particle displacements at different z[:]
        Real del_x[size(z, 1)];
        Real del_z[size(z, 1)];
        //Water particle positions at different z[:]
        Real X[size(z, 1)];
        Real Z[size(z, 1)];
        //Sea Surface elevation corresponding to x_del[:],x_lng[:] and x_lng_lng[:]
        Real eta_del[size(x_del, 1)];
        Real eta_lng[size(x_lng, 1)];
        Real eta_lng_lng[size(x_lng_lng, 1)];
        //Sea surface elevations at x0,x1,x2
        Real eta_x0;
        Real eta_x1;
        Real eta_x2;
        //Total pressure at x0,x1,x2 for at different z_pr[:]
        Real totp_x0[size(z_pr, 1)];
        Real totp_x1[size(z_pr, 1)];
        Real totp_x2[size(z_pr, 1)];
      equation
// Determine the wave profile at points defined in x[:]
        for i in 1:size(x, 1) loop
          eta[i] = gdc.zeta0i[1] * cos(gdc.k[1] * x[i] - gdc.omega[1] * time);
        end for;
// Determine locations to measure wave properties as x0,x1,x2,x3,x4
// Calculate Wave length
        L = 2 * pi / gdc.k[1];
// Coordinate of x-location-1
        x0 = 0;
// Coordinate of x-location-2
        x1 = L / 4;
// Coordinate of x-location-3
        x2 = L / 2;
// Coordinate of x-location-4
        x3 = 3 * L / 4;
// Coordinate of x-location-5
        x4 = L / 3;
// Determine components of wave induced water-particle velocities at x0 to x4
        u0 = pi * gdc.Hr / gdc.Tr * cosh(gdc.k[1] * (z .+ gdc.d)) / cosh(gdc.k[1] * gdc.d) * cos(gdc.k[1] * x0);
        w0 = pi * gdc.Hr / gdc.Tr * sinh(gdc.k[1] * (z .+ gdc.d)) / sinh(gdc.k[1] * gdc.d) * sin(gdc.k[1] * x0);
        u1 = pi * gdc.Hr / gdc.Tr * cosh(gdc.k[1] * (z .+ gdc.d)) / cosh(gdc.k[1] * gdc.d) * cos(gdc.k[1] * x1);
        w1 = pi * gdc.Hr / gdc.Tr * sinh(gdc.k[1] * (z .+ gdc.d)) / sinh(gdc.k[1] * gdc.d) * sin(gdc.k[1] * x1);
        u2 = pi * gdc.Hr / gdc.Tr * cosh(gdc.k[1] * (z .+ gdc.d)) / cosh(gdc.k[1] * gdc.d) * cos(gdc.k[1] * x2);
        w2 = pi * gdc.Hr / gdc.Tr * sinh(gdc.k[1] * (z .+ gdc.d)) / sinh(gdc.k[1] * gdc.d) * sin(gdc.k[1] * x2);
        u3 = pi * gdc.Hr / gdc.Tr * cosh(gdc.k[1] * (z .+ gdc.d)) / cosh(gdc.k[1] * gdc.d) * cos(gdc.k[1] * x3);
        w3 = pi * gdc.Hr / gdc.Tr * sinh(gdc.k[1] * (z .+ gdc.d)) / sinh(gdc.k[1] * gdc.d) * sin(gdc.k[1] * x3);
        u4 = pi * gdc.Hr / gdc.Tr * cosh(gdc.k[1] * (z .+ gdc.d)) / cosh(gdc.k[1] * gdc.d) * cos(gdc.k[1] * x4);
        w4 = pi * gdc.Hr / gdc.Tr * sinh(gdc.k[1] * (z .+ gdc.d)) / sinh(gdc.k[1] * gdc.d) * sin(gdc.k[1] * x4);
// Determine water particle displacements at x0
        del_x = -gdc.Hr / 2 * cosh(gdc.k[1] * (z .+ gdc.d * ones(size(z, 1)))) ./ sinh(gdc.k[1] * gdc.d) * sin(-gdc.omega[1] * time);
        del_z = gdc.Hr / 2 * sinh(gdc.k[1] * (z .+ gdc.d * ones(size(z, 1)))) ./ sinh(gdc.k[1] * gdc.d) * cos(-gdc.omega[1] * time);
// Determine water particle positions
        X = del_x;
        Z = z + del_z;
// Determine SSE using different x[:] for plotting
        eta_del = gdc.Hr / 2 * cos(gdc.k[1] * x_del - gdc.omega[1] * time * ones(size(x_del, 1)));
        eta_lng = gdc.Hr / 2 * cos(gdc.k[1] * x_lng - gdc.omega[1] * time * ones(size(x_lng, 1)));
        eta_lng_lng = gdc.Hr / 2 * cos(gdc.k[1] * x_lng_lng - gdc.omega[1] * time * ones(size(x_lng_lng, 1)));
//Determine SSE at wave crest, wave downcrossing and wave trough for use in pressure calculations
        eta_x0 = gdc.Hr / 2 * cos(gdc.k[1] * x0);
        eta_x1 = gdc.Hr / 2 * cos(gdc.k[1] * x1);
        eta_x2 = gdc.Hr / 2 * cos(gdc.k[1] * x2);
//Determine static and dynamic components of pressure for different z_pr[:] at x0,x1 and x2
        for i in 1:size(z_pr, 1) loop
          if eta_x0 < 0 then
            if z_pr[i] < eta_x0 then
              p_hydx0[i] = -gdc.rho_w * g * z_pr[i];
              p_dynx0[i] = gdc.rho_w * g * eta_x0 * cosh(gdc.k[1] * (z_pr[i] + gdc.d)) / cosh(gdc.k[1] * gdc.d);
            elseif z_pr[i] < 0 then
              p_hydx0[i] = -gdc.rho_w * g * z_pr[i];
              p_dynx0[i] = gdc.rho_w * g * z_pr[i];
            else
              p_hydx0[i] = 0;
              p_dynx0[i] = 0;
            end if;
          else
            if z_pr[i] < 0 then
              p_hydx0[i] = -gdc.rho_w * g * z_pr[i];
              p_dynx0[i] = gdc.rho_w * g * eta_x0 * cosh(gdc.k[1] * (z_pr[i] + gdc.d)) / cosh(gdc.k[1] * gdc.d);
            elseif z_pr[i] < eta_x0 then
              p_hydx0[i] = 0;
              p_dynx0[i] = gdc.rho_w * g * (eta_x0 - z_pr[i]);
            else
              p_hydx0[i] = 0;
              p_dynx0[i] = 0;
            end if;
          end if;
          if eta_x1 < 0 then
            if z_pr[i] < eta_x1 then
              p_hydx1[i] = -gdc.rho_w * g * z_pr[i];
              p_dynx1[i] = gdc.rho_w * g * eta_x1 * cosh(gdc.k[1] * (z_pr[i] + gdc.d)) / cosh(gdc.k[1] * gdc.d);
            elseif z_pr[i] < 0 then
              p_hydx1[i] = -gdc.rho_w * g * z_pr[i];
              p_dynx1[i] = gdc.rho_w * g * z_pr[i];
            else
              p_hydx1[i] = 0;
              p_dynx1[i] = 0;
            end if;
          else
            if z_pr[i] < 0 then
              p_hydx1[i] = -gdc.rho_w * g * z_pr[i];
              p_dynx1[i] = gdc.rho_w * g * eta_x1 * cosh(gdc.k[1] * (z_pr[i] + gdc.d)) / cosh(gdc.k[1] * gdc.d);
            elseif z_pr[i] < eta_x1 then
              p_hydx1[i] = 0;
              p_dynx1[i] = gdc.rho_w * g * (eta_x1 - z_pr[i]);
            else
              p_hydx1[i] = 0;
              p_dynx1[i] = 0;
            end if;
          end if;
          if eta_x2 < 0 then
            if z_pr[i] < eta_x2 then
              p_hydx2[i] = -gdc.rho_w * g * z_pr[i];
              p_dynx2[i] = gdc.rho_w * g * eta_x2 * cosh(gdc.k[1] * (z_pr[i] + gdc.d)) / cosh(gdc.k[1] * gdc.d);
            elseif z_pr[i] < 0 then
              p_hydx2[i] = -gdc.rho_w * g * z_pr[i];
              p_dynx2[i] = gdc.rho_w * g * z_pr[i];
            else
              p_hydx2[i] = 0;
              p_dynx2[i] = 0;
            end if;
          else
            if z_pr[i] < 0 then
              p_hydx2[i] = -gdc.rho_w * g * z_pr[i];
              p_dynx2[i] = gdc.rho_w * g * eta_x2 * cosh(gdc.k[1] * (z_pr[i] + gdc.d)) / cosh(gdc.k[1] * gdc.d);
            elseif z_pr[i] < eta_x2 then
              p_hydx2[i] = 0;
              p_dynx2[i] = gdc.rho_w * g * (eta_x2 - z_pr[i]);
            else
              p_hydx2[i] = 0;
              p_dynx2[i] = 0;
            end if;
          end if;
        end for;
//Calculate total pressures at different z_pr[:]
        totp_x0 = p_hydx0 + p_dynx0;
        totp_x1 = p_hydx1 + p_dynx1;
        totp_x2 = p_hydx2 + p_dynx2;
        annotation(
          Icon(graphics = {Ellipse(origin = {-9, 12}, fillColor = {170, 255, 255}, fillPattern = FillPattern.Sphere, lineThickness = 2, extent = {{-35, 34}, {35, -34}}, endAngle = 360), Line(origin = {43.04, -41.78}, points = {{-27, 27}, {27, -27}, {27, -27}}, thickness = 2.5), Line(origin = {52, 64}, points = {{-18, 16}, {18, -16}, {18, -16}}, color = {255, 0, 0}, thickness = 1), Line(origin = {51, 66}, points = {{19, 18}, {-19, -18}}, color = {255, 0, 0}, thickness = 1), Line(origin = {-60, -62.65}, points = {{-27.9998, -5.35334}, {-11.9998, -21.3533}, {28.0002, 20.6467}}, color = {0, 255, 127}, thickness = 1)}, coordinateSystem(initialScale = 0.1)));
      end RegWavePropProbe;

      model Check_RegularWave
        Connectors.EnvironmentBus environmentBus1 annotation(
          Placement(visible = true, transformation(origin = {0, 0}, extent = {{-10, -10}, {10, 10}}, rotation = 0), iconTransformation(origin = {0, 0}, extent = {{-10, -10}, {10, 10}}, rotation = 0)));
        Regular_Airy_Wave_Test regular_Airy_Wave_Test1(Hr = 1, Tr = 7) annotation(
          Placement(visible = true, transformation(origin = {-50, 0}, extent = {{-10, -10}, {10, 10}}, rotation = 0)));
        OceanEngineering.SampleSimulations.PropertyChk_RegularWave.RegWavePropProbe regWavePropProbe1 annotation(
          Placement(visible = true, transformation(origin = {50, 0}, extent = {{-10, -10}, {10, 10}}, rotation = 0)));
      equation
        connect(environmentBus1.rho_w, regWavePropProbe1.gdc.rho_w) annotation(
          Line(points = {{0, 0}, {40, 0}, {40, 0}, {40, 0}}, thickness = 0.5));
        connect(environmentBus1.d, regWavePropProbe1.gdc.d);
        connect(environmentBus1.Hr, regWavePropProbe1.gdc.Hr);
        connect(environmentBus1.Tr, regWavePropProbe1.gdc.Tr);
        connect(environmentBus1.omega, regWavePropProbe1.gdc.omega);
        connect(environmentBus1.T, regWavePropProbe1.gdc.T);
        connect(environmentBus1.k, regWavePropProbe1.gdc.k);
        connect(environmentBus1.epsilon, regWavePropProbe1.gdc.epsilon);
        connect(environmentBus1.zeta0i, regWavePropProbe1.gdc.zeta0i);
        connect(regular_Airy_Wave_Test1.wdoc.rho_w, environmentBus1.rho_w) annotation(
          Line(points = {{-40, 0}, {-10, 0}, {-10, 0}, {0, 0}}));
        connect(regular_Airy_Wave_Test1.wdoc.d, environmentBus1.d);
        connect(regular_Airy_Wave_Test1.wdoc.Hr, environmentBus1.Hr);
        connect(regular_Airy_Wave_Test1.wdoc.Tr, environmentBus1.Tr);
        connect(regular_Airy_Wave_Test1.wdoc.omega, environmentBus1.omega);
        connect(regular_Airy_Wave_Test1.wdoc.T, environmentBus1.T);
        connect(regular_Airy_Wave_Test1.wdoc.k, environmentBus1.k);
        connect(regular_Airy_Wave_Test1.wdoc.epsilon, environmentBus1.epsilon);
        connect(regular_Airy_Wave_Test1.wdoc.zeta0i, environmentBus1.zeta0i);
        annotation(
          experiment(StartTime = 0, StopTime = 30, Tolerance = 1e-06, Interval = 0.1));
      end Check_RegularWave;

      model Regular_Airy_Wave_Test
        extends Modelica.Blocks.Icons.Block;
        OceanEngineering.SampleSimulations.PropertyChk_RegularWave.WaveDataConnector_PropChk wdoc(n_omega_i = 1) annotation(
          Placement(transformation(extent = {{80, -20}, {120, 20}})));
        constant Real pi = Modelica.Constants.pi "Value of pi";
        constant Real g = Modelica.Constants.g_n "Acceleration due to gravity";
        parameter Modelica.SIunits.Length d = 100 "Water depth";
        parameter Modelica.SIunits.Density rho_w = 1025 "Density of Water";
        parameter Modelica.SIunits.Length Hr = 1 "Wave Height";
        parameter Real Tr = 7 "Wave Period";
        parameter Real Trmp = 20 "Wave ramp time";
        parameter Real Tdel = 0 "Wave delay";
        parameter Integer n_omega_i = 1 "Number of frequency components";
        Real T[size(omega, 1)] "Wave period";
        Modelica.SIunits.AngularFrequency omega[n_omega_i] "Frequency components of time series";
        Modelica.SIunits.Length zeta0i[n_omega_i] "Amplitude of wave component";
        Modelica.SIunits.Length SSE_X0 "Sea surface elevation at x=0 - to plot generated waves";
        Real epsilon[n_omega_i] "Random phase of wave components";
        Real k[n_omega_i](each unit = "m^-1") = OceanEngineering.Functions.waveNumberIterator(d, omega) "Wave number of component waves";
      equation
        for i in 1:n_omega_i loop
          omega[i] = 2 * pi / Tr;
          epsilon[i] = 0;
          if time < Tdel then
            zeta0i[i] = 0;
          elseif time < Tdel + Trmp then
            zeta0i[i] = sin(pi / 2 * (time - Tdel) / Trmp) * Hr / 2;
          else
            zeta0i[i] = Hr / 2;
          end if;
          T[i] = Tr;
        end for;
        SSE_X0 = sum(zeta0i .* cos(omega * time - 2 * pi * epsilon));
        wdoc.d = d;
        wdoc.rho_w = rho_w;
        wdoc.omega = omega;
        wdoc.T = T;
        wdoc.k = k;
        wdoc.epsilon = epsilon;
        wdoc.zeta0i = zeta0i;
        wdoc.SSE_X0 = SSE_X0;
        wdoc.Hr = Hr;
        wdoc.Tr = Tr;
        annotation(
          Icon(graphics = {Text(origin = {-54, 65}, extent = {{-40, 17}, {40, -17}}, textString = "H,T"), Rectangle(origin = {100, 0}, fillColor = {170, 255, 0}, fillPattern = FillPattern.Solid, extent = {{-20, -20}, {20, 20}}), Line(origin = {-15, 0}, points = {{-81, 0}, {81, 0}}), Line(origin = {-14.97, 1.01}, points = {{-81.0281, -1.00867}, {-33.0281, 42.9913}, {36.9719, -43.0087}, {80.9719, -1.00867}}, color = {0, 0, 255}, thickness = 0.5, smooth = Smooth.Bezier)}, coordinateSystem(initialScale = 0.1)),
          experiment(StartTime = 0, StopTime = 400, Tolerance = 1e-06, Interval = 0.5));
      end Regular_Airy_Wave_Test;

      connector WaveDataConnector_PropChk
        parameter Integer n_omega_i = 100 "number of wave components";
        Modelica.Blocks.Interfaces.RealOutput d;
        Modelica.Blocks.Interfaces.RealOutput rho_w;
        Modelica.Blocks.Interfaces.RealOutput omega[n_omega_i];
        Modelica.Blocks.Interfaces.RealOutput T[n_omega_i];
        Modelica.Blocks.Interfaces.RealOutput k[n_omega_i];
        Modelica.Blocks.Interfaces.RealOutput epsilon[n_omega_i];
        Modelica.Blocks.Interfaces.RealOutput zeta0i[n_omega_i];
        Modelica.Blocks.Interfaces.RealOutput SSE_X0;
        Modelica.Blocks.Interfaces.RealOutput Hr;
        Modelica.Blocks.Interfaces.RealOutput Tr;
        annotation(
          Icon(graphics = {Rectangle(fillColor = {170, 255, 0}, fillPattern = FillPattern.Solid, extent = {{-100, 100}, {100, -100}}), Text(origin = {-1, 3}, extent = {{-67, 41}, {67, -41}}, textString = "%name")}, coordinateSystem(initialScale = 0.1)));
      end WaveDataConnector_PropChk;

      connector GeneralDataConnector
        parameter Integer n_omega_i = 1 "number of wave components";
        Modelica.Blocks.Interfaces.RealInput omega[n_omega_i];
        Modelica.Blocks.Interfaces.RealInput T[n_omega_i];
        Modelica.Blocks.Interfaces.RealInput k[n_omega_i];
        Modelica.Blocks.Interfaces.RealInput epsilon[n_omega_i];
        Modelica.Blocks.Interfaces.RealInput zeta0i[n_omega_i];
        Modelica.Blocks.Interfaces.RealInput rho_w;
        Modelica.Blocks.Interfaces.RealInput d;
        Modelica.Blocks.Interfaces.RealInput Hr;
        Modelica.Blocks.Interfaces.RealInput Tr;
        annotation(
          Icon(graphics = {Rectangle(fillColor = {255, 85, 0}, fillPattern = FillPattern.Solid, extent = {{-100, 100}, {100, -100}}), Text(origin = {-1, 3}, extent = {{-67, 41}, {67, -41}}, textString = "%name")}, coordinateSystem(initialScale = 0.1)));
      end GeneralDataConnector;
    end PropertyChk_RegularWave;
  end SampleSimulations;
  annotation(
    Icon(graphics = {Polygon(origin = {0, 31}, lineColor = {0, 170, 255}, fillColor = {0, 170, 255}, fillPattern = FillPattern.Solid, lineThickness = 0, points = {{-140, -59}, {-36, -37}, {4, 9}, {54, 13}, {82, -19}, {36, 3}, {16, -53}, {150, -51}, {-140, -59}}, smooth = Smooth.Bezier), Rectangle(origin = {0, -90}, lineColor = {0, 170, 255}, fillColor = {0, 170, 255}, fillPattern = FillPattern.Solid, lineThickness = 0, extent = {{-100, 70}, {100, -10}})}));
end OceanEngineering;

The OceanEngineering library in this download contains the component-models to simulate waves and depth varying current. The 'SampleSimulations' section contains the models based on which results are presented in the paper titled Modelica Component Models for Oceanic Surface Waves and Depth Varying Current, submitted to the American Modelica Conference 2020.

The regular wave component model is slightly modified to generate and transmit extra information for plotting the results presented in the paper and is named Regular_Airy_Wave_Test. A modified wave data connector is also specified to transmit this information to the environment bus. Another component model named RegWavePropProbe is specified to recieve the generated variables via a GeneralDataConnector and perform the required calculations to plot the various properties such as velocities, pressures, profiles etc.

The outputs from the irregular wave component model and current component model can be plotted directly from the respective component models.

Detailed instructions are given below:

********************************************************************************************************************************************************
**Instructions to plot graphical results discussed in the  above paper **
********************************************************************************************************************************************************
a. Load the OceanEngineering library in OMEdit.

****Figure 8-a
a. Open the component model Check_RegularWave under SampleSimulations>PropertyChk_RegularWave.
b. Set parameters Hr=1 m,Tr=3 s, d=10m, Trmp=10s and Tdel=5 s in the regular_Airy_Wave_Test1 wave component model and run simulation for 30s with time step of 0.1s.
c. In the variables browser, open the regular_Airy_Wave_Test1 results and plot SSE_X0.

****Figure 8-b
a. For the same simulation above, open the regWavePropProbe1 results and plot x Vs eta in an array parameteric plot. 
b. Adjust slider for simulation time of 15,25.9,16.5 and 17.1 s to obtain the respective plots.

****Figure 9-a
a. Re-run the simulation setting Trmp=Tdel=0s.
b. In the variables browser, open regWavePropProbe1 results and plot x_lng_lng Vs eta_lng_lng in an array parameteric plot. 
c. Change parameter Tr under the regular_Airy_Wave_Test1 model to plot the profiles for waves of different periods.

****Figure 9-b
a. In the variables browser, open the regWavePropProbe1 results and plot x_del Vs eta_del in an array parameteric plot to plot the wave profiles for different time periods.
b. In the variables browser, open the regWavePropProbe1 results and in a parametric plot window plot X[1] Vs Z[1]......X[11] Vs Z[11] to plot the water particle trajectories.

****Figure 10-a
a. For the same simulation above, set parameter Tr=6s in regular_Airy_Wave_Test1 and re-simulate.
b. In the variables browser, open regWavePropProbe1 results and plot x_lng Vs eta_lng in an array parameteric plot.

****Figure 10-b to f
a. For the same results above, under regWavePropProbe1 results, u0...u4 and w0...w4 gives the horizontal and vertical water particle velocities for water particles with mean vertical coordinates given under z. These data are used to plot the quiver plots given in 10-b to 10-f.

****Figure 11 a to c
a. For the same results above, under regWavePropProbe1 results, plot p_hydx0 Vs z_pr in an array parametric plot to plot the hydrostatic pressure at x0.
b. plot totp_x0 Vs z_pr to plot the total pressure at x0.
c. repeat for p_hydx1, p_hydx2, totp_x1, totp_x2 Vs z_pr for plotting respective pressures at x1, x2 locations.

**** Figure 12 a to c
a. Open Components>Waves>IrregularWave>IRW_PM_RDFCWI set parameters Hs=1 m, d=100s, omega_min=0.03141 rad/s and omega_max=3.141 rad/s, Tdel=0s, and Trmp=10s, and simulate for 400 s with time step of 0.2s.
b. Plot omega Vs S in an array parametric plot window to get the spectrum plot.
c. Plot SSE_X0 to get the irregular sea surface elevation at x=0.

*** Figure 13
a. Open Components>CurrentProfiles>CurrentProfile_4pt and simulate for 10s with time step of 0.5s.
b. In the results, plot Ucg Vs Zcg in an array parametric plot for t=0,0.8,2.5,5.0 s to get the required plots.
*****************************************************
for clarification, send mail to savinvis@gmail.com
*****************************************************

 
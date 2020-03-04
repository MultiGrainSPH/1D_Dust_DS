DustyWaveMulti.sce is an accompanying material for the paper 'Fast method to simulate dynamics of two_phase medium with intense interphase interaction. I. Gas-dust mixture with polydisperse particles, linear drag, one-dimensional tests' 
DustyWaveMulti.sce is a generator of a reference solution to the problem of  sound wave propagation in the medium of isothermic gas-polydisperse dust. 

To obtain a reference solution one needs: 

1. Download, install (scilab.org) and launch Scilab 
2. Open file 'DustyWaveMulti.sce' in the Scilab SciNotes file editor. 
3. Set the initial data In the file 'DustyWaveMulti.sce'
4. Run the file from Scilab, as a result of its execution should appear a graphic window and output files. 
Depending on the value of the variable 'type_Of_Plot', the graphics window shows: (x) the density and velocity of gas and dust fractions at t = 0 and t = 't_Selected' (t) normalized velocities and densities at the point x = 'x_Selected' (a description of the normalized values is given in files of the gas / dast type (x = 'x_Selected')) 

Constrains on initial data are rather natural:
- numerical values of the initial data, including elements of the arrays 'rho_Arr' and 'ts_Arr', must be positive (the values of 't_Selected' and 'x_Selected' may be zero) 
- value 'type_Of_Plot' can only be equal to 'x' or 't' -for rewriting output files containing data for given values of 'x' or 't' (in their name there is a corresponding variable), you need to set the value 'x' or 't' for the value 'type_Of_Plot', respectively, and the value 'on' for 'overwrite_Files', to disable overwriting, you can set any other value instead of 'on' (for example, 'off') 
- the output file with the name specified in 'name' will be overwritten with the value 'overwrite_Files' = 'on' regardless of the value of 'type_Of_Plot' -when specifying the file name 'name', you don’t need to specify the file type, it is the default set as '* .txt' -when specifying the full path 'path', it is not necessary to write the character '\' at the end, it will be written automatically (e.g. path = 'D: \ articles \ 1DMultigrain \ DustyWaveMultigrain \ ver20191120';) 
- in the array of colors 'color_Arr_Article' the elements must be no less than N + 1, where N is the number of elements in the array 'rho_Arr' 

The contents of the output files: 
-Files with a name of type 'gas (...)', 'dust (...)' contain values of densities and velocities for given values of the variables 'x' or 't' specified in the file names.
-In the file 'InitialData' contains the values of the constant Cs, initial densities, times relaxation of dust fractions, as well as deltas of densities and velocities of gas and dust fractions at t = 0. 
For instance, In the 'InitialData' we have
 
Initial densities 
 g: 1.000000 
 1: 0.100000 
 2: 0.233333 
 3: 0.366667 
 4: 0.500000 

Density deltas at t = 0 
 g: 0.0001*(1.000000*cos(6.283185*x) - 0.000000*sin(6.283185*x)) 
 1: 0.0001*(0.080588*cos(6.283185*x) + 0.048719*sin(6.283185*x)) 
 2: 0.0001*(0.091607*cos(6.283185*x) + 0.134955*sin(6.283185*x)) 
 3: 0.0001*(0.030927*cos(6.283185*x) + 0.136799*sin(6.283185*x)) 
 4: 0.0001*(0.001451*cos(6.283185*x) + 0.090989*sin(6.283185*x)) 

To use the reference solution, you need to set the initial density values as the sum of their initial values and density deltas, initial velocity values as velocity deltas  
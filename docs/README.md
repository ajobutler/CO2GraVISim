<!-- README.md -->

# CO2GraVISim

This is an overarching readme for the **CO<sub>2</sub>** **Gra**vity-controlled, **V**ertically **I**ntegrated **Sim**ulator, a reduced-physics model of the flow of CO<sub>2</sub> in a confined porous reservoir.
The model incorporates spatially variability in ceiling and basement topography, porosity, and permeability, along with residual trapping and dissolution.

This code and model were developed as part of the StrataTrapper project. All comments and queries can be directed towards [ajb278@cam.ac.uk](mailto:ajb278@cam.ac.uk).

---

## Index

1. [Full Setup](#full-setup)
2. [Model Description](#model-description)
3. [File Structure](#file-structure)
4. [Input Formats](#input-formats)
5. [Code Summary](#code-summary)
6. [Code Output](#code-output)
7. [Licence](#licence)

## Full Setup

To fully run the simulator from scratch:

0. Compile the executable files - using `build_CO2GraVISim__Windows.bat` or `build_CO2GraVISim__Linux.sh`
1. Set inputs:
    - grid parameters in `grid_parameters.txt`
    - spatial arrays - `ceil_topo.txt`, `base_topo.txt`, `Permeability.txt`, and `Porosity.txt`,  - e.g. via `CO2GraVISim_InputGen`
    - plot times - e.g. using `Generate_plot_times.py`
    - flow parameters - e.g. using values produced by `Nondimensional_parameters.py`
    - injection locations and injection profile in `Injection_locations.txt` and `Injection_profile.txt` (and check using `Reservoir_preview.py`)
    - boundary conditions in `Boundary_conditions.txt`
2. Run `CO2GraVISim_single_run`
3. Process output - using `Overview_plot.py` and/or `Slices_plot.py`

Once the apropriate executables have been created, Step 0 can be skipped and the input files in Step 1 can be amended as necessary.

By default, the build files attempt to use the `gfortran` compiler to construct the executable files. The files involved in Step 0 can be modified simply to instead use a different compiler, such as `ifort`, if necessary.

On Linux systems, it may be necessary to modify the read-write status of the build file. This can be done by running ```chmod 744 build_CO2GraVISim__Linux.sh``` from within the relevant file path.

---

## Model Description

Assuming that:

1. The flow of the CO<sub>2</sub> is predominantly horizontal,
2. Gravity drives the CO<sub>2</sub> to build up beneath the caprock a short distance after the injection point,

we can take the pore pressure to be hydrostatic in the vertical direction at leading order, and vertically integrate the local continuity equations for the two fluids to arrive at a trio of coupled equations that control

- $h(x,y,t)$ - the thickness of the mobile CO<sub>2</sub> region,
- $h_{res}(x,y,t)$ - the thickness of the residually trapped CO<sub>2</sub> region,
- $P_{a}$ - the deviation from hydrostatic pressure in the ambient.

A full description of the governing equations solved here can be found in `\docs\CO2GraVISim_model_description.pdf`.

### Nondimensionalisation

This code solves the nondimensional version of the mathematical model described above. Vertical length scales are scaled by the representative thickness of the storage interval, $\mathcal{D}$, and the injection flux is scaled by a representative value $\mathcal{Q}$. This then produces corresponding time and pressure scales

$$\mathcal{T} = \frac{\phi_{0} \mathcal{D}^{3} }{\mathcal{Q}}
\, , \quad\quad
\mathcal{P} = \frac{\mu_{c}}{k_{0}}\frac{Q}{\mathcal{D}} \, ,
$$

where $\phi_{0}$ is a representative value of the porosity, $k_{0}$ that of the permeability, and $\mu_{c}$ is the viscosity of the injected CO<sub>2</sub>. This introduces the parameter

$$
\Gamma = \frac{u_{b}}{u_{Q}} =
\frac{ (\Delta \rho g k_{0} / \mu_{c})}{ (\mathcal{Q} / \mathcal{D}^{2} )} \, , 
$$

which is the ratio of the Darcy velocities set by buoyancy and injection respectively, and where $\Delta \rho = \rho_{a} - \rho_{c}$ is the density difference between the ambient and CO<sub>2</sub>.

The python script `Nondimensional_parameters.py` can be used to calculate these dimensional scales, as well as the parameter values used in `Flow_parameters.txt`, from given dimensional values. 

---

## File Structure

```bash
├───build_CO2GraVISim__Linux.sh
├───build_CO2GraVISim__Windows.bat
├───CO2GraVISim_install.zip
│
├───docs
│       README.md
│
├───Input
│       base_topo.txt
│       boundary_conditions.txt
│       ceil_topo.txt
│       flow_parameters.txt
│       grid_parameters.txt
│       injection_locations.txt
│       injection_profile.txt
│       permeability.txt
│       porosity.txt
│       target_plot_times.txt
│
├───Output
│   ├───Current_Pressure
│   │       placeholder_Current_Pressure.txt
│   │
│   ├───Current_Thickness
│   │       placeholder_Current_Thickness.txt
│   │
│   └───Other
│           placeholder_Other.txt
│
├───plots
│   ├───Overview_plot
│   │   └───temp
│   │           placeholder_Overview_plot.txt
│   │
│   ├───Reservoir_preview
│   │       placeholder_Reservoir_preview.txt
│   │
│   └───Slices_plot
│       └───temp
│               placeholder_Slices_plot.txt
│
├───python_scripts
│       Generate_Plot_times.py
│       Nondimensional_parameters.py
│       Overview_plot.py
│       Reservoir_preview.py
│       Slices_plot.py
│
└───src_files
        CO2GraVISim_global.f90
        CO2GraVISim_injection_profiles.f90
        CO2GraVISim_InputGen.f90
        CO2GraVISim_input_parameters.f90
        CO2GraVISim_single_run.f90
        CO2GraVISim_solver.f90
```

---


## Input Formats

### Spatial arrays

The spatial arrays describing the ceiling and basal topography of the reservoir, along with the horizontally varible spatial structure of the permeability and porosity, and read in from `ceil_topo.txt`, `base_topo.txt`, `permeability.txt`, and `porosity.txt`, respectively. These $(n_x \times n_y)$ arrays are stored as $n_y$ rows of $n_x$ real-valued numbers, separated by spaces.

### Grid parameters

Following a comment line specifying the formatting, `grid_parameters.txt` specifies on consecutive lines:

- $n_{x}$ - the number of grid points in $x$,
- $n_{y}$ - the number of grid points in $y$,
- $\Delta x$ - the nondimensional grid spacing in $x$,
- $\Delta y$ - the nondimensional grid spacing in $y$.

The grid spacing values are nondimensionalised by the representative reservoir thickness $\mathcal{D}$, in this model.

### Flow parameters

These values of the relevant flow parameters are set on consecutive lines in `flow_parameters.txt`, following the comment line describing the formatting:

- $M$ - the viscosity ratio,
- $\Gamma = u_b / u_Q$ - the buoyancy-flux ratio,
- $s_{cr}$ - the residual saturation of CO<sub>2</sub>,
- $s_{ai}$ - the irreducible saturation of the ambient,
- $C_{sat}$ - the total amount of CO<sub>2</sub> that can dissolve in the ambient by volume,
- $q_{d}$ - the nondimensional convective flux rate

### Boundary conditions

Dirichlet or Neumann boundary conditions can be specified for $h$ and $P_{a}$ in `boundary_conditions.txt`. Following the comment lines there, these are done on one line for $h$, and the next for $P_{a}$. Each consists of a list of $4$ space-separated numbers, corresponding to the North, East, South, and West boundaries respectively. The value $1$ corresponds to the Dirichley boundary condition $f=0$, while $2$ corresponds to the Neumann boundary condition $df/dn=0$. E.g.

``` bash
2 1 1 1

1 2 2 1
```

would correspond to a Neumann boundary condition for $h$ on the North side of the domain and Dirichlet boundary conditions on the rest, while the pressure $P_{a}$ would have Dirichlet boundary conditions on the North and West sides and Neumann on the East and South sides.

### Injection Locations

Multiple injection locations can be specified via `injection_locations.txt`.
The first line after the comments sets the total number of injection locations, $n_{inj\,locs}$. The subsequent lines then consist of the $x$ and $y$ indices for each location in turn, e.g.

```bash
3
150 75
120 50
180 75
```

specifies $3$ wells, located at the indices $(150,75)$, $(120,50)$, and $(180,75)$, respectively.

### Injection profiles

Injection can be specified in terms of piecwise-constant profiles corresponding to multiple injection locations, via `injection_profiles.txt`. The first line after the comments corresponds to the number of times at which the flux profile changes value, $n_{flux\, vals}$.

There then follows that number of lines. The first value sets the time at which the flux changes, and the subsequent values the corresponding flux value at each of the injection locations specified in `injection_locations.txt` (see [Injection Locations](#injection-locations)).

E.g.

``` bash
6
0 1 0   0
1 0 0   0
2 0 0.5 0
3 0 0   0
6 0 0.5 0.5
7 0 0   0
```

sets $6$ times at which the flux changes. At $t=0$ the first well turns on with flux $1$, while the others remain dormant. At $t=1$ this first well turns off. At $t=2$ the second well turns on with flux $0.5$, while the others remain dormant. At $t=3$ the second well turns off. At $t=6$ both the second and third well turn on, with flux $0.5$, and subsequently turn off at $t=7$.

### Target Plot Times

After the comment line, `target_plot_times.txt` sets first the number, and then the values, of the nondimensional plot times, on consecutive lines. The solver will then aim to record profiles for $h$, $P_{a}$, etc., at these times, though due to the adaptive timestepping routine may output at slightly different values of $t$.

---

## Code Summary

The solver consists of a a set of modules and subroutines written in Fortran90, with accompanying Python scripts for input generation and output interpretation. The details of each of these are listed below.

### CO2GraVISim_global.f90

This module contains fixed parameters that are used by the solver, and not likely to be changed between individual runs. These include

- the numerical working precision, `wp`
- the maximum number of iteration, `itermax`
- the maximum run time, `tmax`
- the initial step size, `dt_init`
- the minimum step size, `dt_min`
- the maximum relative error for the time-stepping routine, `errmax_t_step`
- the maximum relative error for the iterative Pressure solver, `errmax_P_iter`
- the relaxation parameter for the iterative Pressure solver, `P_iter_omega`
- the start time, `t0`

If any of these values are altered, then the codes need to be recompiled [(see Step 0 in the Full Setup)](#full-setup).

### CO2GraVISim_injection_profiles

This module reads in the injection locations and injection profile specified in `/Input/injection_locations.txt` and `/Input/injection_profile.txt`. This is then used by `CO2GraVISim_solver` to calculate the appropriate injection flux array at each timestep.

### CO2GraVISim_InputGen

This program produces demonstration input arrays (ceiling and basement topography, porosity, and permeability) that can then be used for a run of the solver.

The user is first asked to choose from a list of preset porosity distributions:

``` bash
 Choose a Porosity distribution:
 1)Uniform 2)Low-porosity rectangle 3)High-porosity channel 4)Checkerboard
```

1. **Uniform**: This produces a uniform porosity distribution. No other parameters need be specified.
2. **Low-porosity rectangle**: A rectangular region of lower porosity, covering the (nondimensional) region $-3\leq x \leq 3, 2 \leq y \leq 5 $. The user also specifies the ratio of the minimum porosity (inside the rectangle) to the maximum porosity (outside the rectangle).
3. **High-porosity channel**: A channel of higher porosity, covering the (nondimensional) region $2 \leq x \leq 5$. The user also specifies the ratio of the minimum porosity (outside the channel) to the maximum porosity (inside the rectangle).
4. **Checkerboard**: A checkerboard pattern (or series of vertical or horizontal rows) of alternating squares of high and low porosity. The top left square has the high porosity value. The user specifies the ratio of the minimum porosity to the maximum porosity, as well as the number of grid points in $x$ and $y$ that each checkboard tile takes, e.g.

``` xml
You chose Checkerboard (4)
Choose the no. of points in x and y of a checkerboard tile:
grid dimensions are (Nx,Ny) =          200         100
12
12
Checkerboard values are           12          12
Choose a value for (phi_min / phi_max):
0.75
```

The **permeability** distribution is then generated from the chosen porosity via the Kozeny-Carman relation

$$ k = \frac{ (\phi_{0}\phi)^{3} }{ ( 1 - (\phi_{0}\phi) )^{2} } $$

where a nominal value for the scale $\phi_{0}$ of $0.15$ is chosen. This permeability is then rescaled so that mean value is 1:

$$ \bar{k} = \frac{1}{n_{x}*n_{y}}\sum_{i,j} k(x_{i},y_{j}) = 1 \, .  $$

Finally, the reservoir topography is generated. The ceiling topography takes the form of a sloped plane, with the user choosing the (nondimensional) slopes in the $x$ and $y$ directions - positive values of the two slopes result in a shallower ceiling in the 'West' and 'South' directions, respectively.
The basement topography then mirrors the ceiling topography, one non-dimensional unit deeper.

## Code Output

The outputs produced by a single run of `CO2GraVISim` are stored in `\Output\`.

`\Output\Current_Pressure\` contains the saved profiles of the ambient pressure $P_{a}$ at each of the plot times as `P00.txt`, `P01.txt`,...

Similarly, `\Output\Current_thickness\` contains the saved profiles for the mobile current thickness $h$ and the residually trapped current thickness $h_{res}$ as `h00.txt`, `h01.txt`,... , and `h_res00.txt`, `h_res01.txt`,...

Finally, `\Output\Other\` contains the remaining saved data produced during the run. These include a copy of the injection locations used (`injection_locations.txt`); the time values corresponding to the output data (`plot_times.txt`); a list of the various parameter values used in the calculation (`parameter_values.txt`); and various Volumes recorded at each calculation step (`Volumes.txt`).

The recorded parameter values in `parameter_values.txt` are $n_{x}$, $n_{y}$, $\Delta x$, $\Delta y$, $M$, $\Gamma$, $s_{cr}$, $s_{ai}$, $C_{sat}$, and $q_{d}$ [(see Input Formats)](#input-formats).

The data columns stored in `Volumes.txt` correspond to the time of that calculation step; the total volume of the mobile CO<sub>2</sub> phase; the total volume of the trapped CO<sub>2</sub> phase; and the total volume of CO<sub>2</sub> that has been injected up to that point.

---

## Licence

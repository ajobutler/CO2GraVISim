<!-- README.md -->

# CO2GraVISim

This is an overarching readme for the **CO<sub>2</sub>** **Gra**vity-controlled, **V**ertically **I**ntegrated **Sim**ulator, a reduced-physics model of the flow of CO<sub>2</sub> in a confined porous reservoir.
The model incorporates spatially variable ceiling and basement topography, porosity, and permeability, as well as residual trapping and dissolution.

---

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

After the first run, Step 0 can be skipped and the input files in Step 1 can be amended as necessary.

---

## Model Description

Assuming that:

1. The flow of the CO<sub>2</sub> is predominantly horizontal,
2. Gravity drives the CO<sub>2</sub> to build up beneath the caprock a short distance after from the injection point,

we can take the pore pressure to be hydrostatic in the vertical direction at leading order, and vertically integrate the local continuity equations for the two fluids to arrive at a trio of equations that control

- the thickness of the mobile CO<sub>2</sub> region, $h(x,y,t)$,
- the thickness of the residually trapped CO<sub>2</sub> region, $h_{res}(x,y,t)$,
- the pressure in the ambient, $P_{a}$.

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

If any of these values are altered, then the codes need to be recompiled via `build_CO2GraVISim__Windows.bat` or `build_CO2GraVISim__Linux.sh`.

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

```bash
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

## Licence

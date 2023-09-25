<!-- README.md -->

# CO2GraVISim

This is an overarching readme for the **CO<sub>2</sub>** **Gra**vity-controlled, **V**ertically **I**ntegrated **Sim**ulator, which is a reduced-physics model the flow of CO<sub>2</sub> in a confined porous reservoir.

## Summary

This code implements a reduced physics model of CO<sub>2</sub> injected into a porous reservoir. The model incorporates spatially variable ceiling and basement topography, porosity, and permeability, as well as residual trapping and dissolution.

---

## Full Setup

- Compile the executable files - using `build_CO2GraVISim__Windows.bat` or `build_CO2GraVISim__Linux.sh`
- Set inputs:
  - Set grid parameters in `grid_parameters.txt`
  - Set input arrays - `Permeability.txt`, `Porosity.txt`, `ceil_topo.txt`, and `base_topo.txt` - e.g. via `CO2GraVISim_InputGen`
  - Set plot times - e.g. using `Generate_plot_times.py`
  - Set flow parameters - e.g. using values produced by `Nondimensional_parameters.py`
  - Set injection locations and injection profile in `Injection_locations.txt` and `Injection_profile.txt` (and check using `Reservoir_preview.py`)
  - Set boundary conditions in `Boundary_conditions.txt`
- Run `CO2GraVISim_single_run`
- Process output - using `Overview_plot.py` and/or `Slices_plot.py`

---

## Code Summary

Solver codes written in Fortran90, with accompanying Python scripts for input generation and output interpretation.

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

## Licence

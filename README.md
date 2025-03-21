<!-- README.md (19/03/25) -->

# CO2GraVISim

**CO2GraVISim** (**CO2** **Gra**vity-controlled, **V**ertically **I**ntegrated **Sim**ulator) is a reduced-physics simulator for modelling the flow of CO2 injected into sub-surface porous, saturated, reservoirs. Based on the assumption of

- Predominantly horizontal flow
- vertical migration to the caprock being quick relative to the total migration time

the governing equations are vertically integrated across the relevant regions of the reservoir. This reduces the complexity of the governing equations, allowing the behaviour to be simulated much more rapidly in comparison to traditional simulators that solve for the fully three-dimensional flow.

The codebase provided here consists of a number of Fortran files comprising the main simulator, as well as several Python scripts for analysing the inputs and outputs of the simulations.

CO2GraVISim has been developed as part of the **StrataTrapper** project. **Version 1** of CO2GraVISim was developed for an interim milestone of the project, and works with fully non-dimensional input and output data. This has since been further developed into **Version 2**, which now has dimensional inputs and outputs, a more user-friendly input format, extended functionality and added analytical features. As part of this, CO2GraVISim uses FoX, an XML library for Fortran (see [https://fortranwiki.org/fortran/show/FoX](https://fortranwiki.org/fortran/show/FoX)). FoX is licenced under a BSD license.

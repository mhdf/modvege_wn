# Model notes

## Inherited formulation

Daily potential growth is incoming solar radiation × 0.47 (PAR fraction) ×
3 g dry matter/MJ PAR (maximum radiation-use efficiency) × intercepted fraction
`1 - exp(-0.6 × LAI)` × 10 (g/m² to kg/ha). The result is multiplied by a
seasonal term and environmental stress terms. Growth is partitioned between
green vegetative and green reproductive pools. Senescence transfers part of
green biomass to the corresponding dead pool; abscission removes dead biomass.
The simple water balance adds rainfall and subtracts estimated crop ET from
available root-zone water. See the original `legacy/ModVege.py` for the exact
historical parameterization.

## Intentional changes from the legacy file

1. A weather day is simulated once. The old inclusive slices repeated cut dates
   and overwrote their first result. A cut now occurs after that day's growth;
   daily output on the cut date represents pre-cut standing biomass.
2. Cutting yield is returned explicitly as biomass removed above a residual
   height for each of the four pools. The legacy experiment scripts calculated
   yield later from daily state and residual assumptions.
3. On cold days, average pool ages remain defined. The old implementation
   could use unassigned age-change variables when mean temperature was ≤ 0 °C.
4. Senescence and abscission cannot remove more than the biomass present in a
   pool. Digestibility is bounded below by its group-specific minimum.
5. Invalid weather, group fractions, soil properties and cut dates fail with
   actionable errors. The original code often printed a message or continued.
6. The initial state matches the legacy assumptions: full available soil water,
   5 cm residual green and dead vegetative biomass, and initial thermal ages of
   400 and 500 degree-days, respectively. Other initial states are future work.

These changes mean outputs are **not assumed numerically identical** to the
old script. The implementation should be treated as a new, documented version.

Neither the empirical water floor nor NI should be interpreted as a measured
groundwater supply or a mechanistic nitrogen balance.

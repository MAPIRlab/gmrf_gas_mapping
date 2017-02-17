# gmrf_gas_mapping
Implementation of the Time-Variant gas distribution mapping using Gaussian Markov Random Fields.

This algorithm contributes an effective solution to two important problems that have been disregarded so far in the mapping of gases: First, obstacles in the environment (walls, furniture, ...) do affect the gas spatial distribution. Second, when combining odor measurements taken at different instants of time, their 'ages' must be taken into account to model the ephemeral nature of gas distributions. In order to incorporate these two characteristics into the mapping process we propose modeling the spatial distribution of gases as a Gaussian Markov random field (GMRF). This mathematical framework allows us to consider both: (i) the vanishing information of gas readings by means of a time-increasing uncertainty in sensor measurements, and (ii) the influence of objects in the environment by means of correlations among the different areas.

For more info refer to: http://mapir.isa.uma.es/mapirwebsite/index.php/mapir-downloads/papers/193

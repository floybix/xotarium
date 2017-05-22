
* s/pressure/depth/
* s/air/water/


* semi-persistence mechanism - a ratchet? - for e.g. responding to touch.
  - epigenetic marks?
  - could affinity/influence > 1 be enough?


* perf
  - step fn in java/C?


* viz inputs to GRN cells / over time

* test response of GRNs to various inputs
  - (vs a defined reference state)

* separate enviro seed from evo seed


* return cppn mutation info




* movement model
  - do we need phase-off?
  - do we need oscillator? (vs internal cycles)
  - should GRN output adjust frequency (oscillator-freq)


* larger / composite caves

* avoid global state

* time as an input to cppn (instead of phase/freq as output)

* use larger particles for air/water but not for creatures
  - hack LiquidFun to allow non-uniform radius



# Evolution

* behavioural distance?
  - (for complex behaviours, many/all instances will be novel in some respect?)
  - OTOH: if we can group behaviours into qualitatively equivalent sets then ok.

* can run LF worlds in parallel?

* speciation

* crossover at phylogenetically-corresponding points
  - history markers?
  - topo analysis?

* record population measures over time
  - genome sizes
  - genetic diversity
  - behavioural diversity
  - genetic novelty
  - behavioural novelty

* attach mutation parameters to individuals, and evolve them
  - initially, large mutations should be favoured; later, smaller

* transit serialisation



## Behaviour

* set of plant cells eaten

* energy balance / survival time

* coevolution of carnivores and herbivores
  - consume cell by cell?
  - delay cell death to after sensing & messaging


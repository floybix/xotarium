# xotarium

**This is a half-finished idea. It barely works.**

Evolving behaviours of soft creatures in a particle world.

Notable:
* GRNs
* signals between cells in a creature body
* novelty search for interesting movement behaviours
* 2D particle simulation of creature and its (aquatic) environment

## Usage

Kick off a novelty search for 15 generations, using 42 as the random seed.
```
lein run -m org.nfrac.xotarium.evo/run 15 42
```

That last command will have created files named
`beh-archive-seed42-gen14.edn` etc. Each is a snapshot of the
creatures found so far with the most _interesting_ behaviours; for now,
_interesting_ means movement ending up in a different location to any
previously seen creature, after 5 seconds.

You can see the behaviours stored in a file using Replay:
```
lein run -m org.nfrac.xotarium.replay/run beh-archive-seed42-gen14.edn
```

That will open a graphical display window replaying the behaviour of
the first interesting creature. The creature's body is a composite of
red cells and white cells. Red is flexible muscle and white is rigid
bone. Simulation time is shown in the bottom left corner; note the
creature's location at t = 5 seconds.

Press `n` to go on to the next creature found. The first ones
consitute the entire current generation: at the top of the window,
"gen" will be blank. Following that is the archive of interesting
behaviours and these will have a generation number shown.

Press `s` to show inter-cell signals as colours. There are 3
independent signals and these are used as red/green/blue colour
components.

Press `` to open a visualisation of the Gene Regulatory Network in
action.

To continue a run from a snapshot file:

```
lein run -m org.nfrac.xotarium.evo/run-more beh-archive-seed42-gen14.edn 10 42
```



## License

Copyright © 2017-2018 Felix Andrews

Distributed under the Eclipse Public License either version 1.0 or (at
your option) any later version.

#!/bin/bash
# Paired generation: continuum and fine-lattice control, 16 seeds each, kappa 0.618 and 0.757.
set -e
B=$HOME/Proyectos/hyperuniform2D/build/hyper2D
EX=$HOME/Proyectos/hyperuniform2D/examples
for kap in 0.618 0.757; do
  d=k${kap/./}
  mkdir -p $d/cont
  for i in 0 1 2 3; do mkdir -p $d/lat$i; done
  sed -e "s/^kappa: .*/kappa: $kap/" -e 's/^realizations: 1/realizations: 16/' -e 's/^num_threads: 0/num_threads: 6/' $EX/continuum.in > $d/cont/c.in
  for i in 0 1 2 3; do
    sed -e "s/^kappa: .*/kappa: $kap/" -e "s/^seed: .*/seed: $((1729+4*i))\narr_index: $((4*i))/" \
        -e 's/^realizations: 1/realizations: 4/' -e 's/^num_threads: 0/num_threads: 3/' $EX/lattice_control.in > $d/lat$i/l.in
  done
done
for kap in 0618 0757; do
  (cd k$kap/cont && $B c.in > log.txt 2>&1)
done
pids=()
for kap in 0618 0757; do
  for i in 0 1 2 3; do (cd k$kap/lat$i && $B l.in > log.txt 2>&1) & pids+=($!); done
  wait "${pids[@]}"; pids=()
done
echo done

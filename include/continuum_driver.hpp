#pragma once

#include "parametros.hpp"

/// Runs one realization of mode continuum: random start (random sequential addition if there is a hard
/// core), L-BFGS on the weighted spectral residual plus the hard-core term, and the outputs (minimization
/// log, S(k), obstacle file and its .sk companion).
void run_continuum(const Parameters& p, unsigned seed, int arr);

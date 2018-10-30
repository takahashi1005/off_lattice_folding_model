#pragma once

static const int ProteinLength = 46;
static const double energy_per_atom=1.0;
static const double k=1.0e+4;
static const double D=10.0;
static const int dim = 3;

// List of Neighborhood
static const double range_neighborhood=10.0;
static const int term_increment_neighborhood=1000;

// Soft Core Potential
static const double epsilon_h = 1.0;
static const double epsilon_p = 1.0;
static const double sigma = 0.3;

// ATP binding
static const double k_bound = 2.0e-1;
static const double k_unbound = 3.0e-1;
static const int ATPattachmentPoint = 30;

// time
static const double dt = 1.0e-5;
static const double relaxation_time = 3000.0;
static const double tmax=100.0;
static const int term=1000;

#pragma once
#include <vector>
#include "parameter.hpp"
#include "vector_calc.hpp"
//%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
// Define class
//%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

enum color{
  H, P, N
};

enum ligand_state{
  bound, unbound
};

static const color Sequence[] = {
  H, H, H, H, H, H, H, H, H,
  N, N, N,
  P, H, P, H, P, H, P, H,
  N, N, N,
  H, H, H, H, H, H, H, H, H, 
  N, N, N,
  P, H, P, H, P, H, P, H, P, H, 
  P
};

static const std::vector<std::vector<int>> ATP_interaction_pair = {
  {},{},{},{},{},{},{},{},{},{},
  {},{},{},{},{},{},{},{},{},{},
  {},{},{},{},{},{},{},{},{},{},
  {},{},{},{},{},{},{},{},{},{},
  {},{},{},{},{},{}
};

static const std::vector<double> tweeze_position = {0.0,0.0,5.0};

struct monomer {
  std::vector<double> q;
  std::vector<double> p;
  color c;
  std::vector<int> neighbor;
  monomer(){
    q.assign(dim,0.0);
    p.assign(dim,0.0);
  }
};


class MyClass{
  std::vector<monomer> Protein;
  ligand_state ATP;
// Equation of motion
		double bond_interaction(int n, int i);
		double HPN_interaction(int n, int i);
    double ATPinteraction();
		double F(double t, int n, int i);
		double G(double t, int n, int i);
    // angle interaction
    double Theta(std::vector<double> v1, std::vector<double> v2, int i);
    double angle_interaction(int n, int i);
public:
		double ExternalForce;
    double x0;
	 	MyClass();
    void IncrementNeighbor();
    void Switching_ATP_binding(){
      if(ATP == unbound && MT::GetDouble() < k_bound){
        ATP = bound;
      }
      if(ATP == bound && MT::GetDouble() < k_unbound){
        ATP = unbound;
      }
    }
// Output
		double xn() {return Protein[ProteinLength-1].q[0];}
		double temperature();
    double radius_of_gyration();
    double test();
    void output(int n,std::ostream &os) {
      os << Protein[n].q << " " << Protein[n].p << std::endl;
    }
    void output_position(std::ostream &os) {
      for(int n=0;n<ProteinLength;++n){
        os << Protein[n].q << std::endl;
      }
    }
// Method of solving differential equation
    void Euler(double t);
    void LeapFrog_q(double t);
    void LeapFrog_p(double t);
    void RungeKutta(double t);

};

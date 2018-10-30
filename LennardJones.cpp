#include <chrono>
#include <cmath>
#include <fstream>
#include <iostream>
#include <random>
#include <vector>
#include "mt.h"
#include "parameter.hpp"
#include "MyClass.hpp"

//%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
// Main
//%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

int main() {
// clock
	std::chrono::system_clock::time_point start, mid, end;
	start = std::chrono::system_clock::now();
// random number
	std::random_device rnd;
	MT::SetSeed(rnd());

	double t=0.0;
  MyClass gamma;
  double progress;

	std::cout.setf(std::ios::fixed);
	std::cout.precision(3);

	std::ofstream fout("data.txt");  
  fout.setf(std::ios::fixed);
	fout.precision(6);

//  gamma.IncrementNeighbor();
/*
  fout << "#Setting" << std::endl;
  fout << "#ProteinLength = " << ProteinLength << std::endl; 
  fout << "#spring constant = " << k << std::endl;
  fout << "#diffusion coefficient = " << D << std::endl;
  fout << "#range neighborhood = " << range_neighborhood << std::endl;
  fout << "#term increment neighborhood = " << term_increment_neighborhood << std::endl;
  fout << "#epsilon_h = " << epsilon_h << std::endl;
  fout << "#epsilon_p = " << epsilon_p << std::endl;
  fout << "#sigma(HPN interactioln) = " << sigma << std::endl;
  fout << "#ATP attachment point = " << ATPattachmentPoint << std::endl;
  fout << "#time step = " << dt << std::endl;
  fout << "#relaxation time = " << relaxation_time << std::endl;
  fout << "#tmax = " << tmax << std::endl;
  fout << "#term = " << term << std::endl;

  fout << "#time Temp radius_of_gyration" << std::endl;
  fout << t << " " << gamma.temperature() << " " << gamma.radius_of_gyration() << std::endl;
*/
//----------------------------------------------------------------------
// solve using Runge Kutta
//----------------------------------------------------------------------
  
  int i=0;
	while(t<=tmax){
    gamma.Euler(t);
    t+=dt;
    ++i;
/*
    if(++i % term_increment_neighborhood){
     gamma.IncrementNeighbor();
    }
    */

    if(i % term == 0){
 //     fout << t << " " << gamma.temperature() << " " << gamma.radius_of_gyration() << std::endl;

			mid = std::chrono::system_clock::now();
      progress=t/tmax;
      std::cout << 100*progress << "%    " <<  "Estimated remainding time is " << (1.0-progress)*std::chrono::duration_cast<std::chrono::seconds>(mid-start).count()/progress << "seconds     \r" << std::flush;
    }   
  }
  gamma.cdview_output(fout);
	
  fout.close();	
	end = std::chrono::system_clock::now();

	std::cout << "\n" << "Calculation ended. Spended " << std::chrono::duration_cast<std::chrono::seconds>(end-start).count() << " seconds" << std::endl;

	return 0;
}	

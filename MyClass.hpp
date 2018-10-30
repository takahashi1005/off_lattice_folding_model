#include <cmath>
#include <random>
#include "mt.h"
#include "MyClass.hpp"

//----------------------------------------------------------------------
// Initialize
//----------------------------------------------------------------------


MyClass::MyClass(){

	std::random_device rnd;
	MT::SetSeed(rnd());

  std::vector<double> q_previous(dim,0.0);
	double v;
	double ss=0.0;

	ExternalForce=0.0;
  x0=0.0;

  ATP = unbound;

  const monomer test;
  Protein.assign(ProteinLength, test);

// 座標、運動量、配列初期化
  Protein[0].c = Sequence[0];
  for(int n=1;n<ProteinLength;++n){

    q_previous[0] += 1.0;
    Protein[n].q = q_previous;

    for(auto &x :Protein[n].p){
			x=MT::GetGauss();
			ss+=x*x;
		}
    Protein[n].c = Sequence[n];

	}

// 運動量規格化
	v=sqrt(2.0*energy_per_atom*(ProteinLength-1)/ss);
	for(int n=1;n < ProteinLength; ++n){
		for(auto &x : Protein[n].p){
			x=v*x;
		}
	}

}

//----------------------------------------------------------------------
// Increment  Book
//----------------------------------------------------------------------

void MyClass::IncrementNeighbor(){
  int n,m;

  for(n=0;n<ProteinLength;++n){
    Protein[n].neighbor.clear();
    Protein[n].neighbor.shrink_to_fit();
    for(m=n+1;m<ProteinLength;++m){
      if(abs(Protein[n].q - Protein[m].q) < range_neighborhood){
        Protein[n].neighbor.push_back(m);
        Protein[m].neighbor.push_back(n);
      }
    }
  }

}

//----------------------------------------------------------------------
// Temperature
//----------------------------------------------------------------------

double MyClass::temperature(){
	double ss=0.0;

	for(int n=1;n < ProteinLength; ++n){
		for(auto &x : Protein[n].p){
			ss+=x*x;
		}
	}

	return ss/(double)ProteinLength/(double)dim;
}

//----------------------------------------------------------------------
// Radius of Gyration 
//----------------------------------------------------------------------

double MyClass::radius_of_gyration(){
  double R_g2=0.0;
  std::vector<double> q_g(dim,0.0);
  int n,i;

  for(auto &a : Protein){
    q_g = q_g + a.q;
  }

  q_g = q_g / (double)ProteinLength;
  
  for(auto &a : Protein){
    R_g2 += pow(abs(a.q - q_g), 2);
  }

  return R_g2/ProteinLength;
}

//----------------------------------------------------------------------
// Calculate the bond_interaction from neighborhood 
//----------------------------------------------------------------------

double MyClass::bond_interaction(int n, int i){
  
  return - k * (Protein[n].q[i]-Protein[n-1].q[i]) * (1.0 - 1.0 / abs(Protein[n].q - Protein[n-1].q));
  
}

//----------------------------------------------------------------------
// Lenord Jones Potensital 
//----------------------------------------------------------------------

double MyClass::HPN_interaction(int n,int i){
  double force = 0.0;

  for(int m = 0;m < ProteinLength;++m){
    if(m == n) continue;
    double dis = abs(Protein[n].q - Protein[m].q);
    if(Protein[n].c == H && Protein[m].c == H){
      force += 4 * epsilon_h * ( 12 * pow(sigma, 12) / pow(dis, 14) - 6 * pow(sigma, 6) / pow(dis, 8)) * (Protein[n].q[i] - Protein[m].q[i]);
    } else if(Protein[n].c == N || Protein[m].c == N){
      force += 48 * epsilon_h * pow(sigma, 12) / pow(dis, 14) * (Protein[n].q[i] - Protein[m].q[i]);
    } else {
      force += 4 * epsilon_p * ( 12 * pow(sigma, 12) / pow(dis, 14) + 6 * pow(sigma, 6) / pow(dis, 8)) * (Protein[n].q[i] - Protein[m].q[i]);
    }
  }

	return force;
}

//----------------------------------------------------------------------
// Equation of Motion
//----------------------------------------------------------------------

double MyClass::F(double t, int n, int i) {
    return Protein[n].p[i];
}

double MyClass::G( double t, int n, int i) {
  if (n == ProteinLength-1) {
		if(i == 0){
      return  bond_interaction(n,0)/* + HPN_interaction(n,i) + ExternalForce*/ - Protein[n].p[i];
    } else {
      return bond_interaction(n,i)/* + HPN_interaction(n,i)*/ - Protein[n].p[i];
    }
  }	else {
    return bond_interaction(n,i) - bond_interaction(n+1,i)/* + HPN_interaction(n,i)*/ - Protein[n].p[i];
  }
}

//----------------------------------------------------------------------
// Euler
//----------------------------------------------------------------------

void MyClass::Euler(double t){

	double qnext[ProteinLength][dim];
	double pnext[ProteinLength][dim];

  for(int n=1;n<ProteinLength;++n){
    for(int i=0;i<dim;++i){
      qnext[n][i] = Protein[n].q[i] + dt * F(t,n,i);
			pnext[n][i] = Protein[n].p[i] + dt * G(t,n,i) + MT::GetGauss() * sqrt(2*D*dt);
    }
  }
  for(int n=1;n<ProteinLength;++n){
    for(int i=0;i<dim;++i){
      Protein[n].q[i] = qnext[n][i];
			Protein[n].p[i] = pnext[n][i];
    }
  }

}

//----------------------------------------------------------------------
// LeapFrog
//----------------------------------------------------------------------

void MyClass::LeapFrog_q(double t){
  for(int n=0;n<ProteinLength;++n){
 		for(int i=0;i<dim;i++){
			Protein[n].q[i]=Protein[n].q[i]+2*dt*F(t,n,i);
		}
  }
}

void MyClass::LeapFrog_p(double t){
  for(int n=0;n<ProteinLength;++n){
		for(int i=0;i<dim;++i){
			Protein[n].p[i]=Protein[n].p[i]+2*dt*G(t,n,i);
		}
  }
}

//----------------------------------------------------------------------
// RungeKutta
//----------------------------------------------------------------------

void MyClass::RungeKutta(double t){
	double q_next[ProteinLength][dim]={},p_next[ProteinLength][dim]={};
	double q_pre[ProteinLength][dim],p_pre[ProteinLength][dim];
	double q_rk[ProteinLength][dim],p_rk[ProteinLength][dim];
	int n,i;

	for(n=0;n<ProteinLength;++n){
		for(i=0;i<dim;++i){
			q_pre[n][i]=Protein[n].q[i];
			p_pre[n][i]=Protein[n].p[i];
		}
	}

// runge kutta step1
	for(n=0;n<ProteinLength;++n){
		for(i=0;i<dim;++i){
			q_rk[n][i]=q_pre[n][i]+dt*0.5*F(t,n,i);
			p_rk[n][i]=p_pre[n][i]+dt*0.5*G(t,n,i);

			Protein[n].q[i]=q_rk[n][i];
			Protein[n].p[i]=p_rk[n][i];

			q_next[n][i]+=2.0*q_rk[n][i];
			p_next[n][i]+=2.0*p_rk[n][i];

		}
	}

// runge kutta step2
	for(n=0;n<ProteinLength;++n){
		for(i=0;i<dim;++i){
			q_rk[n][i]=q_pre[n][i]+dt*0.5*F(t+0.5*dt,n,i);
			p_rk[n][i]=p_pre[n][i]+dt*0.5*G(t+0.5*dt,n,i);

			Protein[n].q[i]=q_rk[n][i];
			Protein[n].p[i]=p_rk[n][i];
		
			q_next[n][i]+=4.0*q_rk[n][i];
			p_next[n][i]+=4.0*p_rk[n][i];
		}
	}

// runge kutta step3
	for(n=0;n<ProteinLength;++n){
		for(i=0;i<dim;++i){
			q_rk[n][i]=q_pre[n][i]+dt*F(t+0.5*dt,n,i);
			p_rk[n][i]=p_pre[n][i]+dt*G(t+0.5*dt,n,i);

			Protein[n].q[i]=q_rk[n][i];
			Protein[n].p[i]=p_rk[n][i];

			q_next[n][i]+=2.0*q_rk[n][i];
			p_next[n][i]+=2.0*p_rk[n][i];
		}
	}
	
// runge kutta step4
	for(n=0;n<ProteinLength;++n){
		for(i=0;i<dim;++i){
			q_rk[n][i]=-2.0*q_pre[n][i]+dt*F(t+dt,n,i);
			p_rk[n][i]=-2.0*p_pre[n][i]+dt*G(t+dt,n,i);

			q_next[n][i]+=q_rk[n][i];
			p_next[n][i]+=p_rk[n][i];
		}
	}

//  increment
	for(n=0;n<ProteinLength;++n){
		for(i=0;i<dim;++i){
			Protein[n].q[i]=q_next[n][i]/6.0;
			Protein[n].p[i]=p_next[n][i]/6.0;
		}
  }
}

//----------------------------------------------------------------------
// For test
//----------------------------------------------------------------------
double MyClass::test(){
  double Energy=0.0;
  for(int n=0;n<ProteinLength;++n){
    Energy += 0.5 * IP(Protein[n].p, Protein[n].p);
    if(n != ProteinLength - 1){
      Energy += 0.5 * k * pow(abs(Protein[n].q - Protein[n+1].q) - 1.0, 2);
    }
  }
  return Energy/(double)(ProteinLength -1);
}

void MyClass::cdview_output(std::ostream &os){
  for(int i=0;i < ProteinLength; ++i){
    os << i << " ";
    if(Protein[i].c == H) os << "0" << " ";
    else if(Protein[i].c == P) os << "2" << " ";
    else if(Protein[i].c == N) os << "1" << " ";
    os << Protein[i].q << std::endl;
  }
}

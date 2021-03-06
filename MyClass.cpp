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
  double z;
  double phi;

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

    z = 2.0 * MT::GetDouble() - 1.0;
    phi = 2.0 * M_PI * MT::GetDouble();

    q_previous[0] += sqrt(1 - z*z) * cos(phi);
    q_previous[1] += sqrt(1 - z*z) * sin(phi);
    q_previous[2] += z;

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
// Angle
//----------------------------------------------------------------------

double MyClass::angle_interaction(int n, int i){
  double force = 0.0;

  if( n == 1){
    std::vector<double> v2 = Protein[n].q - Protein[n-1].q;
    std::vector<double> v3 = Protein[n+1].q - Protein[n].q;
    std::vector<double> v4 = Protein[n+2].q - Protein[n+1].q;

    double a2 = abs(v2);
    double a3 = abs(v3);
    double a4 = abs(v4);

    double ip23 = - IP(v2, v3);
    double ip34 = - IP(v3, v4);

    double theta23 = acos(ip23 / a2 / a3);
    double theta34 = acos(ip34 / a3 / a4);

    double co23 = k_theta * (theta23 - theta_zero) / sqrt( pow(a2 * a3, 2.0) - ip23 * ip23);
    double co34 = k_theta * (theta34 - theta_zero) / sqrt( pow(a3 * a4, 2.0) - ip34 * ip34);

    force += - co23 * (v3[i] - ip23 / a2 / a2 * (-1.0) * v2[i]);
    force += - co23 * ((-1.0) * v2[i] - ip23 / a3 / a3 * v3[i]);
    force += co34 * (v4[i] - ip34 / a3 / a3 *(-1.0) * v3[i]);

  } else if( n == ProteinLength - 2){
    std::vector<double> v1 = Protein[n-1].q - Protein[n-2].q;
    std::vector<double> v2 = Protein[n].q - Protein[n-1].q;
    std::vector<double> v3 = Protein[n+1].q - Protein[n].q;

    double a1 = abs(v1);
    double a2 = abs(v2);
    double a3 = abs(v3);

    double ip12 = - IP(v1, v2);
    double ip23 = - IP(v2, v3);

    double theta12 = acos(ip12 / a1 / a2);
    double theta23 = acos(ip23 / a2 / a3);

    double co12 = k_theta * (theta12 - theta_zero) / sqrt( pow(a1 * a2, 2.0) - ip12 * ip12);
    double co23 = k_theta * (theta23 - theta_zero) / sqrt( pow(a2 * a3, 2.0) - ip23 * ip23);

    force += co12 * ((-1.0)*v1[i] - ip12 / a2 / a2 * v2[i]);
    force += - co23 * (v3[i] - ip23 / a2 / a2 * (-1.0) * v2[i]);
    force += - co23 * ((-1.0) * v2[i] - ip23 / a3 / a3 * v3[i]);

  } else if( n == ProteinLength - 1){
    std::vector<double> v1 = Protein[n-1].q - Protein[n-2].q;
    std::vector<double> v2 = Protein[n].q - Protein[n-1].q;

    double a1 = abs(v1);
    double a2 = abs(v2);

    double ip12 = - IP(v1, v2);

    double theta12 = acos(ip12 / a1 / a2);

    double co12 = k_theta * (theta12 - theta_zero) / sqrt( pow(a1 * a2, 2.0) - ip12 * ip12);

    force += co12 * ((-1.0)*v1[i] - ip12 / a2 / a2 * v2[i]);

  } else {
    std::vector<double> v1 = Protein[n-1].q - Protein[n-2].q;
    std::vector<double> v2 = Protein[n].q - Protein[n-1].q;
    std::vector<double> v3 = Protein[n+1].q - Protein[n].q;
    std::vector<double> v4 = Protein[n+2].q - Protein[n+1].q;

    double a1 = abs(v1);
    double a2 = abs(v2);
    double a3 = abs(v3);
    double a4 = abs(v4);

    double ip12 = - IP(v1, v2);
    double ip23 = - IP(v2, v3);
    double ip34 = - IP(v3, v4);

    double theta12 = acos(ip12 / a1 / a2);
    double theta23 = acos(ip23 / a2 / a3);
    double theta34 = acos(ip34 / a3 / a4);

    double co12 = k_theta * (theta12 - theta_zero) / sqrt( pow(a1 * a2, 2.0) - ip12 * ip12);
    double co23 = k_theta * (theta23 - theta_zero) / sqrt( pow(a2 * a3, 2.0) - ip23 * ip23);
    double co34 = k_theta * (theta34 - theta_zero) / sqrt( pow(a3 * a4, 2.0) - ip34 * ip34);

    force += co12 * ((-1.0)*v1[i] - ip12 / a2 / a2 * v2[i]);
    force += - co23 * (v3[i] - ip23 / a2 / a2 * (-1.0) * v2[i]);
    force += - co23 * ((-1.0) * v2[i] - ip23 / a3 / a3 * v3[i]);
    force += co34 * (v4[i] - ip34 / a3 / a3 *(-1.0) * v3[i]);

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
  double force = 0;
  if (n == ProteinLength-1) {
		if(i == 0){
      force +=  bond_interaction(n,0);
      force += HPN_interaction(n,i);
      force += ExternalForce;
      force += - Protein[n].p[i];
      force += angle_interaction(n,i);
    } else {
      force += bond_interaction(n,i);
      force += HPN_interaction(n,i);
      force += - Protein[n].p[i];
      force += angle_interaction(n,i);
    }
  }	else {
    force += bond_interaction(n,i);
    force += -bond_interaction(n+1,i);
    force += HPN_interaction(n, i);
    force += - Protein[n].p[i];
    force += angle_interaction(n, i);
  }
  return force;
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

/*
 * pendulate.cc
 *
 *  created on: 21.01.2016
 *      author: rungger
 */

/*
 * information about this example is given in the readme file
 *
 */

#include <array>
#include <iostream>

#include "cuddObj.hh"

#include "SymbolicSet.hh"
#include "SymbolicModelGrowthBound.hh"

#include "TicToc.hh"
#include "RungeKutta4.hh"
#include "FixedPoint.hh"


// bdd output files
std::string bddDir = "./bdd/";
std::string stateSpaceFile = bddDir + "pendulate_ss.bdd";
std::string obstaclesFile = bddDir + "pendulate_obst.bdd";
std::string targetFile = bddDir + "pendulate_target.bdd";
std::string controllerFile = bddDir + "pendulate_controller.bdd";
std::string transitionRelationFile = bddDir + "pendulate_transition_relation.bdd";

/* state space dim */
#define sDIM 3
#define iDIM 1

/* data types for the ode solver */
typedef std::array<double,3> state_type;
typedef std::array<double,1> input_type;

/* sampling time */
const double tau = 0.1;
/* number of intermediate steps in the ode solver */
const int nint=10;
OdeSolver ode_solver(sDIM,nint,tau);

/* system parameters */
const double A = 0.1; // Amplitude of vertical driving motion
const double L = 1; // Length of the pendulum
const double g = 9.81; // Acceleration of gravity
const double b = 0.05; // Linear drag term

const double alpha = g / L;
const double beta = A / L;

/* we integrate the pendulate ode by tau sec (the result is stored in x)  */
auto  pendulate_post = [](state_type &x, input_type &u) -> void {

  /* the ode describing the pendulate */
  auto rhs =[](state_type& dxdt,  const state_type &x, input_type &u) -> void {
      dxdt[0] = x[1];
      dxdt[1] = (beta*u[0]*u[0]*std::cos(x[2]) - alpha)*std::sin(x[0]) - b*x[1];
      dxdt[2] = u[0];
  };
  ode_solver(rhs,x,u);
  
  x[0] = x[0] - 2*M_PI*std::floor(x[0] / (2*M_PI)); // Wrap theta around [0, 2pi]
  x[2] = x[2] - 2*M_PI*std::floor(x[2] / (2*M_PI)); // Wrap phase around [0, 2pi]
};

/* computation of the growth bound (the result is stored in r)  */
auto radius_post = [](state_type &r, input_type &u) -> void {
    r[0] = r[0]+r[2]*std::abs(u[0])*0.3;
    r[1] = r[1]+r[2]*std::abs(u[0])*0.3;
};


/* forward declaration of the functions to setup the state space 
 * and input space of the pendulate example */
scots::SymbolicSet pendulateCreateStateSpace(Cudd &mgr);
scots::SymbolicSet pendulateCreateInputSpace(Cudd &mgr);

int main() {
  /* to measure time */
  TicToc tt;
  /* there is one unique manager to organize the bdd variables */
  Cudd mgr;

  /****************************************************************************/
  /* construct SymbolicSet for the state space */
  /****************************************************************************/
  scots::SymbolicSet ss=pendulateCreateStateSpace(mgr);
  ss.writeToFile(stateSpaceFile.c_str());
  /* write SymbolicSet of obstacles to pendulate_obst.bdd */
  // ss.complement();
  // ss.writeToFile(obstaclesFile.c_str());
  // ss.complement();
  std::cout << "Unfiorm grid details:" << std::endl;
  ss.printInfo(1);

  /****************************************************************************/
  /* the target set */
  /****************************************************************************/
  /* first make a copy of the state space so that we obtain the grid
   * information in the new symbolic set */
  scots::SymbolicSet ts = ss;
  /* define the target set as a symbolic set */
  double H[4*sDIM]={-1, 0, 0, // This matrix lets us simply define bounds on x[0] and x[1]
                    1, 0, 0,
                    0,-1, 0,
                    0, 1, 0};
  /* compute inner approximation of P={ x | H x<= h1 }  */
  double h[4] = {-3,3.3, 0.1,0.1}; // {-lb0, ub0, -lb1, up1}
  ts.addPolytope(4,H,h, scots::INNER);
  ts.writeToFile(targetFile.c_str());


  /****************************************************************************/
  /* construct SymbolicSet for the input space */
  /****************************************************************************/
  scots::SymbolicSet is=pendulateCreateInputSpace(mgr);
  std::cout << std::endl << "Input space details:" << std::endl;
  is.printInfo(1);

  /****************************************************************************/
  /* setup class for symbolic model computation */
  /****************************************************************************/
  /* first create SymbolicSet of post variables 
   * by copying the SymbolicSet of the state space and assigning new BDD IDs */
  scots::SymbolicSet sspost(ss,1);
  /* instantiate the SymbolicModel */
  scots::SymbolicModelGrowthBound<state_type,input_type> abstraction(&ss, &is, &sspost);
  /* compute the transition relation */
  tt.tic();
  abstraction.computeTransitionRelation(pendulate_post, radius_post);
  std::cout << std::endl;
  tt.toc();
  /* get the number of elements in the transition relation */
  std::cout << std::endl << "Number of elements in the transition relation: " << abstraction.getSize() << std::endl;


  /****************************************************************************/
  /* we continue with the controller synthesis */
  /****************************************************************************/
  int verbose=1;
  /* we setup a fixed point object to compute reachabilty controller */
  scots::FixedPoint fp(&abstraction);
  /* the fixed point algorithm operates on the BDD directly */
  BDD T = ts.getSymbolicSet();
  tt.tic();
  BDD C = fp.reach(T,verbose);
  tt.toc();

  /****************************************************************************/
  /* last we store the controller as a SymbolicSet 
   * the underlying uniform grid is given by the Cartesian product of 
   * the uniform gird of the space and uniform gird of the input space */
  /****************************************************************************/
  scots::SymbolicSet controller(ss,is);
  controller.setSymbolicSet(C);
  controller.writeToFile(controllerFile.c_str());

  scots::SymbolicSet tr = abstraction.getTransitionRelation();
  tr.writeToFile(transitionRelationFile.c_str());
  
  return 0;
}

scots::SymbolicSet pendulateCreateStateSpace(Cudd &mgr) {

  const double maxOmega = 11;

  /* setup the workspace of the synthesis problem and the uniform grid */
  /* lower bounds of the hyper rectangle */
  double lb[sDIM]={0,-maxOmega, 0};
  /* upper bounds of the hyper rectangle */
  double ub[sDIM]={2*M_PI,maxOmega, 2*M_PI}; 
  /* grid node distance diameter */
  double eta[sDIM]={.1, .5, .2};
  // double eta[sDIM]={1, 5, 1};



  /* eta is added to the bound so as to ensure that the whole
   * [0,10]x[0,10]x[-pi-eta,pi+eta] is covered by the cells */

  scots::SymbolicSet ss(mgr,sDIM,lb,ub,eta);

  /* add the grid points to the SymbolicSet ss */
  ss.addGridPoints();

 return ss;
}

scots::SymbolicSet pendulateCreateInputSpace(Cudd &mgr) {

  /* lower bounds of the hyper rectangle */
  double lb[sDIM]={1};  
  /* upper bounds of the hyper rectangle */
  double ub[sDIM]={50}; 
  /* grid node distance diameter */
  double eta[sDIM]={1};   
  // double eta[sDIM]={5};   

  scots::SymbolicSet is(mgr,iDIM,lb,ub,eta);
  is.addGridPoints();

  return is;
}


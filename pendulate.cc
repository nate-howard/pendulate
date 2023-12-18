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
std::string bddDir = "./bdd/stay_";
std::string stateSpaceFile = bddDir + "pendulate_ss.bdd";
std::string targetFile = bddDir + "pendulate_target.bdd";
std::string controllerFile = bddDir + "pendulate_controller.bdd";
std::string transitionRelationFile = bddDir + "pendulate_transition_relation.bdd";

/* state space dim */
#define sDIM 2
#define iDIM 1

/* data types for the ode solver */
typedef std::array<double,2> state_type;
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
// const double b = 0.05; // Linear drag term

const double alpha = g / L;
const double beta = A * A / (4 * L * L);

/* we integrate the pendulate ode by tau sec (the result is stored in x)  */
auto  pendulate_post = [](state_type &x, input_type &u) -> void {

  /* the ode describing the pendulate */
  auto rhs =[](state_type& dxdt,  const state_type &x, input_type &u) -> void {
      dxdt[0] = x[1];
      dxdt[1] = -alpha*std::sin(x[0]) - beta*u[0]*u[0]*std::sin(2*x[0]);
  };
  ode_solver(rhs,x,u);
  
  x[0] = x[0] - 2*M_PI*std::floor(x[0] / (2*M_PI)); // Wrap theta around [0, 2pi]
};

/* computation of the growth bound (the result is stored in r)  */
auto radius_post = [](state_type &r, state_type &x, input_type &u) -> void {
  bool print = ((double) std::rand() / (RAND_MAX)) < (1/500000.0);
  if(print) std::cout << "r0=[" << r[0] << "," << r[1] << "], x=[" << x[0] << "," << x[1] << "], u = " << u[0] << std::endl;
  // double bw = beta * u[0];
  double new_r1 = r[1] + r[0]*abs(alpha*std::cos(x[0]) + 2*beta*u[0]*u[0]*std::cos(2*x[0]))*tau;
  r[0] = r[0] + new_r1*tau;
  r[1] = new_r1;
  if(print) std::cout << "r1 = [" << r[0] << ", " << r[1] << "]" << std::endl;
};


/* functions to setup the state space 
 * and input space of the pendulate example */
scots::SymbolicSet pendulateCreateStateSpace(Cudd &mgr) {

  const double maxOmega = 8;

  /* setup the workspace of the synthesis problem and the uniform grid */
  /* lower bounds of the hyper rectangle */
  double lb[sDIM]={0,-maxOmega};
  /* upper bounds of the hyper rectangle */
  double ub[sDIM]={2*M_PI,maxOmega}; 
  /* grid node distance diameter */
  double eta[sDIM]={.05, .05};
  // double eta[sDIM]={.2, .5, .5};

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

int main2() {
  state_type x = {2.5, 0};
  input_type u = {50};

  int iter = 300;
  for(int i=0; i<iter; i++) {
    pendulate_post(x, u);
    std::cout << "[" << x[0] << ", " << x[1] << ", "<< x[2] << "]" << std::endl;
  }

  return 0;
}

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

  std::cout << "Unfiorm grid details:" << std::endl;
  ss.printInfo(1);

  /****************************************************************************/
  /* the target set */
  /****************************************************************************/
  /* first make a copy of the state space so that we obtain the grid
   * information in the new symbolic set */

  scots::SymbolicSet ts = ss;
  /* define the target set as a symbolic set */
  double H[4*sDIM]={-1, 0, // This matrix lets us simply define bounds on x[0] and x[1]
                    1, 0,
                    0,-1,
                    0, 1};
  /* compute inner approximation of P={ x | H x<= h1 }  */
  double h[4] = {-2.5,2.9, 0.2,0.2}; // {-lb0, ub0, -lb1, up1}
  ts.addPolytope(4,H,h, scots::INNER);
  ts.writeToFile(targetFile.c_str());
  std::cout << std::endl << "Target set info:" << std::endl;
  ts.printInfo(1);


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
  // int verbose=1;
  /* we setup a fixed point object to compute reachabilty controller */
  scots::FixedPoint fp(&abstraction);
  /* the fixed point algorithm operates on the BDD directly */
  BDD T = ts.getSymbolicSet();

  // Reach synthesis
  // tt.tic();
  // BDD C = fp.safe(T,1);
  // BDD C = fp.pre(T);
  // for(int i=1; i<=100; i++) {
  //   BDD temp = fp.pre(C);
  //   if(temp.IsZero()) {
  //     std::cout << "Set is empty after " << i+1 << " iterations!!!" << std::endl;
  //     break;
  //   } else {
  //     C = temp;
  //   }
  // }
  // std::cout << "Empty: " << C.IsZero() << std::endl;
  // C = fp.pre(C);
  // std::cout << "Empty: " << C.IsZero() << std::endl;
  // C = fp.pre(C);
  // std::cout << "Empty: " << C.IsZero() << std::endl;
  // // C = fp.pre(C);
  // // std::cout << "Empty: " << C.IsZero() << std::endl;
  // tt.toc();

  // Reach and stay synthesis
  tt.tic();
  size_t i,j;
  /* outer fp*/
  BDD X=mgr.bddOne();
  BDD XX=mgr.bddZero();
  /* inner fp*/
  BDD Y=mgr.bddZero();
  BDD YY=mgr.bddOne();
  /* the controller */
  BDD C=mgr.bddZero();
  BDD U=is.getCube();
  /* as long as not converged */
  for(i=1; XX != X; i++) {
    X=XX;
    BDD preX=fp.pre(X);
    /* init inner fp */
    YY = mgr.bddOne();
    for(j=1; YY != Y; j++) {
      Y=YY;
      YY= ( fp.pre(Y) & T ) | preX;
    }
    XX=YY;
    std::cout << "Iterations inner: " << j << std::endl;
    /* remove all (state/input) pairs that have been added
     * to the controller already in the previous iteration * */
    BDD N = XX & (!(C.ExistAbstract(U)));
    /* add the remaining pairs to the controller */
    C=C | N;
    //std::cout << C.CountMinterm(17) << std::endl;
  }
  std::cout << "Iterations outer: " << i << std::endl;
  tt.toc();

  /****************************************************************************/
  /* last we store the controller as a SymbolicSet 
   * the underlying uniform grid is given by the Cartesian product of 
   * the uniform gird of the space and uniform gird of the input space */
  /****************************************************************************/
  scots::SymbolicSet controller(ss,is);
  controller.setSymbolicSet(C);
  std::cout << "Domain size: " << controller.getSize() << std::endl;
  controller.writeToFile(controllerFile.c_str());

  scots::SymbolicSet tr = abstraction.getTransitionRelation();
  tr.writeToFile(transitionRelationFile.c_str());
  
  return 0;
}

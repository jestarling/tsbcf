#include <RcppArmadillo.h>

#include <cmath>
#include "funs.h"
#include "rng.h"
#include <map>
#ifdef MPIBART
#include "mpi.h"
#endif

using Rcpp::Rcout;
using namespace arma;
using namespace Rcpp;


//-------------------------------------------------------------
// Squared Exponential Covariance Function for two time vectors (x,y).
//-------------------------------------------------------------
mat cov_se(vec t1, vec t2, double ls, double var)
{
   //-------------------------------------------------------------
   // INPUTS:	   x,y   = two vectors from the same space.
   //				ls    = b 		= length (of period)
   //				var   = tau1.sq = variance of function
   //-------------------------------------------------------------
   // OUTPUT:	The squared exponential covariance matrix.
   //-------------------------------------------------------------
   double arg;
   int n1 = t1.size();
   int n2 = t2.size();
   mat C(n1,n2);
   for(int i = 0; i < n1; i++) {
      for(int j=0; j < n2; j++) {
         arg = (t1[i] - t2[j])/ls;
         C(i,j) = var*exp(-0.5*arg*arg);
      }
   }

   // Add jitter to diagonal to ensure is pos semidef.
   C.diag() = C.diag() + .000001;

   return C;
}

//-------------------------------------------------------------
// Utility function for calculating posterior MVN params for N(Phi^(-1)*m, Phi^(-1))
//-------------------------------------------------------------
List mvn_post_util(double sigma, vec mu0, mat Prec0, vec n_vec, vec sy_vec){

   // Initialize matrices and vectors.
   double s2 = sigma * sigma;
   mat Lam = diagmat(n_vec) / s2;
   vec m = (1/s2) * sy_vec + Prec0 * mu0;  // K = Prec0

   return List::create(
      _["Phi"] = Prec0 + Lam,
      _["m"]   = m
   ) ;
}

List mvn_post_util_het(vec mu0, mat Prec0, vec n0_vec, vec n_vec, vec sy_vec){
   //n0_vec for het is vector of nt for each time (not weighted).
   //n_vec for het is vector of sum(w_i) precisions at each time.
   //sy_vec for het is vector of sum(y_i*w_i) at each time.

   // Initialize matrices and vectors.
   mat Lam = diagmat(n_vec) + Prec0;
   vec m = sy_vec + Prec0 * mu0;

   return List::create(
      _["Phi"] = Lam,
      _["m"]   = m
   ) ;
}

//-------------------------------------------------------------
// Generates realizations from multivariate normal.
//-------------------------------------------------------------
mat rmvnormArma(int n, vec mu, mat sigma) {
   //-------------------------------------------------------------
   // INPUTS:	   n = sample size
   //				   mu = vector of means
   //				   sigma = covariance matrix
   //-------------------------------------------------------------
   // OUTPUT:	n realizations of the specified MVN.
   //-------------------------------------------------------------
   int ncols = sigma.n_cols;
   mat Y = randn(n, ncols);
   mat result = (repmat(mu, 1, n).t() + Y * chol(sigma)).t();
   return result;
}

//-------------------------------------------------------------------------------
// log of the integrated likelihood for tsbart, for a given tree/leaf.
//-------------------------------------------------------------------------------
double lil_ts(vec nt, vec sy_vec, double sy2, double sigma, vec mu0, mat Prec0){
   // nt = vector of number of obs in each time point for the given tree/leaf. nt = [nl_{t=1}, ..., nl_{t=T}]
   // sy = vector of sums of y's at each time point.  sy = [ sum(y in t=1), ..., sum(y in t=T) ]
   // sigma = error sd, sqrt of sigma^2
   // mu0 = vector of prior means.
   // Prec0 = prior precision matrix for means (from sq exp kernel)

   // For computational efficiency, we leave out the -.5*(mu_0^T K mu_0) term, since
   // we let mu_0=0.  Add this term if mu0 != 0.

   double sig2 = sigma*sigma;    //sigma^2
   double nl = sum(nt);          // Total number of yl in leaf.

   // Precache a few terms to make calculation faster.
   vec b = sy_vec/sig2 + Prec0*mu0;
   mat C = Prec0;
   C.diag() = C.diag() + nt/sig2;

   // Calculate log-likelihood.  Note: mu0.t()*K*mu0 excluded as mu0=0.
   double ll = -.5*nl*log(2*PI*sig2) + .5*log(det(Prec0)) - .5*log(det(C)) -
                .5*as_scalar(sy2/sig2 - b.t()*C.i()*b);

   return(ll);
}

// For het variances.
double lilhet_ts(double n0, double n, vec n_vec, vec sy_vec, double sy2, vec mu0, mat Prec0){
   // n0 = number of obs in given tree/leaf.
   // n = sum of log-precisions phi in given tree/leaf.
   // n_vec = vector of sum of phi's for het, at each time point, for given tree/leaf.
   // sy = vector of sums of y's * phi's at each time point.  sy = [ sum(phi*y in t=1), ..., sum(phi*y in t=T) ]
   // sy2 = scalar, sum of y*y*phi for all obs in given tree/leaf.
   // mu0 = vector of prior means.
   // Prec0 = prior precision matrix for means (from sq exp kernel)

   // For computational efficiency, we leave out the -.5*(mu_0^T K mu_0) term, since
   // we let mu_0=0.  Add this term if mu0 != 0.

   // Precache a few terms to make calculation faster.
   mat C = Prec0; // K = Prec0
   C.diag() += n_vec; // Add sums of precisions to diagonal.
   vec b = sy_vec + Prec0 * mu0;

   // Calculate log-likelihood.  Note: mu0.t()*K*mu0 excluded as cancels in ratios.
   double ll = - .5*n0*log(2*PI)
               + .5*n // This is the .5 * log(det(Lambda)) term, where Lambda=diag(w).
               + .5*log(det(Prec0))
               - .5*log(det(C))
               - as_scalar(.5*(sy2 - b.t()*C.i()*b));

   return(ll);
}


//-------------------------------------------------------------------------------
// NEW: tsbart: get sufficients stats for all bottom nodes
//-------------------------------------------------------------------------------

void allsuff_ts(tree& x, xinfo& xi, dinfo& di, tree::npv& bnv, std::vector<sinfo>& sv)
{
   // Bottom nodes are written to bnv.
   // Suff stats for each bottom node are written to elements (each of class sinfo) of sv.
   // Initialize data structures
   tree::tree_cp tbn; //the pointer to the bottom node for the current observations.  tree_cp bc not modifying tree directly.
   size_t ni;         //the  index into vector of the current bottom node
   double *xx;        //current x
   double y;          //current y
   double t;          //current t

   bnv.clear();      // Clear the bnv variable if any value is already saved there.
   x.getbots(bnv);   // Save bottom nodes for x to bnv variable.

   typedef tree::npv::size_type bvsz;  // Is a better C way to set type.  (tree::npv::size_type) will resolve to an integer,
                                       // or long int, etc.  We don't have to know that ahead of time by using this notation.
   bvsz nb = bnv.size();   // Initialize new var nb of type bvsz for number of bottom nodes, then...
   sv.resize(nb);          // Re-sizing suff stat vector to have same size as bottom nodes.

   // Resize vectors within sufficient stats to have di.tlen length.
   for(size_t i = 0; i < nb; ++i){
      sv[i].n_vec.resize(di.tlen);
      sv[i].sy_vec.resize(di.tlen);
   }

   // bnmap is a tuple (lookups, like in Python).  Want to index by bottom nodes.
   std::map<tree::tree_cp,size_t> bnmap;
   for(bvsz i=0;i!=bnv.size();i++) bnmap[bnv[i]]=i;  // bnv[i]
   //map looks like
   // bottom node 1 ------ 1
   // bottom node 2 ------ 2

   // Sum the y values (sy) and the y^2 values (sy2) for each node and store in sv.
   // Loop through each observation.  Push each obs x down the tree and find its bottom node,
   // then index into the suff stat for the bottom node corresponding to that obs.

   for(size_t i=0;i<di.n;i++) {
      xx = di.x + i*di.p;  //Index value: di.x is pointer to first element of n*p data vector.  Iterates through each element.
      y=di.y[i];           // Resolves to r.
      t = di.t[i];         // NEW: Resolves to current t

      tbn = x.bn(xx,xi); // Find bottom node for this observation.
      ni = bnmap[tbn];   // Map bottom node to integer index

      // Update the sufficient stats for the bottom node to which that obs belongs.
      ++(sv[ni].n);
      sv[ni].sy += y;
      sv[ni].sy2 += y*y;

      // Find index of matching time point.  Then adjust that entry for n_vec and sy_vec.
      uvec id = find(di.tref == t); // Idx of current obs t value.
      sv[ni].n_vec(id) += 1;
      sv[ni].sy_vec(id) += y;

   } // End obs loop.
}

// For het variances.
void allsuffhet_ts(tree& x, xinfo& xi, dinfo& di, double* phi, tree::npv& bnv, std::vector<sinfo>& sv)
{
   // phi are precisions for each observation.

   // Bottom nodes are written to bnv.
   // Suff stats for each bottom node are written to elements (each of class sinfo) of sv.
   // Initialize data structures
   tree::tree_cp tbn; //the pointer to the bottom node for the current observations.  tree_cp bc not modifying tree directly.
   size_t ni;         //the  index into vector of the current bottom node
   double *xx;        //current x
   double y;          //current y
   double t;          //current t

   bnv.clear();      // Clear the bnv variable if any value is already saved there.
   x.getbots(bnv);   // Save bottom nodes for x to bnv variable.

   typedef tree::npv::size_type bvsz;  // Is a better C way to set type.  (tree::npv::size_type) will resolve to an integer,
   // or long int, etc.  We don't have to know that ahead of time by using this notation.
   bvsz nb = bnv.size();   // Initialize new var nb of type bvsz for number of bottom nodes, then...
   sv.resize(nb);          // Re-sizing suff stat vector to have same size as bottom nodes.

   // Resize vectors within sufficient stats to have di.tlen length.
   for(size_t i = 0; i < nb; ++i){
      sv[i].n0_vec.resize(di.tlen);
      sv[i].n_vec.resize(di.tlen);
      sv[i].sy_vec.resize(di.tlen);

      // Fill with zeros, in case.
      sv[i].n0_vec.fill(0);
      sv[i].n_vec.fill(0);
      sv[i].sy_vec.fill(0);
   }

   // bnmap is a tuple (lookups, like in Python).  Want to index by bottom nodes.
   std::map<tree::tree_cp,size_t> bnmap;
   for(bvsz i=0;i!=bnv.size();i++) bnmap[bnv[i]]=i;  // bnv[i]
   //map looks like
   // bottom node 1 ------ 1
   // bottom node 2 ------ 2

   // Sum the y values (sy) and the y^2 values (sy2) for each node and store in sv.
   // Loop through each observation.  Push each obs x down the tree and find its bottom node,
   // then index into the suff stat for the bottom node corresponding to that obs.

   for(size_t i=0;i<di.n;i++) {
      xx = di.x + i*di.p;  //Index value: di.x is pointer to first element of n*p data vector.  Iterates through each element.
      y=di.y[i];           // Resolves to r.
      t = di.t[i];         // NEW: Resolves to current t

      tbn = x.bn(xx,xi); // Find bottom node for this observation.
      ni = bnmap[tbn];   // Map bottom node to integer index

      // Update the sufficient stats for the bottom node to which that obs belongs.
      sv[ni].n0 += 1;
      sv[ni].n += log(phi[i]); //sv[ni].n += phi[i];
      sv[ni].sy += phi[i]*y;
      sv[ni].sy2 += phi[i]*y*y;

      // Find index of matching time point.  Then adjust that entry for n_vec and sy_vec.
      uvec id = find(di.tref == t); // Idx of current obs t value.
      sv[ni].n0_vec(id) +=1;
      sv[ni].n_vec(id) += phi[i];
      sv[ni].sy_vec(id) += phi[i]*y;

   } // End obs loop.
}

//-------------------------------------------------------------------------------
// NEW: tsbart: get sufficient stats for children of node nx in tree x
// (for birth proposal)
//-------------------------------------------------------------------------------
// Birth proposal, homog variances.
void getsuff_ts(tree& x, tree::tree_cp nx, size_t v, size_t c, xinfo& xi, dinfo& di, sinfo& sl, sinfo& sr, size_t tlen)
{
   double *xx;//current x
   double y;  //current y
   double t;  //current t

   sl.n=0;sl.sy=0.0;sl.sy2=0.0;
   sr.n=0;sr.sy=0.0;sr.sy2=0.0;

   sl.n_vec = zeros(tlen); sl.sy_vec = zeros(tlen);
   sr.n_vec = zeros(tlen); sr.sy_vec = zeros(tlen);

   for(size_t i=0;i<di.n;i++) {
      xx = di.x + i*di.p;
      if(nx==x.bn(xx,xi)) { //does the bottom node = xx's bottom node

         y = di.y[i];   // extract current yi.  resolves to r.
         t = di.t[i];   // extract current ti

         if(xx[v] < xi[v][c]) { // Update left.
            sl.n++;
            sl.sy += y;
            sl.sy2 += y*y;

            // Find index of matching time point.  Then adjust that entry for n_vec and sy_vec.
            uvec id = find(di.tref == t); // Idx of current obs t value.
            sl.n_vec(id) += 1;
            sl.sy_vec(id) += y;

         } else { //Update right.
            sr.n++;
            sr.sy += y;
            sr.sy2 += y*y;

            // Find index of matching time point.  Then adjust that entry for n_vec and sy_vec.
            uvec id = find(di.tref == t); // Idx of current obs t value.
            sr.n_vec(id) += 1;
            sr.sy_vec(id) += y;

         }
      }
   }
}

// Birth proposal, het variances.
void getsuffhet_ts(tree& x, tree::tree_cp nx, size_t v, size_t c, xinfo& xi, dinfo& di, double* phi, sinfo& sl, sinfo& sr, size_t tlen)
{
   double *xx;//current x
   double y;  //current y
   double t;  //current t

   sl.n0=0;sl.n=0;sl.sy=0.0;sl.sy2=0.0;
   sr.n0=0;sr.n=0;sr.sy=0.0;sr.sy2=0.0;

   sl.n0_vec = zeros(tlen); sl.n_vec = zeros(tlen); sl.sy_vec = zeros(tlen);
   sr.n0_vec = zeros(tlen); sr.n_vec = zeros(tlen); sr.sy_vec = zeros(tlen);

   for(size_t i=0;i<di.n;i++) {
      xx = di.x + i*di.p;
      if(nx==x.bn(xx,xi)) { //does the bottom node = xx's bottom node

         y = di.y[i];   // extract current yi.
         t = di.t[i];   // extract current ti

         if(xx[v] < xi[v][c]) { // Update left.
            sl.n0 += 1;
            sl.n += log(phi[i]); //+= phi[i];
            sl.sy += phi[i]*y;
            sl.sy2 += phi[i]*y*y;

            // Find index of matching time point.  Then adjust that entry for n_vec and sy_vec.
            uvec id = find(di.tref == t); // Idx of current obs t value.
            sl.n0_vec(id) += 1;
            sl.n_vec(id) += phi[i];
            sl.sy_vec(id) += phi[i]*y;

         } else { //Update right.
            sr.n0 += 1;
            sr.n += log(phi[i]); //phi[i];
            sr.sy += phi[i]*y;
            sr.sy2 += phi[i]*y*y;

            // Find index of matching time point.  Then adjust that entry for n_vec and sy_vec.
            uvec id = find(di.tref == t); // Idx of current obs t value.
            sr.n0_vec(id) += 1;
            sr.n_vec(id) += phi[i];
            sr.sy_vec(id) += phi[i]*y;
         }
      }
   }
}

//--------------------------------------------------
//get sufficient stats for pair of bottom children nl(left) and nr(right) in tree x
// (for death proposal)
//--------------------------------------------------
// Death proposal, homog variance.
void getsuff_ts(tree& x, tree::tree_cp nl, tree::tree_cp nr, xinfo& xi, dinfo& di, sinfo& sl, sinfo& sr, size_t tlen)
{
   double *xx;//current x
   double y;  //current y
   double t;  //current t

   sl.n=0;sl.sy=0.0;sl.sy2=0.0;
   sr.n=0;sr.sy=0.0;sr.sy2=0.0;

   sl.n_vec = zeros(tlen); sl.sy_vec = zeros(tlen);
   sr.n_vec = zeros(tlen); sr.sy_vec = zeros(tlen);

   for(size_t i=0;i<di.n;i++) {
      xx = di.x + i*di.p;
      tree::tree_cp bn = x.bn(xx,xi);

      y = di.y[i];   // extract current yi.
      t = di.t[i];   // extract current ti

      uvec id = find(di.tref == t); // Idx of current obs t value.

      if(bn==nl) {
         sl.n++;
         sl.sy += y;
         sl.sy2 += y*y;

         // Find index of matching time point.  Then adjust that entry for n_vec and sy_vec.
         sl.n_vec(id) += 1;
         sl.sy_vec(id) += y;
      }

      if(bn==nr) {
         sr.n++;
         sr.sy += y;
         sr.sy2 += y*y;

         // Find index of matching time point.  Then adjust that entry for n_vec and sy_vec.
         sr.n_vec(id) += 1;
         sr.sy_vec(id) += y;
      }
   }
}

// For death proposal, het variances.
void getsuffhet_ts(tree& x, tree::tree_cp nl, tree::tree_cp nr, xinfo& xi, dinfo& di, double* phi, sinfo& sl, sinfo& sr, size_t tlen)
{
   double *xx;//current x
   double y;  //current y
   double t;  //current t

   sl.n0=0;sl.n=0;sl.sy=0.0;sl.sy2=0.0;
   sr.n0=0;sr.n=0;sr.sy=0.0;sr.sy2=0.0;

   sl.n0_vec = zeros(tlen); sl.n_vec = zeros(tlen); sl.sy_vec = zeros(tlen);
   sr.n0_vec = zeros(tlen); sr.n_vec = zeros(tlen); sr.sy_vec = zeros(tlen);

   for(size_t i=0;i<di.n;i++) {
      xx = di.x + i*di.p;
      tree::tree_cp bn = x.bn(xx,xi);

      y = di.y[i];   // extract current yi.
      t = di.t[i];   // extract current ti

      if(bn==nl) {
         y = di.y[i];
         sl.n0 += 1;
         sl.n += log(phi[i]); //phi[i];
         sl.sy += phi[i]*y;
         sl.sy2 += phi[i]*y*y;

         // Find index of matching time point.  Then adjust that entry for n_vec and sy_vec.
         uvec id = find(di.tref == t); // Idx of current obs t value.
         sl.n0_vec(id) += 1;
         sl.n_vec(id) += phi[i];
         sl.sy_vec(id) += phi[i]*y;
      }

      if(bn==nr) {
         y = di.y[i];
         sr.n0 += 1;
         sr.n += log(phi[i]); //phi[i];
         sr.sy += phi[i]*y;
         sr.sy2 += phi[i]*y*y;

         // Find index of matching time point.  Then adjust that entry for n_vec and sy_vec.
         uvec id = find(di.tref == t); // Idx of current obs t value.
         sr.n0_vec(id) += 1;
         sr.n_vec(id) += phi[i];
         sr.sy_vec(id) += phi[i]*y;
      }
   }
}

//--------------------------------------------------
// draw all the bottom node mu's

// For homog variances.
void drmu(tree& t, xinfo& xi, dinfo& di, pinfo& pi, RNG& gen)
{
   tree::npv bnv;
   std::vector<sinfo> sv(di.tlen);
   allsuff_ts(t,xi,di,bnv,sv);

   List post_pars;
   mat Phi;
   vec m;
   vec mu_draw;

   for(tree::npv::size_type i=0;i!=bnv.size();i++) {

      // Draw new mu value from MVN.
      post_pars = mvn_post_util(pi.sigma, pi.mu0, pi.Prec0, sv[i].n_vec, sv[i].sy_vec);

      Phi = Rcpp::as<arma::mat>(post_pars["Phi"]);
      m = Rcpp::as<arma::vec>(post_pars["m"]);
      mu_draw = rmvnorm_post(m, Phi);

      // Assign botton node values to new mu draw.
      bnv[i] -> setm(mu_draw);

      // Check for NA result.
      if(sum(bnv[i]->getm() == bnv[i]->getm()) == 0) {
          Rcpp::stop("drmu failed");
      }
   }
}

// For heterogeneous variances.
void drmuhet(tree& t, xinfo& xi, dinfo& di, double* phi, pinfo& pi, RNG& gen)
{
   tree::npv bnv;
   std::vector<sinfo> sv(di.tlen);
   allsuffhet_ts(t,xi,di,phi,bnv,sv);

   List post_pars;
   mat Phi;
   vec m;
   vec mu_draw;

   for(tree::npv::size_type i=0;i!=bnv.size();i++) {

      // Draw new mu value from MVN.
      post_pars = mvn_post_util_het(pi.mu0, pi.Prec0, sv[i].n0_vec, sv[i].n_vec, sv[i].sy_vec);

      Phi = Rcpp::as<arma::mat>(post_pars["Phi"]);
      m = Rcpp::as<arma::vec>(post_pars["m"]);
      mu_draw = rmvnorm_post(m, Phi);

      // Assign botton node values to new mu draw.
      bnv[i] -> setm(mu_draw);

      // Check for NA result.
      if(sum(bnv[i]->getm() == bnv[i]->getm()) == 0) {
         for(size_t i=0; i<di.n; ++i) Rcout << *(di.x + i*di.p) <<" "; //*(x + p*i+j)
         Rcpp::stop("drmu failed");
      }
   }
}

//--------------------------------------------------
// normal density N(x, mean, variance)
double pn(double x, double m, double v)
{
	double dif = x-m;
	return exp(-.5*dif*dif/v)/sqrt(2*PI*v);
}

//--------------------------------------------------
// draw from discrete distributin given by p, return index
int rdisc(double *p, RNG& gen)
{

	double sum;
	double u = gen.uniform();

    int i=0;
    sum=p[0];
    while(sum<u) {
		i += 1;
		sum += p[i];
    }
    return i;
}

//--------------------------------------------------
//evalute tree tr on grid given by xi and write to os
void grm(tree& tr, xinfo& xi, std::ostream& os)
{
	size_t p = xi.size();
	if(p!=2) {
		cout << "error in grm, p !=2\n";
		return;
	}
	size_t n1 = xi[0].size();
	size_t n2 = xi[1].size();
	tree::tree_cp bp; //pointer to bottom node
	double *x = new double[2];
	for(size_t i=0;i!=n1;i++) {
		for(size_t j=0;j!=n2;j++) {
			x[0] = xi[0][i];
			x[1] = xi[1][j];
			bp = tr.bn(x,xi);
			os << x[0] << " " << x[1] << " " << bp->getm() << " " << bp->nid() << endl;
		}
	}
	delete[] x;
}

//--------------------------------------------------
//does this bottom node n have any variables it can split on.
bool cansplit(tree::tree_p n, xinfo& xi)
{
	int L,U;
	bool v_found = false; //have you found a variable you can split on
	size_t v=0;
	while(!v_found && (v < xi.size())) { //invar: splitvar not found, vars left
		L=0; U = xi[v].size()-1;
		n->rg(v,&L,&U);
		if(U>=L) v_found=true;
		v++;
	}
	return v_found;
}

//--------------------------------------------------
//compute prob of a birth, goodbots will contain all the good bottom nodes
double getpb(tree& t, xinfo& xi, pinfo& pi, tree::npv& goodbots)
{
	double pb;  //prob of birth to be returned
	tree::npv bnv; //all the bottom nodes
	t.getbots(bnv);
	for(size_t i=0;i!=bnv.size();i++)
		if(cansplit(bnv[i],xi)) goodbots.push_back(bnv[i]);
	if(goodbots.size()==0) { //are there any bottom nodes you can split on?
		pb=0.0;
	} else {
		if(t.treesize()==1) pb=1.0; //is there just one node?
		else pb=pi.pb;
	}
	return pb;
}

//--------------------------------------------------
//find variables n can split on, put their indices in goodvars
void getgoodvars(tree::tree_p n, xinfo& xi,  std::vector<size_t>& goodvars)
{
	int L,U;
	for(size_t v=0;v!=xi.size();v++) {//try each variable
		L=0; U = xi[v].size()-1;
		n->rg(v,&L,&U);
		if(U>=L) goodvars.push_back(v);
	}
}

//--------------------------------------------------
//get prob a node grows, 0 if no good vars, else alpha/(1+d)^beta
double pgrow(tree::tree_p n, xinfo& xi, pinfo& pi)
{
	if(cansplit(n,xi)) {
		return pi.alpha/pow(1.0+n->depth(),pi.beta);
	} else {
		return 0.0;
	}
}

//--------------------------------------------------

//get counts for all bottom nodes
std::vector<int> counts(tree& x, xinfo& xi, dinfo& di, tree::npv& bnv)
{
  tree::tree_cp tbn; //the pointer to the bottom node for the current observations
	size_t ni;         //the  index into vector of the current bottom node
	double *xx;        //current x
	double y;          //current y

	bnv.clear();
	x.getbots(bnv);

	typedef tree::npv::size_type bvsz;
//	bvsz nb = bnv.size();

  std::vector<int> cts(bnv.size(), 0);

	std::map<tree::tree_cp,size_t> bnmap;
	for(bvsz i=0;i!=bnv.size();i++) bnmap[bnv[i]]=i;

	for(size_t i=0;i<di.n;i++) {
		xx = di.x + i*di.p;
		y=di.y[i];

		tbn = x.bn(xx,xi);
		ni = bnmap[tbn];

    cts[ni] += 1;
	}
  return(cts);
}

void update_counts(int i, std::vector<int>& cts, tree& x, xinfo& xi,
                   dinfo& di,
                   tree::npv& bnv, //vector of pointers to bottom nodes
                   int sign)
{
  tree::tree_cp tbn; //the pointer to the bottom node for the current observations
  size_t ni;         //the  index into vector of the current bottom node
	double *xx;        //current x
	double y;          //current y

	typedef tree::npv::size_type bvsz;
//	bvsz nb = bnv.size();

	std::map<tree::tree_cp,size_t> bnmap;
	for(bvsz ii=0;ii!=bnv.size();ii++) bnmap[bnv[ii]]=ii; // bnmap[pointer] gives linear index

	xx = di.x + i*di.p;
	y=di.y[i];

	tbn = x.bn(xx,xi);
	ni = bnmap[tbn];

  cts[ni] += sign;
}

void update_counts(int i, std::vector<int>& cts, tree& x, xinfo& xi,
                   dinfo& di,
                   std::map<tree::tree_cp,size_t>& bnmap,
                   int sign)
{
  tree::tree_cp tbn; //the pointer to the bottom node for the current observations
  size_t ni;         //the  index into vector of the current bottom node
  double *xx;        //current x
	double y;          //current y
  /*
	typedef tree::npv::size_type bvsz;
	bvsz nb = bnv.size();

	std::map<tree::tree_cp,size_t> bnmap;
	for(bvsz ii=0;ii!=bnv.size();ii++) bnmap[bnv[ii]]=ii; // bnmap[pointer] gives linear index
	*/
	xx = di.x + i*di.p;
	y=di.y[i];

	tbn = x.bn(xx,xi);
	ni = bnmap[tbn];

  cts[ni] += sign;
}


void update_counts(int i, std::vector<int>& cts, tree& x, xinfo& xi,
                   dinfo& di,
                   std::map<tree::tree_cp,size_t>& bnmap,
                   int sign,
                   tree::tree_cp &tbn
                   )
{
  //tree::tree_cp tbn; //the pointer to the bottom node for the current observations
  size_t ni;         //the  index into vector of the current bottom node
  double *xx;        //current x
  double y;          //current y
  /*
	typedef tree::npv::size_type bvsz;
	bvsz nb = bnv.size();

	std::map<tree::tree_cp,size_t> bnmap;
	for(bvsz ii=0;ii!=bnv.size();ii++) bnmap[bnv[ii]]=ii; // bnmap[pointer] gives linear index
	*/
	xx = di.x + i*di.p;
	y=di.y[i];

	tbn = x.bn(xx,xi);
	ni = bnmap[tbn];

  cts[ni] += sign;
}

bool min_leaf(int minct, std::vector<tree>& t, xinfo& xi, dinfo& di) {
  bool good = true;
  tree::npv bnv;
  std::vector<int> cts;
  int m = 0;
  for (size_t tt=0; tt<t.size(); ++tt) {
    cts = counts(t[tt], xi, di, bnv);
    m = std::min(m, *std::min_element(cts.begin(), cts.end()));
    if(m<minct) {
      good = false;
      break;
    }
  }
  return good;
}

//--------------------------------------------------
//fit for multiple data points, not by reference.
void fit(tree& t, xinfo& xi, dinfo& di, vec& fv)
{
   double *xx;
   tree::tree_cp bn;
   fv.resize(di.n);

   arma::uvec id; //idx of current obs t value.

   for(size_t i=0;i<di.n;i++) {
      xx = di.x + i*di.p;
      bn = t.bn(xx,xi);

      // Find index of mu (time point) corresponding to each obs.  Use this mu.
      fv[i] = bn->getm(di.t[i], di.tref);
   }
}

//--------------------------------------------------
//fit for multiple data points, by reference.
void fit(tree& t, xinfo& xi, dinfo& di, double* fv)
{
   double *xx;
   tree::tree_cp bn;

   arma::uvec id; //idx of current obs t value.


   for(size_t i=0;i<di.n;i++) {
      xx = di.x + i*di.p;
      bn = t.bn(xx,xi);

      // Find index of mu (time point) corresponding to each obs.  Use this mu.
      fv[i] = bn->getm(di.t[i], di.tref);
   }
}

//--------------------------------------------------
//partition
void partition(tree& t, xinfo& xi, dinfo& di, std::vector<size_t>& pv)
{
	double *xx;
	tree::tree_cp bn;
	pv.resize(di.n);
	for(size_t i=0;i<di.n;i++) {
		xx = di.x + i*di.p;
		bn = t.bn(xx,xi);
		pv[i] = bn->nid();
	}
}

//--------------------------------------------------
//write cutpoint information to screen
void prxi(xinfo& xi)
{
	cout << "xinfo: \n";
	for(size_t v=0;v!=xi.size();v++) {
		cout << "v: " << v << endl;
		for(size_t j=0;j!=xi[v].size();j++) cout << "j,xi[v][j]: " << j << ", " << xi[v][j] << endl;
	}
	cout << "\n\n";
}

//--------------------------------------------------
//make xinfo = cutpoints
void makexinfo(size_t p, size_t n, double *x, xinfo& xi, size_t nc)
{
	double xinc;

	//compute min and max for each x
	std::vector<double> minx(p,INFINITY);
	std::vector<double> maxx(p,-INFINITY);
	double xx;
	for(size_t i=0;i<p;i++) {
		for(size_t j=0;j<n;j++) {
			xx = *(x+p*j+i);
			if(xx < minx[i]) minx[i]=xx;
			if(xx > maxx[i]) maxx[i]=xx;
		}
	}
	//make grid of nc cutpoints between min and max for each x.
	xi.resize(p);
	for(size_t i=0;i<p;i++) {
		xinc = (maxx[i]-minx[i])/(nc+1.0);
		xi[i].resize(nc);
		for(size_t j=0;j<nc;j++) xi[i][j] = minx[i] + (j+1)*xinc;
	}
}
// get min/max needed to make cutpoints
void makeminmax(size_t p, size_t n, double *x, std::vector<double> &minx, std::vector<double> &maxx)
{
	double xx;

	for(size_t i=0;i<p;i++) {
		for(size_t j=0;j<n;j++) {
			xx = *(x+p*j+i);
			if(xx < minx[i]) minx[i]=xx;
			if(xx > maxx[i]) maxx[i]=xx;
		}
	}
}
//make xinfo = cutpoints give the minx and maxx vectors
void makexinfominmax(size_t p, xinfo& xi, size_t nc, std::vector<double> &minx, std::vector<double> &maxx)
{
	double xinc;
	//make grid of nc cutpoints between min and max for each x.
	xi.resize(p);
	for(size_t i=0;i<p;i++) {
		xinc = (maxx[i]-minx[i])/(nc+1.0);
		xi[i].resize(nc);
		for(size_t j=0;j<nc;j++) xi[i][j] = minx[i] + (j+1)*xinc;
	}
}

// Check if a vector is sorted.  For checking z and zpred for causal funbart.
bool is_sort(arma::vec x) {
     int n=x.n_elem;
     for (int i=0; i<n-1; ++i)
         if (x[i] < x[i+1]) return false;
     return true;
}

// //--------------------------------------------------
// //log of the integrated likelihood
// double lil(double n, double sy, double sy2, double sigma, double tau)
// {
//    double yb,yb2,S,sig2,d;
//    double sum, rv;
//
//    yb = sy/n;
//    yb2 = yb*yb;
//    S = sy2 - (n*yb2);
//    sig2 = sigma*sigma;
//    d = n*tau*tau + sig2;
//    sum = S/sig2 + (n*yb2)/d;
//    rv = -(n*LTPI/2.0) - (n-1)*log(sigma) -log(d)/2.0;
//    rv = rv -sum/2.0;
//    return rv;
// }
//
// double lilhet(double n, double sy, double sy2, double sigma, double tau)
// {
//    double d = 1/(tau*tau) + n;// n is \sum phi_i for het
//
//    double out = -log(tau) - 0.5*log(d);
//    out += 0.5*sy*sy/d - 0.5*sy2;
//    return out;
// }
//
// double lilprec(double n, double sy, double sy2, double sigma, double tau)
// {
//    double LN2 = 0.693147180559945309417232121458;
//    //gamma prior
//    //double rv = tau*log(tau) - (n*LTPI/2.0) - lgamma(tau);
//    //rv += lgamma(tau + 0.5*n) - (tau+0.5*n)*log(tau + 0.5*sy2);
//
//    //mixture prior
//    double rv = tau*log(tau) - (n*LTPI/2.0) - lgamma(tau) - LN2;
//    double loga = lgamma(tau + 0.5*n) - (tau+0.5*n)*log(tau + 0.5*sy2); //nc for gamma
//    double logb = gig_norm(0.5*n-tau, 2.0*tau, sy2);
//    rv += logsumexp(loga, logb);
//    return rv;
// }


// //get sufficients stats for all bottom nodes (sy, sy2)
// void allsuff(tree& x, xinfo& xi, dinfo& di, tree::npv& bnv, std::vector<sinfo>& sv)
// {
//    // Bottom nodes are written to bnv.
//    // Suff stats for each bottom node are written to elements (each of class sinfo) of sv.
//    // Initialize data structures
//    tree::tree_cp tbn; //the pointer to the bottom node for the current observations.  tree_cp bc not modifying tree directly.
//    size_t ni;         //the  index into vector of the current bottom node
//    double *xx;        //current x
//    double y;          //current y
//
//    bnv.clear();      // Clear the bnv variable if any value is already saved there.
//    x.getbots(bnv);   // Save bottom nodes for x to bnv variable.
//
//    // Not sure what this part here is doing.
//    typedef tree::npv::size_type bvsz;  // Is a better C way to set type.  (tree::npv::size_type) will resolve to an integer,
//    // or long int, etc.  We don't have to know that ahead of time by using this notation.
//    bvsz nb = bnv.size();   // Initialize new var nb of type bvsz for number of bottom nodes, then...
//    sv.resize(nb);          // Re-sizing suff stat vector to have same size as bottom nodes.
//
//    // bnmap is a tuple (lookups, like in Python).  Want to index by bottom nodes.
//    std::map<tree::tree_cp,size_t> bnmap;
//    for(bvsz i=0;i!=bnv.size();i++) bnmap[bnv[i]]=i;  // bnv[i]
//    //map looks like
//    // bottom node 1 ------ 1
//    // bottom node 2 ------ 2
//
//    // Sum the y values (sy) and the y^2 values (sy2) for each node and store in sv.
//    // Loop through each observation.  Push each obs x down the tree and find its bottom node,
//    // then index into the suff stat for the bottom node corresponding to that obs.
//    for(size_t i=0;i<di.n;i++) {
//       xx = di.x + i*di.p;  //Index value: di.x is pointer to first element of n*p data vector.  Iterates through each element.
//       y=di.y[i];           // Resolves to r.
//
//       tbn = x.bn(xx,xi); // Find bottom node for this observation.
//       ni = bnmap[tbn];   // Map bottom node to integer index
//
//       // Update the sufficient stats for the bottom node to which that obs belongs.
//       ++(sv[ni].n);
//       sv[ni].sy += y;
//       sv[ni].sy2 += y*y;
//    }
// }


// //--------------------------------------------------
// //get sufficient stats for children (v,c) of node nx in tree x
// void getsuff(tree& x, tree::tree_cp nx, size_t v, size_t c, xinfo& xi, dinfo& di, sinfo& sl, sinfo& sr)
// {
//    double *xx;//current x
//    double y;  //current y
//    sl.n=0;sl.sy=0.0;sl.sy2=0.0;
//    sr.n=0;sr.sy=0.0;sr.sy2=0.0;
//
//    for(size_t i=0;i<di.n;i++) {
//       xx = di.x + i*di.p;
//       if(nx==x.bn(xx,xi)) { //does the bottom node = xx's bottom node
//          y = di.y[i];
//          if(xx[v] < xi[v][c]) {
//             sl.n++;
//             sl.sy += y;
//             sl.sy2 += y*y;
//          } else {
//             sr.n++;
//             sr.sy += y;
//             sr.sy2 += y*y;
//          }
//       }
//    }
// }

// //--------------------------------------------------
// //get sufficient stats for pair of bottom children nl(left) and nr(right) in tree x
// void getsuff(tree& x, tree::tree_cp nl, tree::tree_cp nr, xinfo& xi, dinfo& di, sinfo& sl, sinfo& sr)
// {
//    double *xx;//current x
//    double y;  //current y
//    sl.n=0;sl.sy=0.0;sl.sy2=0.0;
//    sr.n=0;sr.sy=0.0;sr.sy2=0.0;
//
//    for(size_t i=0;i<di.n;i++) {
//       xx = di.x + i*di.p;
//       tree::tree_cp bn = x.bn(xx,xi);
//       if(bn==nl) {
//          y = di.y[i];
//          sl.n++;
//          sl.sy += y;
//          sl.sy2 += y*y;
//       }
//       if(bn==nr) {
//          y = di.y[i];
//          sr.n++;
//          sr.sy += y;
//          sr.sy2 += y*y;
//       }
//    }
// }

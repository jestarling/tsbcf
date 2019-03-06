#include <iostream>

#include "info.h"
#include "tree.h"
#include "bd.h"
#include "funs.h"

using std::cout;
using std::endl;

/*
notation: (as in old code): going from state x to state y (eg, incoming tree is x).

note: rather than have x and making a tree y
we just figure out what we need from x, the drawn bottom node,the drawn (v,c).
note sure what the right thing to do is.
Could make y (using a birth) and figure stuff out from y.
That is how the old code works.
*/

double bd(tree& x, xinfo& xi, dinfo& di, pinfo& pi, RNG& gen)
{
   tree::npv goodbots;                    //nodes we could birth at (split on)
   double PBx = getpb(x,xi,pi,goodbots);  //prob of a birth at x

   //Rcpp::Rcout << "PBx: " << PBx << endl;

   // If statement for selecting birth or death proposal.
   if(gen.uniform() < PBx) {
      //--------------------------------------------------
      // BIRTH PROPOSAL
      //--------------------------------------------------

      //--------------------------------------------------
      //draw proposal

      //draw bottom node, choose node index ni from list in goodbots
      size_t ni = floor(gen.uniform()*goodbots.size());
      tree::tree_p nx = goodbots[ni]; //the bottom node we might birth at

      //draw v,  the variable
      std::vector<size_t> goodvars; //variables nx can split on
      getgoodvars(nx,xi,goodvars);
      size_t vi = floor(gen.uniform()*goodvars.size()); //index of chosen split variable
      size_t v = goodvars[vi];

      //draw c, the cutpoint
      int L,U;
      L=0; U = xi[v].size()-1;
      nx->rg(v,&L,&U);
      size_t c = L + floor(gen.uniform()*(U-L+1)); //U-L+1 is number of available split points

      //--------------------------------------------------
      //compute things needed for metropolis ratio

      double Pbotx = 1.0/goodbots.size(); //proposal prob of choosing nx
      size_t dnx = nx->depth();
      double PGnx = pi.alpha/pow(1.0 + dnx,pi.beta); //prior prob of growing at nx

      double PGly, PGry; //prior probs of growing at new children (l and r) of proposal
      if(goodvars.size()>1) { //know there are variables we could split l and r on
         PGly = pi.alpha/pow(1.0 + dnx+1.0,pi.beta); //depth of new nodes would be one more
         PGry = PGly;
      } else { //only had v to work with, if it is exhausted at either child need PG=0
         if((int)(c-1)<L) { //v exhausted in new left child l, new upper limit would be c-1
            PGly = 0.0;
         } else {
            PGly = pi.alpha/pow(1.0 + dnx+1.0,pi.beta);
         }
         if(U < (int)(c+1)) { //v exhausted in new right child r, new lower limit would be c+1
            PGry = 0.0;
         } else {
            PGry = pi.alpha/pow(1.0 + dnx+1.0,pi.beta);
         }
      }

      double PDy; //prob of proposing death at y
      if(goodbots.size()>1) { //can birth at y because splittable nodes left
         PDy = 1.0 - pi.pb;
      } else { //nx was the only node you could split on
         if((PGry==0) && (PGly==0)) { //cannot birth at y
            PDy=1.0;
         } else { //y can birth at either l or r
            PDy = 1.0 - pi.pb;
         }
      }

      double Pnogy; //death prob of choosing the nog node at y
      size_t nnogs = x.nnogs();
      tree::tree_cp nxp = nx->getp();
      if(nxp==0) { //no parent, nx is the top and only node
         Pnogy=1.0;
      } else {
         //if(nxp->ntype() == 'n') { //if parent is a nog, number of nogs same at x and y
         if(nxp->isnog()) { //if parent is a nog, number of nogs same at x and y
            Pnogy = 1.0/nnogs;
         } else { //if parent is not a nog, y has one more nog.
           Pnogy = 1.0/(nnogs+1.0);
         }
      }

      //--------------------------------------------------
      //compute sufficient statistics
      sinfo sl,sr; //sl for left from nx and sr for right from nx (using rule (v,c))
      getsuff_ts(x,nx,v,c,xi,di,sl,sr, di.tlen);

      //--------------------------------------------------
      //compute alpha

      double alpha=0.0,alpha1=0.0,alpha2=0.0;
      double lill=0.0,lilr=0.0,lilt=0.0;

      //Rcpp::Rcout << "birth alpha sl.n: " << sl.n << ", sr.n: " << sr.n << endl;
      //if((sl.n>=5) && (sr.n>=5)) { //cludge?
      if((sl.n>=2) && (sr.n>=2)) { //cludge?

         lill = lil_ts(sl.n_vec, sl.sy_vec, sl.sy2, pi.sigma, pi.mu0, pi.Prec0);                        // USED lil(...)
         lilr = lil_ts(sr.n_vec, sr.sy_vec, sr.sy2, pi.sigma, pi.mu0, pi.Prec0);                        // USED lil(...)
         lilt = lil_ts(sl.n_vec + sr.n_vec, sl.sy_vec + sr.sy_vec, sl.sy2 + sr.sy2, pi.sigma, pi.mu0, pi.Prec0);  // USED lil(...)

         alpha1 = (PGnx*(1.0-PGly)*(1.0-PGry)*PDy*Pnogy)/((1.0-PGnx)*PBx*Pbotx);
         alpha2 = alpha1*exp(lill+lilr-lilt);
         alpha = std::min(1.0,alpha2);

       } else {
         alpha=0.0;
      }

      /*
      cout << "sigma, tau: " << pi.sigma << ", " << pi.tau << endl;
      cout << "birth prop: node, v, c: " << nx->nid() << ", " << v << ", " << c << "," << xi[v][c] << endl;
      cout << "L,U: " << L << "," << U << endl;
      cout << "PBx, PGnx, PGly, PGry, PDy, Pnogy,Pbotx:" <<
         PBx << "," << PGnx << "," << PGly << "," << PGry << "," << PDy <<
         ", " << Pnogy << "," << Pbotx << endl;
      cout << "left ss: " << sl.n << ", " << sl.sy << ", " << sl.sy2 << endl;
      cout << "right ss: " << sr.n << ", " << sr.sy << ", " << sr.sy2 << endl;
      cout << "lill, lilr, lilt: " << lill << ", " << lilr << ", " << lilt << endl;
      cout << "alphas: " << alpha1 << ", " << alpha2 << ", " << alpha << endl;
      */

      //--------------------------------------------------
      //finally ready to try metrop
      //--------------------------------------------------

      vec mul, mur;
      mul = zeros(di.tlen);
      mur = zeros(di.tlen);
      if(gen.uniform() < alpha) {

         //--------------------------------------------------
         // do birth:
         // Set mul and mur to zero vectors, since we will immediately
         // fill them in the MCMC by using drmu.  This saves computation cost.

			x.birth(nx->nid(),v,c,mul,mur);
         return alpha;
      } else {
         return alpha+10;
      }
   } else {
      //--------------------------------------------------
      // DEATH PROPOSAL
      //--------------------------------------------------

       //--------------------------------------------------
      //draw proposal

      //draw nog node, any nog node is a possibility
      tree::npv nognds; //nog nodes
      x.getnogs(nognds);
      size_t ni = floor(gen.uniform()*nognds.size());
      tree::tree_p nx = nognds[ni]; //the nog node we might kill children at

      //--------------------------------------------------
      //compute things needed for metropolis ratio

      double PGny; //prob the nog node grows
      size_t dny = nx->depth();
      PGny = pi.alpha/pow(1.0+dny,pi.beta);

      //better way to code these two?
      double PGlx = pgrow(nx->getl(),xi,pi);
      double PGrx = pgrow(nx->getr(),xi,pi);

      double PBy;  //prob of birth move at y
      //if(nx->ntype()=='t') { //is the nog node nx the top node
      if(!(nx->p)) { //is the nog node nx the top node
         PBy = 1.0;
      } else {
         PBy = pi.pb;
      }

      double Pboty;  //prob of choosing the nog as bot to split on when y
      int ngood = goodbots.size();
      if(cansplit(nx->getl(),xi)) --ngood; //if can split at left child, lose this one
      if(cansplit(nx->getr(),xi)) --ngood; //if can split at right child, lose this one
      ++ngood;  //know you can split at nx
      Pboty=1.0/ngood;

      double PDx = 1.0-PBx; //prob of a death step at x
      double Pnogx = 1.0/nognds.size();

      //--------------------------------------------------
      //compute sufficient statistics
      sinfo sl(di.tlen),sr(di.tlen); //sl for left from nx and sr for right from nx (using rule (v,c))

      // Test if problem is getl and getr.
      getsuff_ts(x,nx->getl(),nx->getr(),xi,di,sl,sr,di.tlen);

      //--------------------------------------------------
      //compute alpha

      double lill = lil_ts(sl.n_vec, sl.sy_vec, sl.sy2, pi.sigma, pi.mu0, pi.Prec0);
      double lilr = lil_ts(sr.n_vec, sr.sy_vec, sr.sy2, pi.sigma, pi.mu0, pi.Prec0);
      double lilt = lil_ts(sl.n_vec + sr.n_vec, sl.sy_vec + sr.sy_vec, sl.sy2 + sr.sy2, pi.sigma, pi.mu0, pi.Prec0);

      double alpha1 = ((1.0-PGny)*PBy*Pboty)/(PGny*(1.0-PGlx)*(1.0-PGrx)*PDx*Pnogx);
      double alpha2 = alpha1*exp(lilt - lill - lilr);
      double alpha = std::min(1.0,alpha2);

      /*
      cout << "death prop: " << nx->nid() << endl;
      cout << "nognds.size(), ni, nx: " << nognds.size() << ", " << ni << ", " << nx << endl;
      cout << "depth of nog node: " << dny << endl;
      cout << "PGny: " << PGny << endl;
      cout << "PGlx: " << PGlx << endl;
      cout << "PGrx: " << PGrx << endl;
      cout << "PBy: " << PBy << endl;
      cout << "Pboty: " << Pboty << endl;
      cout << "PDx: " << PDx << endl;
      cout << "Pnogx: " << Pnogx << endl;
      cout << "left ss: " << sl.n << ", " << sl.sy << ", " << sl.sy2 << endl;
      cout << "right ss: " << sr.n << ", " << sr.sy << ", " << sr.sy2 << endl;
      cout << "lill, lilr, lilt: " << lill << ", " << lilr << ", " << lilt << endl;
      cout << "sigma: " << pi.sigma << endl;
      cout << "tau: " << pi.tau << endl;
      cout << "alphas: " << alpha1 << ", " << alpha2 << ", " << alpha << endl;
      */

      //--------------------------------------------------
      //finally ready to try metrop

      // All values wrong, but ok bc get reset immediately in the draw mu.
      if(gen.uniform()<alpha) {
         //draw mu for nog (which will be bot)
         vec mu;
         mu = zeros(di.tlen);

         //do death
			x.death(nx->nid(),mu);
         return -alpha;
      } else {
         return -alpha-10;
      }
   }
}

// Het variances.
double bdhet(tree& x, xinfo& xi, dinfo& di, double* phi, pinfo& pi, RNG& gen)
{
   tree::npv goodbots;                    //nodes we could birth at (split on)
   double PBx = getpb(x,xi,pi,goodbots);  //prob of a birth at x
   //Rcpp::Rcout << "PBx: " << PBx << endl;

   // If statement for selecting birth or death proposal.
   if(gen.uniform() < PBx) {

      //--------------------------------------------------
      // BIRTH PROPOSAL
      //--------------------------------------------------

       //--------------------------------------------------
      //draw proposal

      //draw bottom node, choose node index ni from list in goodbots
      size_t ni = floor(gen.uniform()*goodbots.size());
      tree::tree_p nx = goodbots[ni]; //the bottom node we might birth at

       //draw v,  the variable
      std::vector<size_t> goodvars; //variables nx can split on
      getgoodvars(nx,xi,goodvars);
      size_t vi = floor(gen.uniform()*goodvars.size()); //index of chosen split variable
      size_t v = goodvars[vi];

      //draw c, the cutpoint
      int L,U;
      L=0; U = xi[v].size()-1;
      nx->rg(v,&L,&U);
      size_t c = L + floor(gen.uniform()*(U-L+1)); //U-L+1 is number of available split points

      //--------------------------------------------------
      //compute things needed for metropolis ratio

      double Pbotx = 1.0/goodbots.size(); //proposal prob of choosing nx
      size_t dnx = nx->depth();
      double PGnx = pi.alpha/pow(1.0 + dnx,pi.beta); //prior prob of growing at nx

      double PGly, PGry; //prior probs of growing at new children (l and r) of proposal
      if(goodvars.size()>1) { //know there are variables we could split l and r on
         PGly = pi.alpha/pow(1.0 + dnx+1.0,pi.beta); //depth of new nodes would be one more
         PGry = PGly;
      } else { //only had v to work with, if it is exhausted at either child need PG=0
         if((int)(c-1)<L) { //v exhausted in new left child l, new upper limit would be c-1
            PGly = 0.0;
         } else {
            PGly = pi.alpha/pow(1.0 + dnx+1.0,pi.beta);
         }
         if(U < (int)(c+1)) { //v exhausted in new right child r, new lower limit would be c+1
            PGry = 0.0;
         } else {
            PGry = pi.alpha/pow(1.0 + dnx+1.0,pi.beta);
         }
      }

      double PDy; //prob of proposing death at y
      if(goodbots.size()>1) { //can birth at y because splittable nodes left
         PDy = 1.0 - pi.pb;
      } else { //nx was the only node you could split on
         if((PGry==0) && (PGly==0)) { //cannot birth at y
            PDy=1.0;
         } else { //y can birth at either l or r
            PDy = 1.0 - pi.pb;
         }
      }

      double Pnogy; //death prob of choosing the nog node at y
      size_t nnogs = x.nnogs();
      tree::tree_cp nxp = nx->getp();
      if(nxp==0) { //no parent, nx is the top and only node
         Pnogy=1.0;
      } else {
         //if(nxp->ntype() == 'n') { //if parent is a nog, number of nogs same at x and y
         if(nxp->isnog()) { //if parent is a nog, number of nogs same at x and y
            Pnogy = 1.0/nnogs;
         } else { //if parent is not a nog, y has one more nog.
            Pnogy = 1.0/(nnogs+1.0);
         }
      }

      //--------------------------------------------------
      //compute sufficient statistics
      sinfo sl,sr; //sl for left from nx and sr for right from nx (using rule (v,c))
      getsuffhet_ts(x,nx,v,c,xi,di,phi,sl,sr,di.tlen);

      //--------------------------------------------------
      //compute alpha

      double alpha=0.0,alpha1=0.0,alpha2=0.0;
      double lill=0.0,lilr=0.0,lilt=0.0;

      //if((sl.n>=5) && (sr.n>=5)) {
      if((abs(sl.n)>=2) && (abs(sr.n)>=2)) { //cludge?

         lill = lilhet_ts(sl.n0, sl.n, sl.n_vec, sl.sy_vec, sl.sy2, pi.mu0, pi.Prec0);
         lilr = lilhet_ts(sr.n0, sr.n, sr.n_vec, sr.sy_vec, sr.sy2, pi.mu0, pi.Prec0);
         lilt = lilhet_ts(sl.n0+sr.n0, sl.n+sr.n, sl.n_vec+sr.n_vec,
                          sl.sy_vec+sr.sy_vec,
                          sl.sy2+sr.sy2, pi.mu0, pi.Prec0);

         alpha1 = (PGnx*(1.0-PGly)*(1.0-PGry)*PDy*Pnogy)/((1.0-PGnx)*PBx*Pbotx);
         alpha2 = alpha1*exp(lill+lilr-lilt);
         alpha = std::min(1.0,alpha2);
      } else {
         alpha=0.0;
      }

      /*
      cout << "sigma, tau: " << pi.sigma << ", " << pi.tau << endl;
      cout << "birth prop: node, v, c: " << nx->nid() << ", " << v << ", " << c << "," << xi[v][c] << endl;
      cout << "L,U: " << L << "," << U << endl;
      cout << "PBx, PGnx, PGly, PGry, PDy, Pnogy,Pbotx:" <<
         PBx << "," << PGnx << "," << PGly << "," << PGry << "," << PDy <<
            ", " << Pnogy << "," << Pbotx << endl;
      cout << "left ss: " << sl.n << ", " << sl.sy << ", " << sl.sy2 << endl;
      cout << "right ss: " << sr.n << ", " << sr.sy << ", " << sr.sy2 << endl;
      cout << "lill, lilr, lilt: " << lill << ", " << lilr << ", " << lilt << endl;
      cout << "alphas: " << alpha1 << ", " << alpha2 << ", " << alpha << endl;
      */

      //--------------------------------------------------
      //finally ready to try metrop
      //--------------------------------------------------

      vec mul, mur;
      mul = zeros(di.tlen);
      mur = zeros(di.tlen);
      if(gen.uniform() < alpha) {

         //--------------------------------------------------
         // do birth:
         // Set mul and mur to zero vectors, since we will immediately
         // fill them in the MCMC by using drmu.  This saves computation cost.

         x.birth(nx->nid(),v,c,mul,mur);
         return alpha;
      } else {
         return alpha+10;
      }
   } else {
      //--------------------------------------------------
      // DEATH PROPOSAL
      //--------------------------------------------------

      //--------------------------------------------------
      //draw proposal

      //draw nog node, any nog node is a possibility
      tree::npv nognds; //nog nodes
      x.getnogs(nognds);
      size_t ni = floor(gen.uniform()*nognds.size());
      tree::tree_p nx = nognds[ni]; //the nog node we might kill children at

      //--------------------------------------------------
      //compute things needed for metropolis ratio

      double PGny; //prob the nog node grows
      size_t dny = nx->depth();
      PGny = pi.alpha/pow(1.0+dny,pi.beta);

      //better way to code these two?
      double PGlx = pgrow(nx->getl(),xi,pi);
      double PGrx = pgrow(nx->getr(),xi,pi);

      double PBy;  //prob of birth move at y
      //if(nx->ntype()=='t') { //is the nog node nx the top node
      if(!(nx->p)) { //is the nog node nx the top node
         PBy = 1.0;
      } else {
         PBy = pi.pb;
      }

      double Pboty;  //prob of choosing the nog as bot to split on when y
      int ngood = goodbots.size();
      if(cansplit(nx->getl(),xi)) --ngood; //if can split at left child, lose this one
      if(cansplit(nx->getr(),xi)) --ngood; //if can split at right child, lose this one
      ++ngood;  //know you can split at nx
      Pboty=1.0/ngood;

      double PDx = 1.0-PBx; //prob of a death step at x
      double Pnogx = 1.0/nognds.size();

      //--------------------------------------------------
      //compute sufficient statistics
      sinfo sl(di.tlen),sr(di.tlen); //sl for left from nx and sr for right from nx (using rule (v,c))

      // Test if problem is getl and getr.
      getsuffhet_ts(x,nx->getl(),nx->getr(),xi,di,phi,sl,sr,di.tlen);

      //--------------------------------------------------
      //compute alpha

      double lill = lilhet_ts(sl.n0, sl.n, sl.n_vec, sl.sy_vec, sl.sy2, pi.mu0, pi.Prec0);
      double lilr = lilhet_ts(sr.n0, sr.n, sr.n_vec, sr.sy_vec, sr.sy2, pi.mu0, pi.Prec0);
      double lilt = lilhet_ts(sl.n0+sr.n0, sl.n+sr.n, sl.n_vec+sr.n_vec,
                              sl.sy_vec+sr.sy_vec,
                              sl.sy2+sr.sy2, pi.mu0, pi.Prec0);

      double alpha1 = ((1.0-PGny)*PBy*Pboty)/(PGny*(1.0-PGlx)*(1.0-PGrx)*PDx*Pnogx);
      double alpha2 = alpha1*exp(lilt - lill - lilr);
      double alpha = std::min(1.0,alpha2);

      /*
      cout << "death prop: " << nx->nid() << endl;
      cout << "nognds.size(), ni, nx: " << nognds.size() << ", " << ni << ", " << nx << endl;
      cout << "depth of nog node: " << dny << endl;
      cout << "PGny: " << PGny << endl;
      cout << "PGlx: " << PGlx << endl;
      cout << "PGrx: " << PGrx << endl;
      cout << "PBy: " << PBy << endl;
      cout << "Pboty: " << Pboty << endl;
      cout << "PDx: " << PDx << endl;
      cout << "Pnogx: " << Pnogx << endl;
      cout << "left ss: " << sl.n << ", " << sl.sy << ", " << sl.sy2 << endl;
      cout << "right ss: " << sr.n << ", " << sr.sy << ", " << sr.sy2 << endl;
      cout << "lill, lilr, lilt: " << lill << ", " << lilr << ", " << lilt << endl;
      cout << "sigma: " << pi.sigma << endl;
      cout << "tau: " << pi.tau << endl;
      cout << "alphas: " << alpha1 << ", " << alpha2 << ", " << alpha << endl;
      */

      //--------------------------------------------------
      //finally ready to try metrop

      // All values wrong, but ok bc get reset immediately in the draw mu.
      if(gen.uniform()<alpha) {
         //draw mu for nog (which will be bot)
         vec mu;
         mu = zeros(di.tlen);

         //do death
         x.death(nx->nid(),mu);

         //return true;
         return -alpha;
      } else {
         return -alpha-10;
      }
   }
}

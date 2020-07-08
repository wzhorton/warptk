#include <RcppArmadillo.h>
using namespace Rcpp;

// [[Rcpp::depends(RcppArmadillo)]]

// ##############################################################
// ##############################################################
//              UTILITY FUNCTIONS
// ##############################################################
// ##############################################################

arma::vec qform(arma::mat y, arma::vec mu, arma::mat covprec, bool is_prec){
  arma::mat y_shift;
  y_shift.copy_size(y.t());
  for(int i = 0; i <  y.n_cols; i++){
    y_shift.row(i) = (y.col(i) - mu).t();
  }
  if (is_prec) {
    return sum((y_shift * covprec) % y_shift, 1);
  } else {
    return sum(solve(covprec,y_shift.t()).t() % y_shift, 1);
  }
}

arma::vec qformI(arma::mat y, arma::vec mu){
  arma::mat y_shift;
  y_shift.copy_size(y.t());
  for(int i = 0; i <  y.n_cols; i++){
    y_shift.row(i) = (y.col(i) - mu).t();
  }
  return sum(y_shift % y_shift, 1);
}

arma::vec rep4k(arma::vec u){
  int n = u.n_elem;
  double tmp3;
  double ui;
  arma::vec out(n);
  for(int i=0; i<n; i++){
    ui = u(i);
    if(ui > 1.0){
      out(i) = 0.0;
    } else if(ui >= 0.0){
      tmp3 = 1.0 - ui;
      out(i) = tmp3*tmp3*tmp3;
    } else {
      out(i) = 0.0;
    }
  }
  return out;
}

arma::vec rep3k(arma::vec u){
  int n = u.n_elem;
  double ui;
  arma::vec out(n);
  for(int i=0; i<n; i++){
    ui = u(i);
    if(ui > 2.0){
      out(i) = 0.0;
    } else if(ui > 1.0){
      out(i) = 2.0 - 3.0*ui + 3.0/2.0*ui*ui - 1.0/4.0*ui*ui*ui;
    } else if(ui > 0.0){
      out(i) = 3.0*ui - 9.0/2.0*ui*ui + 7.0/4.0*ui*ui*ui;
    } else {
      out(i) = 0.0;
    }
  }
  return out;
}

arma::vec rep2k(arma::vec u){
  int n = u.n_elem;
  double ui;
  arma::vec out(n);
  for(int i=0; i<n; i++){
    ui = u(i);
    if(ui > 3.0){
      out(i) = 0.0;
    } else if(ui > 2.0){
      out(i) = 9.0/2.0 - 9.0/2.0*ui + 3.0/2.0*ui*ui - 1.0/6.0*ui*ui*ui;
    } else if(ui > 1.0){
      out(i) = -3.0/2.0 + 9.0/2.0*ui - 3.0*ui*ui + 7.0/12.0*ui*ui*ui;
    } else if(ui > 0.0){
      out(i) = 3.0/2.0*ui*ui - 11.0/12.0*ui*ui*ui;
    } else {
      out(i) = 0.0;
    }
  }
  return out;
}

arma::vec evenk(arma::vec u){
  int n = u.n_elem;
  double ui;
  arma::vec out(n);
  for(int i=0; i<n; i++){
    ui = u(i);
    if(ui > 4.0){
      out(i) = 0.0;
    } else if(ui > 3.0){
      out(i) = -1.0/6.0 * (ui-4)*(ui-4)*(ui-4);
    } else if(ui > 2.0){
      out(i) = -22.0/3.0 + 10.0*ui - 4.0*ui*ui + 1.0/2.0*ui*ui*ui;
    } else if(ui > 1.0){
      out(i) = 2.0/3.0 - 2.0*ui + 2.0*ui*ui - 1.0/2.0*ui*ui*ui;
    } else if(ui > 0.0){
      out(i) = 1.0/6.0*ui*ui*ui;
    } else {
      out(i) = 0.0;
    }
  }
  return out;
}

// [[Rcpp::export(".bs_even_C")]]
arma::mat bs_even(arma::vec time, int nk){
  arma::vec u;
  double dnk = nk - 1.0;
  u = dnk / time.max() * (time - time.min());
  int nbasis = nk + 2;
  int ni = nk - 2;
  int neven_basis = ni - 2;
  arma::mat out(time.n_elem, nbasis);
  out.col(0) = rep4k(u);
  out.col(1) = rep3k(u);
  out.col(2) = rep2k(u);
  arma::vec revu = 3 + neven_basis - u;
  out.col(nbasis - 1) = rep4k(revu);
  out.col(nbasis - 2) = rep3k(revu);
  out.col(nbasis - 3) = rep2k(revu);
  for(int i=0; i < neven_basis; i++){
    out.col(3 + i) = evenk(u - i);
  }
  return out;
}

arma::vec rmnorm(arma::vec mu, arma::mat covprec, bool is_prec) {//slower but more stable
  arma::vec z = rnorm(mu.n_elem);
  arma::vec evals;
  arma::mat evecs;
  if (is_prec) {
    return mu + solve(chol(covprec), z);
  } else {
    //return mu + chol(covprec).t() * z;
    arma::vec evals;
    arma::mat evecs;
    arma::eig_sym(evals, evecs, covprec);
    arma::mat Dhalf = diagmat(sqrt(arma::clamp(evals, 0.0, evals.max())));
    //arma::mat Dhalf = diagmat(sqrt(evals));
    return mu + evecs*Dhalf*evecs.t()*z;
  }
}

// [[Rcpp::export(".det_sympd_C")]]
double det_sympd (arma::mat x, bool Log = false) {
  arma::mat cholx = chol(x);
  arma::vec y = log(cholx.diag());
  if( Log ) {
    return 2 * sum(y);
  } else {
    return exp(2 * sum(y));
  }
}

// [[Rcpp::export(".dist_C")]]
arma::mat dist(arma::vec x, arma::vec y){
  arma::mat dists(x.n_elem, y.n_elem);
  for(int i=0; i<y.n_elem; i++){
    dists.col(i) = abs(x - y(i));
  }
  return dists;
}

bool is_increasing(arma::vec x){
  for(int i = 1; i < x.n_elem; i++){
    if(x(i) <= x(i-1)){
      return false;
    }
  }
  return true;
}

double bound(double x, double mn, double mx){
  return std::min(mx, std::max(x, mn));
}

// [[Rcpp::export(".monotonize_C")]]
void curve_monotonize(arma::vec x, arma::vec &y){

  int n = x.n_elem;
  int i, ii, eloop;
  double mn;
  y = x;
  for(i=0; i<(n-1); i++){
    if(y(i) > y(i+1)){
      mn = arma::mean(y.subvec(i, i+1));
      y(i) = mn;
      y(i+1) = mn;
      ii = i;
      eloop = 1;
      while(eloop){
        if(i > 0){
          if(y(ii) < y(ii-1)){
            y.subvec(ii-1, i+1).fill(arma::mean(y.subvec(ii-1, i+1)));
            ii = ii - 1;
            if(ii == 0) eloop = 0;
          } else {
            eloop = 0;
          }
        } else {
          eloop = 0;
        }
      }
    }
  }
  for(i=0; i<(n-1); i++){
    y(i) = bound(y(i), 0, 1);
  }
}

// ##############################################################
// ##############################################################
//              REGISTRATION FUNCTIONS
// ##############################################################
// ##############################################################

// [[Rcpp::export(".two_step_warp_C")]]
List two_step_warp(arma::mat ymat, arma::mat wtime, arma::mat P, int niter, int nburn, int nthin,
                   double a_eps, double b_eps, double a_a, double b_a, double a_c, double b_c,
                   double a_tau, double b_tau){
  int p = P.n_cols;
  int n = ymat.n_cols;
  int m = ymat.n_rows;
  int nrun = niter + nburn;
  int save_index = 0;

  double sig2_eps = 1;
  arma::vec sig2e_chain(niter/nthin);

  double tau2 = 1;
  arma::vec tau2_chain(niter/nthin);

  arma::vec a(n);
  a.zeros();
  arma::mat a_chain(n,niter/nthin);

  arma::vec c(n);
  c.ones();
  arma::mat c_chain(n,niter/nthin);


  arma::vec ac0(2);
  ac0(0) = 0;
  ac0(1) = 1;

  arma::vec beta = rnorm(p);
  arma::mat beta_chain(p,niter/nthin);

  double sig2_a = 1;
  arma::vec sig2a_chain(niter/nthin);

  double sig2_c = 1;
  arma::vec sig2c_chain(niter/nthin);



  arma::cube H(m,p,n);
  for(int i=0; i<n; i++){
    H.slice(i) = bs_even(wtime.col(i),p-2);
  }

  arma::mat mu(m,n);
  arma::cube mu_chain(m,n,niter/nthin);

  double ssq_sig2;
  arma::mat ssq_bV(p,p);
  arma::vec ssq_bm(p);

  arma::mat W(m,2);
  W.ones();

  arma::vec ya(n);
  arma::mat cH(n,p);
  arma::vec ac(2);

  arma::mat VVb(p,p);
  arma::mat VVac;
  arma::mat Viac(2,2);
  Viac.zeros();

  for(int iter=0; iter<nrun; iter++){
    ssq_sig2 = 0;
    ssq_bV.zeros();
    ssq_bm.zeros();
    Viac(0,0) = 1/sig2_a;
    Viac(1,1) = 1/sig2_c;
    for(int i=0; i<n; i++){
      W.col(1) = H.slice(i)*beta;
      VVac = Viac + 1/sig2_eps*W.t()*W;
      ac = rmnorm(solve(VVac, Viac*ac0 + 1/sig2_eps*W.t()*ymat.col(i)), VVac, true);
      a(i) = ac(0);
      c(i) = ac(1);
    }
    a -= mean(a);
    c -= mean(c) - 1;
    for(int i=0; i<n; i++){
      cH = c(i)*H.slice(i);
      mu.col(i) = a(i) + cH*beta;
      ya = ymat.col(i) - a(i);

      ssq_sig2 += dot(ymat.col(i) - mu.col(i), ymat.col(i) - mu.col(i));
      ssq_bV += cH.t()*cH;
      ssq_bm += cH.t()*ya;
    }

    sig2_eps = 1/rgamma(1, a_eps + 0.5*n*m, 1/(b_eps + 0.5*ssq_sig2))(0);
    VVb = 1/sig2_eps*ssq_bV + 1/tau2*P;
    beta = rmnorm(solve(VVb, 1/sig2_eps*ssq_bm), VVb, true);

    sig2_a = 1/rgamma(1, a_a + n/2, 1/(b_a + 0.5*dot(a,a)))(0);
    sig2_c = 1/rgamma(1, a_c + n/2, 1/(b_c + 0.5*dot(c-1,c-1)))(0);
    tau2 = 1/rgamma(1, a_tau + p/2, 1/(b_tau + 0.5*arma::as_scalar(beta.t()*P*beta)))(0);

    if(iter >= nburn && (iter % nthin) == 0){
      sig2e_chain(save_index) = sig2_eps;
      sig2a_chain(save_index) = sig2_a;
      sig2c_chain(save_index) = sig2_c;
      tau2_chain(save_index) = tau2;
      beta_chain.col(save_index) = beta;
      mu_chain.slice(save_index) = mu;
      a_chain.col(save_index) = a;
      c_chain.col(save_index) = c;
      save_index++;
    }
  }

  List chains = List::create(a_chain, c_chain, beta_chain, sig2e_chain, sig2a_chain, sig2c_chain, tau2_chain, wtime, mu_chain);
  return chains;
}





// [[Rcpp::export(".bhcr_warp_C")]]
List bhcr_warp(arma::mat ymat, arma::vec time, arma::mat P, arma::mat Q, arma::vec U, int niter, int nburn, int nthin,
                   double a_eps, double b_eps, double a_a, double b_a, double a_c, double b_c,
                   double a_tau, double b_tau, double a_lam, double b_lam){
  int p = P.n_cols;
  int q = Q.n_cols;
  int n = ymat.n_cols;
  int m = ymat.n_rows;
  int nrun = niter + nburn;
  int save_index = 0;


  double sig2_eps = 1;
  arma::vec sig2e_chain(niter/nthin);

  double tau2 = 1;
  arma::vec tau2_chain(niter/nthin);

  double lam2 = 1;
  arma::vec lam2_chain(niter/nthin);

  double sig2_a = 1;
  arma::vec sig2a_chain(niter/nthin);

  double sig2_c = 1;
  arma::vec sig2c_chain(niter/nthin);



  arma::vec a(n);
  a.randu();
  arma::mat a_chain(n,niter/nthin);

  arma::vec c(n);
  c.randu();
  arma::mat c_chain(n,niter/nthin);

  arma::vec beta = rnorm(p);
  arma::mat beta_chain(p,niter/nthin);

  arma::mat phi(q,n);
  arma::cube phi_chain(q, n, niter/nthin);
  for(int i=0; i<n; i++){
    phi.col(i) = U;
  }


  arma::mat wtime(m,n);
  arma::cube wtime_chain(m,n,niter/nthin);
  for(int i=0; i<n;i++){
    wtime.col(i) = time;
  }

  arma::cube H(m,p,n);
  for(int i=0; i<n; i++){
    H.slice(i) = bs_even(wtime.col(i),p-2);
  }
  arma::mat H_q = bs_even(time,q-2);

  arma::mat mu(m,n);
  arma::cube mu_chain(m,n,niter/nthin);
  for(int i=0; i<n;i++){
    mu.col(i) = a(i) + c(i)*H.slice(i)*beta;
  }

  double ssq_sig2;
  arma::mat ssq_bV(p,p);
  arma::vec ssq_bm(p);
  double ssq_lam2;

  arma::vec ac0(2);
  ac0(0) = 0;
  ac0(1) = 1;
  double tune = 0.005;

  arma::mat W(m,2);
  W.ones();

  arma::vec cand_phi = U;
  double cand_lpost;
  double cand_elem;
  arma::vec cand_mu;
  arma::mat cand_H;
  arma::vec curr_lpost(n);
  for(int i=0; i<n; i++){
    curr_lpost(i) = -0.5/sig2_eps*dot(ymat.col(i) - mu.col(i), ymat.col(i) - mu.col(i)) -
        0.5/lam2*arma::as_scalar((phi.col(i)-U).t()*Q*(phi.col(i)-U));
  }

  arma::vec ya(n);
  arma::mat cH(n,p);
  arma::vec ac(2);

  arma::mat VVb(p,p);
  arma::mat VVac;
  arma::mat Viac(2,2);
  Viac.zeros();

  for(int iter=0; iter<nrun; iter++){
    ssq_sig2 = 0;
    ssq_lam2 = 0;
    ssq_bV.zeros();
    ssq_bm.zeros();
    Viac(0,0) = 1/sig2_a;
    Viac(1,1) = 1/sig2_c;
    for(int i=0; i<n; i++){
      //wtime.col(i) = time;
      //H.slice(i) = bs_even(wtime.col(i),p-2);
      cand_phi = phi.col(i);
      curr_lpost(i) = -0.5/sig2_eps*dot(ymat.col(i) - mu.col(i), ymat.col(i) - mu.col(i)) -
        0.5/lam2*arma::as_scalar((phi.col(i)-U).t()*Q*(phi.col(i)-U));

      for(int j=1; j<(q-1); j++){
        cand_elem = cand_phi(j) + runif(1,-tune,tune)(0);
        //cand_elem = runif(1,cand_phi(j-1),cand_phi(j+1))(0);
        if(cand_elem > cand_phi(j-1) && cand_elem < cand_phi(j+1)){
          cand_phi(j) = cand_elem;
          cand_H = bs_even(H_q*cand_phi,p-2);
          cand_mu = a(i) + c(i)*cand_H*beta;
          cand_lpost = -0.5/sig2_eps*dot(ymat.col(i) - cand_mu, ymat.col(i) - cand_mu) -
                0.5/lam2*arma::as_scalar((cand_phi-U).t()*Q*(cand_phi-U));
          if(log(runif(1)(0)) < cand_lpost - curr_lpost(i)){
            curr_lpost(i) = cand_lpost;
          } else {
            cand_phi(j) = phi.col(i)(j);
          }
        }
      }
      phi.col(i) = cand_phi;
      wtime.col(i) = H_q*phi.col(i);
      H.slice(i) = bs_even(wtime.col(i),p-2);

      W.col(1) = H.slice(i)*beta;
      VVac = Viac + 1/sig2_eps*W.t()*W;
      ac = rmnorm(solve(VVac, Viac*ac0 + 1/sig2_eps*W.t()*ymat.col(i)), VVac, true);
      a(i) = ac(0);
      c(i) = ac(1);
    }
    a -= mean(a);
    c -= mean(c) - 1;

    for(int i=0; i<n; i++){
      cH = c(i)*H.slice(i);
      mu.col(i) = a(i) + cH*beta;
      ya = ymat.col(i) - a(i);

      ssq_sig2 += dot(ymat.col(i) - mu.col(i), ymat.col(i) - mu.col(i));
      ssq_bV += cH.t()*cH;
      ssq_bm += cH.t()*ya;
      ssq_lam2 += arma::as_scalar((phi.col(i) - U).t()*Q*(phi.col(i) - U));
    }

    sig2_eps = 1/rgamma(1, a_eps + 0.5*n*m, 1/(b_eps + 0.5*ssq_sig2))(0);
    VVb = 1/sig2_eps*ssq_bV + 1/tau2*P;
    beta = rmnorm(solve(VVb, 1/sig2_eps*ssq_bm), VVb, true);

    sig2_a = 1/rgamma(1, a_a + n/2, 1/(b_a + 0.5*dot(a,a)))(0);
    sig2_c = 1/rgamma(1, a_c + n/2, 1/(b_c + 0.5*dot(c-1,c-1)))(0);
    tau2 = 1/rgamma(1, a_tau + p/2, 1/(b_tau + 0.5*arma::as_scalar(beta.t()*P*beta)))(0);
    lam2 = 1/rgamma(1, a_lam + q*n/2, 1/(b_lam + 0.5*ssq_lam2))(0);

    if(iter >= nburn && (iter % nthin) == 0){
      sig2e_chain(save_index) = sig2_eps;
      sig2a_chain(save_index) = sig2_a;
      sig2c_chain(save_index) = sig2_c;
      tau2_chain(save_index) = tau2;
      lam2_chain(save_index) = lam2;
      beta_chain.col(save_index) = beta;
      mu_chain.slice(save_index) = mu;
      wtime_chain.slice(save_index) = wtime;
      phi_chain.slice(save_index) = phi;
      a_chain.col(save_index) = a;
      c_chain.col(save_index) = c;
      save_index++;
    }
  }

  List chains = List::create(a_chain, c_chain, beta_chain, sig2e_chain, sig2a_chain, sig2c_chain, tau2_chain, lam2_chain, phi_chain, wtime_chain, mu_chain);
  return chains;
}







// [[Rcpp::export(".template_warp_C")]]
List template_warp(arma::mat ymat, arma::vec time, arma::mat wtime_init, arma::mat lmk_time, arma::vec ref_time, arma::mat P,
                   int niter, int nburn, int nthin,
                   double a_eps, double b_eps, double a_a, double b_a, double a_c, double b_c,
                   double a_tau, double b_tau, double a_lam, double b_lam){
  int p = P.n_cols;
  int n = ymat.n_cols;
  int m = ymat.n_rows;
  int l = ref_time.n_elem;
  int nrun = niter + nburn;
  int save_index = 0;


  double sig2_eps = 1;
  arma::vec sig2e_chain(niter/nthin);

  double tau2 = 1;
  arma::vec tau2_chain(niter/nthin);

  double lam2 = 1;
  arma::vec lam2_chain(niter/nthin);

  double sig2_a = 1;
  arma::vec sig2a_chain(niter/nthin);

  double sig2_c = 1;
  arma::vec sig2c_chain(niter/nthin);

  double alpha = 5;
  arma::vec alpha_chain(niter/nthin);

  double eta = .005;
  arma::vec eta_chain(niter/nthin);



  arma::vec a(n);
  a.randu();
  arma::mat a_chain(n,niter/nthin);

  arma::vec c(n);
  c.randu();
  arma::mat c_chain(n,niter/nthin);

  arma::vec beta = rnorm(p);
  arma::mat beta_chain(p,niter/nthin);

  arma::mat D2time = arma::square(dist(time, time));
  arma::cube D2feat(l,l,n);
  arma::cube D2timefeat(m,l,n);
  for(int i=0; i<n; i++){
    D2timefeat.slice(i) = arma::square(dist(time, lmk_time.col(i)));
    D2feat.slice(i) = arma::square(dist(lmk_time.col(i), lmk_time.col(i)));
  }
  arma::mat M_time = arma::exp(-alpha/2*D2time);
  M_time.diag() += eta; // only employ nugget on observed data, not these interpolation points (only in generation!)
  //M_time(0,0) = 1; // might be more correct, but this is more stable
  //M_time(m-1,m-1) = 1;
  arma::cube M_feat(l,l,n);
  arma::cube M_timefeat(m,l,n);
  for(int i=0; i<n; i++){
    M_timefeat.slice(i) = arma::exp(-alpha/2*D2timefeat.slice(i));
    M_feat.slice(i) = arma::exp(-alpha/2*D2feat.slice(i));
    M_feat.slice(i).diag() += eta;
    //M_feat.slice(i)(0,0) = 1;
    //M_feat.slice(i)(l-1,l-1) = 1;
  }

  double lam2_gen = .00005;
  arma::mat M_time_gen = arma::exp(-50/2*D2time);
  arma::cube M_feat_gen(l,l,n);
  arma::cube M_timefeat_gen(m,l,n);
  for(int i=0; i<n; i++){
    M_timefeat_gen.slice(i) = arma::exp(-50/2*D2timefeat.slice(i));
    M_feat_gen.slice(i) = arma::exp(-50/2*D2feat.slice(i));
    //M_feat_gen.slice(i).diag() += .0005;
    M_feat_gen.slice(i)(0,0) = 1; // may not be stable
    M_feat_gen.slice(i)(l-1,l-1) = 1;
  }

  arma::mat wtime = wtime_init;
  arma::cube wtime_chain(m,n,niter/nthin);
  /*for(int i=0; i<n;i++){
    wtime.col(i) = time + M_timefeat_gen.slice(i)*solve(M_feat_gen.slice(i), ref_time - lmk_time.col(i));
  }*/

  arma::cube H(m,p,n);
  for(int i=0; i<n; i++){
    H.slice(i) = bs_even(wtime.col(i),p-2);
  }

  arma::mat mu(m,n);
  arma::cube mu_chain(m,n,niter/nthin);
  for(int i=0; i<n;i++){
    mu.col(i) = a(i) + c(i)*H.slice(i)*beta;
  }

  double ssq_sig2;
  arma::mat ssq_bV(p,p);
  arma::vec ssq_bm(p);
  double cand_ssq_lam2;
  arma::vec ssq_wtime = qform(wtime, time, M_time, false); //COULD CAUSE ISSUES TAKING THE INVERSE. SHOULD IT HAVE A NUGGET AND JUST NOT WHEN GENERATING DATA? yes.
  double ssq_lam2 = sum(ssq_wtime);
  arma::vec cand_ssq_wtime(n);

  arma::vec cand_ae(2);
  arma::mat tune_ae(2,2);
  tune_ae(0,0) = 2;
  tune_ae(1,1) = .0001;
  tune_ae(1,0) = 0;
  tune_ae(0,1) = 0;
  arma::vec ae(2);
  ae(0) = alpha;
  ae(1) = eta;
  double curr_lpost_ae;
  double cand_lpost_ae;
  arma::mat candM;

  arma::vec ac0(2);
  ac0(0) = 0;
  ac0(1) = 1;


  arma::mat W(m,2);
  W.ones();

  arma::mat cH(n,p);
  arma::vec ac(2);

  arma::mat VVb(p,p);
  arma::mat VVac;
  arma::mat Viac(2,2);
  Viac.zeros();

  arma::vec cand_wtime;
  double curr_lpost_wtime;
  double cand_lpost_wtime;
  arma::mat candH;
  arma::vec cand_mu;

  double ldet_val;
  double ldet_sign;
  arma::mat M_time_inv = arma::inv_sympd(M_time);
  arma::vec mon_warp(n);


  for(int iter=0; iter<nrun; iter++){

    ////////// Update Lambda^2 //////////

    // Update Step
    lam2 = 1/rgamma(1, a_lam + m*n/2, 1/(b_lam + 0.5*ssq_lam2))(0);

    ////////// Update Alpha, Eta //////////

    //Update Step
    log_det(ldet_val, ldet_sign, lam2*M_time); //return List::create(save_index);
    curr_lpost_ae = -0.5*n*ldet_val - 0.5/lam2*ssq_lam2;
    cand_ae = abs(rmnorm(ae, tune_ae, false)); //clean up conditioning, why are negatives getting in
    //cand_ae(0) = rnorm(1, alpha, tune_ae(0,0))(0);
    //cand_ae(1) = rnorm(1, log(eta), tune_ae(1,1))(0); // log transform of eta, maybe revert?
    cand_ae(1) = eta; // fix to see if alpha becomes stable.
    if(cand_ae(0) > 0 && cand_ae(0) < m && cand_ae(1) > 0 && cand_ae(1) < 5){ // change m to 50 or something
      candM = arma::exp(-cand_ae(0)/2*D2time);
      candM.diag() += cand_ae(1);// exp from transform, see initialization, removed
      //candM(0,0) = 1; // supurfluous and causes instability, removed.
      //candM(m-1,m-1) = 1;
      cand_ssq_wtime = qform(wtime, time, candM, false);
      cand_ssq_lam2 = sum(cand_ssq_wtime);
      log_det(ldet_val, ldet_sign, lam2*candM);
      cand_lpost_ae =  -0.5*n*ldet_val - 0.5/lam2*cand_ssq_lam2;// + cand_ae(1) - log(5.0);//jacobian, removed
      //Rcpp::Rcout << "ldet: " << ldet_val << ", dsym: " << det_sympd(lam2*candM, true) << std::endl;
      //Rcpp::Rcout << "cand-curr: " << cand_lpost_ae - curr_lpost_ae << std::endl;
      if(log(runif(1)(0)) < cand_lpost_ae - curr_lpost_ae){
        // Primary Updates
        ae = cand_ae;
        alpha = ae(0);
        eta = ae(1);    //Rcpp::Rcout << "eta: " << eta << std::endl; // transform, removed
        // Peripheral Updates
        ssq_lam2 = cand_ssq_lam2;
        ssq_wtime = cand_ssq_wtime;
        M_time = candM;
        M_time_inv = arma::inv_sympd(M_time);
        for(int i=0; i<n; i++){
          M_timefeat.slice(i) = arma::exp(-alpha/2*D2timefeat.slice(i));
          M_feat.slice(i) = arma::exp(-alpha/2*D2feat.slice(i));
          M_feat.slice(i).diag() += eta;
          //M_feat.slice(i)(0,0) = 1;
          //M_feat.slice(i)(m-1,m-1) = 1; // supurfluous with stability concerns, removed
        }
      }
    }


    ////////// Update warp_time //////////



    //Rcpp::Rcout << "Flag 1" << std::endl;o
    for(int i=0; i<n; i++){ //Could speed up by cutting out the middle loop. Remake qform function (or make another)
      curr_lpost_wtime = -.5/sig2_eps*dot(ymat.col(i)-mu.col(i), ymat.col(i)-mu.col(i)) - .5/lam2*ssq_wtime(i);
      //return List::create(M_timefeat_gen, M_feat_gen, M_time_gen);///////////skdjfhskdjfhskdjhsdf
      cand_wtime = rmnorm(wtime.col(i),// + M_timefeat_gen.slice(i)*solve(M_feat_gen.slice(i), ref_time - lmk_time.col(i)), //0 out after shift.
                          lam2_gen*(M_time_gen - M_timefeat_gen.slice(i)*solve(M_feat_gen.slice(i), M_timefeat_gen.slice(i).t())), false);
      if(i==1){
        //Rcpp::Rcout << (M_time_gen - M_timefeat_gen.slice(i)*solve(M_feat_gen.slice(i), M_timefeat_gen.slice(i).t()))(34,34) << std::endl;
      }
      //for(int j=0; j<l; j++){
      //  cand_wtime(lmk_inds(j,i)) = ref_time(j);
      //}
      /// THIS ISN"T CENTERED ABOUT THE PREVIOUS CURVES AND IT SHOULDN"T HAVE THE NUGGET> (but liks shouls) START FROM D2 SCRATCH
      candH = bs_even(cand_wtime,p-2); // NEED TO PROPOSE VALUES NOT USING ETA AND ALPHA IN M!!
      cand_mu = a(i) + c(i)*candH*beta;
      cand_ssq_wtime(i) = qform(cand_wtime, time, M_time_inv, true)(0);
      cand_lpost_wtime = -.5/sig2_eps*dot(ymat.col(i)-cand_mu, ymat.col(i)-cand_mu) - .5/lam2*cand_ssq_wtime(i);
      //Rcpp::Rcout << "cand-curr: " << cand_lpost_wtime-curr_lpost_wtime <<std::endl;

      if(iter < nburn){
         if(is_increasing(cand_wtime) && log(runif(1)(0)) < cand_lpost_wtime - curr_lpost_wtime){
           // Primary Updates
           wtime.col(i) = cand_wtime;
           H.slice(i) = candH;

           // Peripheral Updates
           ssq_wtime(i) = cand_ssq_wtime(i);
         }
       } else {
         if(log(runif(1)(0)) < cand_lpost_wtime - curr_lpost_wtime){
           // Primary Updates
           curve_monotonize(cand_wtime, mon_warp);
           wtime.col(i) = mon_warp;
           H.slice(i) = bs_even(mon_warp,p-2);

           // Peripheral Updates
           ssq_wtime(i) = qform(mon_warp, time, M_time_inv, true)(0);
         }
       }
    }
    ssq_lam2 = sum(ssq_wtime);


    ////////// Update a_i and c_i for all i //////////

    // Peripherals
    Viac(0,0) = 1/sig2_a;
    Viac(1,1) = 1/sig2_c;

    // Update Step
    for(int i=0; i<n; i++){
      W.col(1) = H.slice(i)*beta;
      VVac = Viac + 1/sig2_eps*W.t()*W;
      ac = rmnorm(solve(VVac, Viac*ac0 + 1/sig2_eps*W.t()*ymat.col(i)), VVac, true);
      a(i) = ac(0);
      c(i) = ac(1);
    }
    a -= mean(a);
    c -= mean(c) - 1;


   ////////// Update Beta //////////

    //Peripherals
    ssq_bV.zeros();
    ssq_bm.zeros();
    for(int i=0; i<n; i++){
      cH = c(i)*H.slice(i);
      mu.col(i) = a(i) + cH*beta;
      ssq_bV += cH.t()*cH;
      ssq_bm += cH.t()*(ymat.col(i) - a(i));
    }
    ssq_sig2 = arma::accu(arma::square(ymat - mu));
    VVb = 1/sig2_eps*ssq_bV + 1/tau2*P;

    //Update Step
    beta = rmnorm(solve(VVb, 1/sig2_eps*ssq_bm), VVb, true);

    ////////// Update Variance Parameters //////////
    sig2_eps = 1/rgamma(1, a_eps + 0.5*n*m, 1/(b_eps + 0.5*ssq_sig2))(0);
    sig2_a = 1/rgamma(1, a_a + n/2, 1/(b_a + 0.5*dot(a,a)))(0);
    sig2_c = 1/rgamma(1, a_c + n/2, 1/(b_c + 0.5*dot(c-1,c-1)))(0);
    tau2 = 1/rgamma(1, a_tau + p/2, 1/(b_tau + 0.5*arma::as_scalar(beta.t()*P*beta)))(0);


    if(iter >= nburn && (iter % nthin) == 0){
      sig2e_chain(save_index) = sig2_eps;
      sig2a_chain(save_index) = sig2_a;
      sig2c_chain(save_index) = sig2_c;
      tau2_chain(save_index) = tau2;
      lam2_chain(save_index) = lam2;
      alpha_chain(save_index) = alpha;
      eta_chain(save_index) = eta;
      beta_chain.col(save_index) = beta;
      mu_chain.slice(save_index) = mu;
      wtime_chain.slice(save_index) = wtime;
      a_chain.col(save_index) = a;
      c_chain.col(save_index) = c;
      save_index++;
    }
  }

  List chains = List::create(a_chain, c_chain, beta_chain, sig2e_chain, sig2a_chain,
                             sig2c_chain, tau2_chain, lam2_chain, alpha_chain, eta_chain, wtime_chain, mu_chain);
  return chains;
}




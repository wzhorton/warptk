#include <RcppArmadillo.h>
using namespace Rcpp;

// [[Rcpp::depends(RcppArmadillo)]]

// ##############################################################
// ##############################################################
//              UTILITY FUNCTIONS
// ##############################################################
// ##############################################################

// [[Rcpp::export(".dmat_C")]]
arma::mat dmat(arma::vec x) {
  int n = x.n_elem;
  arma::mat outmat(n, n);
  double dist;

  for(int i = 0; i < n; i++){
    for(int j = 0; j <= i; j++){
      dist = abs(x(i) - x(j));
      outmat(i,j) = dist;
      outmat(j,i) = dist;
    }
  }
  return outmat;
}

arma::vec qform(arma::mat &y, arma::rowvec &mu, arma::mat &covprec, bool is_prec){
  int nrows = y.n_rows;
  arma::mat y_shift;
  y_shift.copy_size(y);
  for(int i = 0; i < nrows; i++){
    y_shift.row(i) = y.row(i) - mu;
  }
  if (is_prec) {
    return sum((y_shift * covprec) % y_shift, 1);
  } else {
    return sum((y_shift * inv_sympd(covprec)) % y_shift, 1);
  }
}

arma::vec rep4k(arma::vec &u){
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

arma::vec rep3k(arma::vec &u){
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

arma::vec rep2k(arma::vec &u){
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

arma::vec rmnorm(arma::vec mu, arma::mat covprec, bool is_prec) {
  arma::vec z = rnorm(mu.n_elem);
  if (is_prec) {
    return mu + solve(chol(covprec), z);
  } else {
    return mu + chol(covprec).t() * z;
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

    if(iter >= nburn & (iter % nthin) == 0){
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

  List chains = List::create(a_chain, c_chain, beta_chain, sig2e_chain, sig2a_chain, sig2c_chain, tau2_chain, mu_chain);
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
        if(cand_elem > cand_phi(j-1) & cand_elem < cand_phi(j+1)){
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

    if(iter >= nburn & (iter % nthin) == 0){
      sig2e_chain(save_index) = sig2_eps;
      sig2a_chain(save_index) = sig2_a;
      sig2c_chain(save_index) = sig2_c;
      tau2_chain(save_index) = tau2;
      lam2_chain(save_index) = lam2;
      beta_chain.col(save_index) = beta;
      mu_chain.slice(save_index) = mu;
      phi_chain.slice(save_index) = phi;
      a_chain.col(save_index) = a;
      c_chain.col(save_index) = c;
      save_index++;
    }
  }

  List chains = List::create(a_chain, c_chain, beta_chain, sig2e_chain, sig2a_chain, sig2c_chain, tau2_chain, lam2_chain, phi_chain, mu_chain);
  return chains;
}




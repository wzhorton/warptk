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
List two_step_warp(arma::mat ymat, arma::mat wtime, arma::mat P, int niter, int nburn,
                   double a_eps, double b_eps, double a_a, double b_a,
                   double a_c, double b_c, double a_tau, double b_tau){
  int p = P.n_cols;
  int n = ymat.n_cols;
  int m = ymat.n_rows;
  int nrun = niter + nburn;
  int save_index = -nburn-1;

  double sig2_eps = 1;
  arma::vec sig2e_chain(niter);

  double tau2 = 1;
  arma::vec tau2_chain(niter);

  arma::vec a(n);
  a.zeros();
  arma::mat a_chain(n,niter);

  arma::vec c(n);
  c.ones();
  arma::mat c_chain(n,niter);

  arma::vec beta = rnorm(p);
  arma::mat beta_chain(p,niter);

  double sig2_a = 1;
  arma::vec sig2a_chain(niter);

  double sig2_c = 1;
  arma::vec sig2c_chain(niter);

  //arma::mat wtime(m,n);
  //for(int i=0; i<n;i++){
  //  wtime.col(i) = time;
  //}

  arma::cube H(m,p,n);
  for(int i=0; i<n; i++){
    H.slice(i) = bs_even(wtime.col(i),p-2);
  }

  arma::mat mu(m,n);
  for(int i=0; i<n; i++){
    mu.col(i) = a(i) + c(i)*H.slice(i)*beta;
  }
  arma::cube mu_chain(m,n,niter);

  double ssq_sig2;
  arma::mat ssq_bV(p,p);
  arma::vec ssq_bm(p);

  arma::vec HB(n);
  arma::vec ya(n);
  arma::mat cH(n,p);

  arma::mat VV(p,p);

  for(int iter=0; iter<nrun; iter++){
    ssq_sig2 = 0;
    ssq_bV.zeros();
    ssq_bm.zeros();
    for(int i=0; i<n; i++){
      //wtime.col(i) = time;
      //H.slice(i) = bs_even(wtime.col(i),p-2);
      HB = H.slice(i)*beta;
      cH = c(i)*H.slice(i);

      mu.col(i) = a(i) + c(i)*H.slice(i)*beta;

      a(i) = rnorm(1,(arma::sum(ymat.col(i) - c(i)*HB)/sig2_eps)/(m/sig2_eps + 1/sig2_a), 1/(m/sig2_eps + 1/sig2_a))(0);
      ya = ymat.col(i) - a(i);
      c(i) = rnorm(1,(arma::dot(ya, HB)/sig2_eps+1/sig2_c)/(dot(HB,HB)/sig2_eps + 1/sig2_c), 1/(dot(HB,HB)/sig2_eps + 1/sig2_c))(0);

      ssq_sig2 += dot(ya - c(i)*HB, ya - c(i)*HB);
      ssq_bV += cH.t()*cH;
      ssq_bm += cH.t()*ya;
    }

    sig2_eps = 1/rgamma(1, a_eps + 0.5*n*m, b_eps + 0.5*ssq_sig2)(0);
    VV = 1/sig2_eps*ssq_bV + 1/tau2*P;
    beta = rmnorm(solve(VV, 1/sig2_eps*ssq_bm), VV, true);

    sig2_a = 1/rgamma(1, a_a + n/2, b_a + 0.5*dot(a,a))(0);
    sig2_c = 1/rgamma(1, a_c + n/2, b_c + 0.5*dot(c-1,c-1))(0);
    tau2 = 1/rgamma(1, a_tau + p/2, b_tau + 0.5*arma::as_scalar(beta.t()*P*beta))(0);

    save_index++;
    if(save_index >= 0){
      sig2e_chain(save_index) = sig2_eps;
      sig2a_chain(save_index) = sig2_a;
      sig2c_chain(save_index) = sig2_c;
      tau2_chain(save_index) = tau2;
      beta_chain.col(save_index) = beta;
      mu_chain.slice(save_index) = mu;
      for(int i=0; i<n; i++){
        a_chain.col(save_index) = a;
        c_chain.col(save_index) = c;
      }
    }
  }

  List chains = List::create(a_chain, c_chain, beta_chain, sig2e_chain, sig2a_chain, sig2c_chain, tau2_chain, mu_chain);
  return chains;
}


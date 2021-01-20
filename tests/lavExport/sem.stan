functions{

  vector[] sem_mean(vector[] alpha, real[,,] B, real[,,] gamma, int[] g, int k, int Ng, int gamind, real[,] meanx){
    matrix[k,k] iden;
    vector[k] evlv[Ng];

    iden = diag_matrix(rep_vector(1.0, k));

    for(j in 1:Ng){
      if(gamind == 1){
        evlv[j] = inverse(iden - to_matrix(B[,,j])) * (alpha[j] + to_matrix(gamma[,,j]) * to_vector(meanx[,j]));

      } else {
        evlv[j] = inverse(iden - to_matrix(B[,,j])) * alpha[j];
      }
    }

    return evlv;
  }

  real sem_lv_lpdf(matrix x, real[,,] alpha, real[,,] B, real[,,] psi, real[,,] gamma, int gamind, real[,] meanx, int[] g, int k, int N, int Ng, int diagpsi, int fullbeta, int nlv, int[] lvind, int nlvno0){
    real ldetcomp[Ng];
    matrix[k,k] iden;
    vector[k] alpha2[Ng];
    vector[k] psivecinv[Ng];
    matrix[k,k] psimatinv[Ng];
    matrix[k,k] psimat[Ng];
    matrix[k,k] siginv[Ng];
    vector[k] xvec;
    vector[k] evlv[Ng];
    int idx[(k-nlv+nlvno0)];
    real xvectm;
    real ldetsum;
    int nov;
    int nidx;

    nov = k - nlv;
    nidx = nov + nlvno0;

    iden = diag_matrix(rep_vector(1.0, k));

    if(nlvno0 > 0){
      idx[1:nlvno0] = lvind;
    }
    if(nov > 0){
      for(j in 1:nov){
        idx[nlvno0+j] = nlv + j; //nlvno0 + j?
      }
    }

    for(j in 1:Ng){
      alpha2[j] = to_vector(alpha[,1,j]);
    }

    evlv = sem_mean(alpha2, B, gamma, g, k, Ng, gamind, meanx);

    if(diagpsi){
      for(j in 1:Ng){
        for(i in 1:nidx){
          psivecinv[j,idx[i]] = 1/psi[idx[i],idx[i],j];
        }
        psimatinv[j] = diag_matrix(psivecinv[j]);

        siginv[j,1:nidx,1:nidx] = (iden[idx,idx] - to_matrix(B[idx,idx,j])') * psimatinv[j,idx,idx] * (iden[idx,idx] - to_matrix(B[idx,idx,j]));

	if(fullbeta){
	  ldetcomp[j] = log_determinant(iden[idx,idx] - to_matrix(B[idx,idx,j]));
	  ldetcomp[j] = -2 * ldetcomp[j] + sum(log(diagonal(to_matrix(psi[idx,idx,j]))));
	} else {
          ldetcomp[j] = sum(log(diagonal(to_matrix(psi[idx,idx,j]))));
  	}
      }
    } else {
      for(j in 1:Ng){
	psimat[j] = to_matrix(psi[,,j]) + to_matrix(psi[,,j])' - diag_matrix(diagonal(to_matrix(psi[,,j])));

	ldetcomp[j] = log_determinant(psimat[j,idx,idx]);
	if(fullbeta){
	  ldetcomp[j] = ldetcomp[j] - 2 * log_determinant(iden[idx,idx] - to_matrix(B[idx,idx,j]));
	}

	psimatinv[j] = psimat[j];
	psimatinv[j,1:nidx,1:nidx] = inverse_spd(psimat[j,idx,idx]);
        siginv[j,1:nidx,1:nidx] = (iden[idx,idx] - to_matrix(B[idx,idx,j])') * psimatinv[j,1:nidx,1:nidx] * (iden[idx,idx] - to_matrix(B[idx,idx,j]));
      }
    }

    xvectm = 0;
    ldetsum = 0;
    for(i in 1:N){
      xvec = x[i,]';
      xvectm = xvectm + (xvec[idx] - evlv[g[i],idx])' * siginv[g[i],1:nidx,1:nidx] * (xvec[idx] - evlv[g[i],idx]);
      ldetsum = ldetsum + ldetcomp[g[i]];
    }

    return -0.5 * (ldetsum + xvectm);
  }

  matrix fill_lower(matrix x){
    matrix[rows(x),cols(x)] newx;

    newx = x;
    for(i in 1:(rows(x) - 1)){
      for(j in (i+1):rows(x)){
        newx[j,i] = x[i,j];
      }
    }
    return newx;
  }
}

data{
  int N;
  int g[N];
  int lvind[2];
  int etaind[2];
  real sampmean[2,2];
  real meanx[1,2];
  int dummyov[2];
  int dummylv[2];
  vector[2] x[N];
  real lambdaframe[2,2,2];
  real thetaframe[2,2,2];
  real psiframe[2,2,2];
  real betaframe[2,2,2];
  real nuframe[2,1,2];
  real alphaframe[2,1,2];
}

parameters{
  vector<lower=0>[2] psifree;
  vector[2] betafree;
  vector[2] alphafree;
}

transformed parameters{
  real lambda[2,2,2];
  real theta[2,2,2];
  matrix[0,0] thetld[2];
  real psi[2,2,2];
  real beta[2,2,2];
  real nu[2,1,2];
  real alpha[2,1,2];
  real mu[N,2];
  matrix[N,2] eta;

  eta = rep_matrix(0, N, 2);

  lambda = lambdaframe;
  theta = thetaframe;
  psi = psiframe;
  beta = betaframe;
  nu = nuframe;
  alpha = alphaframe;

  beta[1,2,1] = betafree[1];
  psi[1,1,1] = pow(psifree[1],-1);
  psi[2,2,1] = 1.14305737120731;
  alpha[1,1,1] = alphafree[1];
  alpha[2,1,1] = -0.0128301239307355;
  beta[1,2,2] = betafree[2];
  psi[1,1,2] = pow(psifree[2],-1);
  psi[2,2,2] = 1.12735433727174;
  alpha[1,1,2] = alphafree[2];
  alpha[2,1,2] = 0.00190004874826493;

  // mu definitions
  for(i in 1:N) {
    eta[i,1:2] = x[i]';
  }

  for(j in 1:2){
  }

}

model {
  for(i in 1:N) {
  }

  eta ~ sem_lv(alpha, beta, psi, beta, 0, meanx, g, 2, N, 2, 1, 0, 0, etaind, 0);

  // Priors
  target += normal_lpdf(betafree[1] | 0,1);
  target += gamma_lpdf(psifree[1] | 1,.5);
  target += normal_lpdf(alphafree[1] | 0,1000^.5);
  target += normal_lpdf(betafree[2] | 0,1);
  target += gamma_lpdf(psifree[2] | 1,.5);
  target += normal_lpdf(alphafree[2] | 0,1000^.5);
}
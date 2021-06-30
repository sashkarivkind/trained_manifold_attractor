function o=num_expectation_gauss_uni(fun,n_gauss,n_uni, n_tries)
mcgauss = randn(n_tries,n_gauss);
mcuni  = rand(n_tries,n_uni);
o = mean(fun(mcgauss,mcuni),1);



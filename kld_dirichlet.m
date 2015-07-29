function [ d ] = kld_dirichlet( alpha_p, alpha_q )

    tq = sum(alpha_q);
    tp = sum(alpha_p);
    qmp = alpha_q - alpha_p;
    d = gammaln( tq ) - gammaln( tp ) + sum( gammaln(alpha_p) - gammaln(alpha_q) ) + sum( qmp .* psi(alpha_q) ) - psi(tq) * sum(qmp);

end

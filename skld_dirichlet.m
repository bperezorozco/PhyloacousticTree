function [ d ] = skld_dirichlet( alpha_p, alpha_q )
%KLD_DIRICHLET Summary of this function goes here
%   Detailed explanation goes here

d = ( kld_dirichlet( alpha_p, alpha_q ) + kld_dirichlet( alpha_q, alpha_p ) ) / 2;

end




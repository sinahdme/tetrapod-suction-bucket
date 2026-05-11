function [P1]=pycurve2(fi,gamma,d,z,do_b,l_b,n,n_l,A_c,m_c)

Sigma_p_0 = gamma.*z;
nu = (1-sind(fi))/(2-sind(fi));
E_oed = m_c*100;%*sqrt(100*Sigma_p_0);
E_50 = (1-nu-2*nu^2)/(1-nu)* E_oed;
% Effective overburden pressure
P_u = (415*Sigma_p_0/(fi/do_b) + 169)*A_c*tan(pi/n); %kN/m

% Shape coefficients
beta_1 = (1.59*10^(-8)*(E_50*do_b)+0.57);
beta_2 = (-1.18*10^(-5)*(E_50*do_b)+32.2)*(4/n/tan(pi/n))^3;
beta_3 = (3.93*10^(-8)*(E_50*do_b)+0.52);
beta_4 = (-1.31*10^(-5)*(E_50*do_b)+23.5)*(4/n/tan(pi/n))^3;

P = P_u.*(beta_1 * (tanh(beta_2 * (d./do_b)))^(1/3) + beta_3 * (tanh(beta_4 * (d./do_b)))^(1/3))*(l_b / n_l); %kN

P1=P;
end
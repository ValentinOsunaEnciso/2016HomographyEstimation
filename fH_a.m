function f = fH_a(x, X,d,T_noise_squared)
f = [];
% bandera=0;
% while ~bandera
%         f(1,1:d)=round(min_range+(max_range-min_range).*rand(1,d));
%         bandera=validateMSS_homography(X, f(in1,1:d));
% end
[Theta, ~] = estimate_homography(X, x(1,1:d));
E=error_homography(Theta, X, [], []);
CS = (E <= T_noise_squared);
f(1)=-length(find(CS));
% f(2)=get_consensus_set_rank(CS, E, 'MSAC', T_noise_squared, 1, d);
f(2) = T_noise_squared;%length(find(CS));
% f(2) = sum(E);%-get_consensus_set_rank(CS, E, 'MSAC', T_noise_squared, 1, d);
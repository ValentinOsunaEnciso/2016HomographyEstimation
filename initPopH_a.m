% Valentin Osuna-Enciso, CUCEI-UDG, Enero, 2014.%%%%%%%%%%%%%%%%%%%%%%%%%%%
function f = initPopH_a(Np, d, min_range, max_range,X)
% Initialize and Evaluates a population vectors of Real numbers 
% Np - Population size
% M  - Number of objective functions
% d  - Number of decision variables
% min_range, max_range - limits of feasible space
K = 2 + d; f=[];
%% Initialize each individual: %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
for in1=1:Np    
    bandera=0;
    while ~bandera
        f(in1,1:d)=min_range+(max_range-min_range).*rand(1,d);
        bandera=validateMSS_homography(X, round(f(in1,1:d-1)));
    end
    f(in1,d + 1:K) = fH_a(round(f(in1,1:d-1)),X,d-1,f(in1,d));
%     [Theta, ~] = estimate_homography(X, f(in1,1:d));
%     E=error_homography(Theta, X, [], []);
%     CS = (E <= T_noise_squared);
%     f(in1,d + 1:K) = -get_consensus_set_rank(CS, E, 'MSAC', T_noise_squared, 1, d);
%     f(in1,d + 2) = length(find(CS));
end
% -------------------------------------------------------------------------
%
% File : run_optimal_unitary_inversion.m
%
% Discription : 
% Code to obtain the maximal fidelity of transforming 'n' calls of any 
% 'd'-dimensional unitary operations into its inverse map.
% The interpreter cvx is used to conduct SDP.
% Please set the parameters 'd', 'n', 'protocol', and 'isComplex'.
%
% -------------------------------------------------------------------------

clear

% ------------------------------------------------------------------
%                   Start of setting the parameters
% ------------------------------------------------------------------

d = 2;                  % Dimension of the input unitary operation
n = 4;                  % Number of calls
protocol = 2;           % 1 for parallel, 2 for sequential
isComplex = 0;          % Set 0 for real Choi matrix

% ------------------------------------------------------------------
%                   End of setting the parameters
% ------------------------------------------------------------------

tic;

switch protocol
    case 1 
        maxF = deterministic_parallel_unitary_inversion(d,n,isComplex);
    case 2
        maxF = deterministic_sequential_unitary_inversion(d,n,isComplex);
end

maximal_fidelity = maxF

total_time_in_seconds = toc
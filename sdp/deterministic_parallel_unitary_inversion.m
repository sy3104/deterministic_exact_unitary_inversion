% -------------------------------------------------------------------------
%
% File: deterministic_parallel_unitary_inversion.m
%
% Description :
% Code to obtain the maximal fidelity of transforming parallel 'n' calls of
% any 'd'-dimensional unitary operations into its inverse map
%
% -------------------------------------------------------------------------

function maxF = deterministic_parallel_unitary_inversion(d,n,isComplex)

N = n+1;

run('young_diagrams.m')

cvx_solver sdpt3

cvx_begin SDP quiet

% ------------------------------------------------------------------
%          List of young diagrams with depth < d
% ------------------------------------------------------------------
for i = 1:N
    list_young_diagrams{i}={};
    for alpha = 1:num_young_diagrams{i}
        if depth{i}{alpha}<=d
            list_young_diagrams{i}{end+1} = alpha;
        end
    end
end

% ------------------------------------------------------------------
%          Calculate young diagrams by adding or removing a box
% ------------------------------------------------------------------
for i = 1:N
    list_child{i}={};
    list_parent{i}={};
    for index_alpha = 1:length(list_young_diagrams{i})
        alpha = list_young_diagrams{i}{index_alpha};
        if i<N
            list_child{i}{alpha}={};
            for index_mu=1:length(list_young_diagrams{i+1})
                mu = list_young_diagrams{i+1}{index_mu};
                if sum(X{i+1}{alpha}{mu},'all') ~= 0
                    list_child{i}{alpha}{end+1} = mu;
                end
            end
        end
        if i>1
            list_parent{i}{alpha}={};
            for index_gamma=1:length(list_young_diagrams{i-1})
                gamma = list_young_diagrams{i-1}{index_gamma};
                if sum(X{i}{gamma}{alpha},'all') ~= 0
                    list_parent{i}{alpha}{end+1} = gamma;
                end
            end
        end
    end
end

% ------------------------------------------------------------------
%          Declare SDP variables 
% ------------------------------------------------------------------
for index_mu = 1:length(list_young_diagrams{N})
    mu = list_young_diagrams{N}{index_mu};
    for index_nu = 1:length(list_young_diagrams{N})
        nu = list_young_diagrams{N}{index_nu};
        C{mu}{nu} = semidefinite(dim{N}{mu}*dim{N}{nu},isComplex);
    end
end

% ------------------------------------------------------------------
%          Conditions for the fidelity
% ------------------------------------------------------------------
for index_mu = 1:length(list_young_diagrams{N})
    mu = list_young_diagrams{N}{index_mu};
    PI = Tensor(perm{N}{mu},eye(dim{N}{mu}));
    Omega{mu} = PI*pure_to_mixed(MaxEntangled(dim{N}{mu}))*PI'*dim{N}{mu}/(d^2*mult{N}{mu});
end
F=0;
for index_mu = 1:length(list_young_diagrams{N})
    mu = list_young_diagrams{N}{index_mu};
    F = F + real(trace(C{mu}{mu}*Omega{mu}));
end

% ------------------------------------------------------------------
%          Parallel comb conditions
% ------------------------------------------------------------------
for index_alpha = 1:length(list_young_diagrams{N-1})
    alpha = list_young_diagrams{N-1}{index_alpha};
    D=0;
    for index_mu = 1:length(list_child{N-1}{alpha})
        mu = list_child{N-1}{alpha}{index_mu};
        for index_nu = 1:length(list_young_diagrams{N})
            nu = list_young_diagrams{N}{index_nu};
            XI = Tensor(X{N}{alpha}{mu},eye(dim{N}{nu}));
            D = D + PartialTrace(XI*C{mu}{nu}*XI', [2], [dim{N-1}{alpha} dim{N}{nu}]);
        end
    end
    for index_nu = 1:length(list_young_diagrams{N})
        LHS = 0;
        nu = list_young_diagrams{N}{index_nu};
        for index_mu = 1:length(list_child{N-1}{alpha})
            mu = list_child{N-1}{alpha}{index_mu};
            XI = Tensor(X{N}{alpha}{mu},eye(dim{N}{nu}));
            LHS = LHS + XI*C{mu}{nu}*XI'/mult{N}{nu};
        end
        LHS == Tensor(D,eye(dim{N}{nu}))/(d^(n+1));
    end
end
norm = 0;
for index_mu = 1:length(list_young_diagrams{N})
    mu = list_young_diagrams{N}{index_mu};
    for index_nu = 1:length(list_young_diagrams{N})
        nu = list_young_diagrams{N}{index_nu};
        norm = norm +trace(C{mu}{nu});
    end
end
norm == d^(n+1);

maximise F

cvx_end
cvx_status

maxF=F;

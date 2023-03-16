% -------------------------------------------------------------------------
%
% File: deterministic_sequential_unitary_inversion.m
%
% Description :
% Code to obtain the maximal fidelity of transforming sequential 'n' calls
% of any 'd'-dimensional unitary operations into its inverse map
%
% -------------------------------------------------------------------------

function maxF = deterministic_sequential_unitary_inversion(d,n,isComplex)

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
    permutation = Tensor(perm{N}{mu},eye(dim{N}{mu}));
    Omega{mu} = permutation*pure_to_mixed(MaxEntangled(dim{N}{mu}))*permutation'*dim{N}{mu}/(d^2*mult{N}{mu});
end
F=0;
for index_mu = 1:length(list_young_diagrams{N})
    mu = list_young_diagrams{N}{index_mu};
    F = F + real(trace(C{mu}{mu}*Omega{mu}));
end

% ------------------------------------------------------------------
%          Quantum comb conditions
% ------------------------------------------------------------------
for index_mu = 1:length(list_young_diagrams{N})
    mu = list_young_diagrams{N}{index_mu};
    for index_nu = 1:length(list_young_diagrams{N})
        nu = list_young_diagrams{N}{index_nu};
        D{N}{mu}{nu} = C{mu}{nu};
    end
end

for i = N:-1:2
    for index_gamma = 1:length(list_young_diagrams{i-1})
        gamma = list_young_diagrams{i-1}{index_gamma};
        for index_delta = 1:length(list_young_diagrams{i-1})
            delta = list_young_diagrams{i-1}{index_delta};
            D{i-1}{gamma}{delta} = 0;
            for index_alpha = 1:length(list_child{i-1}{gamma})
                alpha = list_child{i-1}{gamma}{index_alpha};
                for index_beta = 1:length(list_child{i-1}{delta})
                    beta = list_child{i-1}{delta}{index_beta};
                    XX = Tensor(X{i}{gamma}{alpha},X{i}{delta}{beta});
                    D{i-1}{gamma}{delta} = D{i-1}{gamma}{delta} + XX*D{i}{alpha}{beta}*XX'/d;
                end
            end
        end
    end
    for index_gamma = 1:length(list_young_diagrams{i-1})
        gamma = list_young_diagrams{i-1}{index_gamma};
        for index_beta = 1:length(list_young_diagrams{i})
            beta = list_young_diagrams{i}{index_beta};
            LHS = 0;
            RHS = 0;
            for index_alpha = 1:length(list_child{i-1}{gamma})
                alpha = list_child{i-1}{gamma}{index_alpha};
                XI = Tensor(X{i}{gamma}{alpha},eye(dim{i}{beta}));
                LHS = LHS + XI*D{i}{alpha}{beta}*XI'/mult{i}{beta};
            end
            for index_delta = 1:length(list_parent{i}{beta})
                delta = list_parent{i}{beta}{index_delta};
                IX = Tensor(eye(dim{i-1}{gamma}),X{i}{delta}{beta});
                RHS = RHS + IX'*D{i-1}{gamma}{delta}*IX/mult{i-1}{delta};
            end
            LHS == RHS;
        end
    end
end

trace(D{1}{1}{1}) == d;

maximise F

cvx_end
cvx_status

maxF=F;
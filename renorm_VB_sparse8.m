function [J, A, beta] = renorm_VB_sparse8(B, G, A, A0, g0, SG)

N = size(G, 2);
M = size(G, 1);
T = size(B, 2);

%% J Step

trtime(N, 0);   A1 = inv(A); trtime(N);
trtime(N, 1);   SG1 = inv(SG); trtime(N);
trtime(N, 2);   A1GT = A1 * G'; trtime(N); %#ok<MINV>
trtime(N, 3);   A1GTSB = A1GT / (G * A1GT + SG1); trtime(N);

trtime(N, 4);   J = A1GTSB * B; trtime(N);

% mean(abs(J))

% trtime(N, 2);   AGT = A * G'; trtime(N); 
% trtime(N, 3);   SB = sparse(G * AGT + SG); trtime(N);  % 0.12
% trtime(N, 4);   L = AGT / SB; trtime(N); 
% trtime(5);      J1 = L * B; trtime(N);

%fprintf(' THE DIFF : %e %e %e\n', min(abs(J)), min(abs(J1)), max(abs(J - J1)));

% beta
%mean(abs(B))
%mean(abs(G*J))
trtime(N, 51);  ERR = B - G*J; trtime(N);
%trtime(N, 52);  B1 = ERR' * SG * ERR; trtime(N);
trtime(N, 52);  B1 = sum(ERR'.*(SG * ERR)',2); trtime(N);
%trtime(N, 53);  B2 = J' * A * J; trtime(N);
trtime(N, 52);  B2 = sum(J'.*(A * J)',2); trtime(N);

gb = (M * T)/2;
% size(B1)
% size(B2)
trtime(N, 54);  beta = gb / (0.5 * sum (B1 + B2,1)); trtime(N);

%% Alpha Step

ga = g0 + T/2;

% diagonais do alpha
% I = speye(size(G,2));
% trtime(N, 61);  IA1GTSBG = I - A1GTSB * G; trtime(N);
% trtime(N, 62);  A1IA1GTSBG1 = sum(A1.*IA1GTSBG',2); trtime(N); 

trtime(N, 61);  A1A1GTSB = A1 * A1GTSB; trtime(N);
trtime(N, 62);  A1A1GTSBG = sum(A1A1GTSB.*G',2); trtime(N);
trtime(N, 63);  A1IA1GTSBG = diag(A1) - A1A1GTSBG ; trtime(N); 

% mean(abs(A1IA1GTSBG1))
% mean(abs(A1IA1GTSBG))
% max(A1IA1GTSBG1 - A1IA1GTSBG)

trtime(N, 7); 
for p = 1:N

    % alpha
%     (g0/a0)
%     J(p,:)
    
    A(p,p) = ga/ ( (g0 / A0(p,p)) + (beta/2) * sum( J(p,:).^2 ,2) + (T/2) * A1IA1GTSBG(p) );

end
trtime(N); 

function trtime(N, varargin)
    if N > 100000
        if nargin > 1
            s = '%d ';
            if nargin > 2
                s = '!%d ';
            end
            fprintf (s, varargin{1}); drawnow ; tic;
        else
            toc;
        end
    end
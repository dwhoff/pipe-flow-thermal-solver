function [G,lambda_sq] = get_eigenvalues(n,shape,Re)
% retrieves the eigenvalues and coefficients for the unheated starting
% length solution in a pipe.
% Currently only laminar round/square pipes are supported.

switch shape
    case 'round'
        if Re<2300
            if n<6
                lambda_sq_list = [7.312 44.62 113.8 215.2 348.5];
                G_list = [0.749 0.544 0.463 0.414 0.382];
                lambda_sq = lambda_sq_list(n);
                G = G_list(n);
            else
                lambda_sq = (4*(n-1)+8/3)^2;
                G = 1.01276*lambda_sq^(-1/6);
            end
        else
            error('Flow is turbulent: Eigenvalues not available')
        end
    case 'rectangle'
        if Re<2300
            load('eigenvalues_square.mat')
            lambda_sq_list = lambda_sq;
            G_list = G;
            if n<=numel(G_list)
                lambda_sq = lambda_sq_list(n);
                G = G_list(n);
            else
                error([sprintf('Warning: Not enough eigenvalues available.\n') ...
                    sprintf('Please decrease the number of elements to improve accuracy.')])
            end
        else
            error('Flow is turbulent: Eigenvalues not available')
        end
end

end
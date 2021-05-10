function [stepSize,breakCond] = backTrackingLineSearch(fun,startPos,gradient,initStepSize,ctrlParam)
%backTrackingLineSearch Backtracking line search
%   Method based on the Armijo-Goldstein condition to find an adequate step
%   size for gradient descent minimization
% fun:          Objective function to minimize
% startPos:     Initial guess 
% gradient:     Gradient of objective function evaluated at startPos
% initStepSize: Initial step size
% ctrlParam:    Control parameter (backtracking parameter [0.1,0.8])
    stepSize = initStepSize; % initial stepsize
    f1 = fun(startPos);
    f2 = fun(startPos-stepSize*gradient);
    i = 1;
    breakCond = 0;
    while f1-f2 < ctrlParam*stepSize*norm(gradient(:))^2
        stepSize = 0.5*stepSize;
        f2 = fun(startPos-stepSize*gradient);
        i = i + 1;
        if i == 50
            breakCond = 1;
            break;
%             error('Error: backtracking doesn''t seem to find a good step size. Maximum steps reached.')
        elseif stepSize < 1e-12
            breakCond = 1;
            break;
%             error('Error: backtracking doesn''t seem to find a good step size. Estimation may be already good enough.')
        end
    end
end


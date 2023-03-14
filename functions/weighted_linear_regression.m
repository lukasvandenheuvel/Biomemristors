function [b,S] = weighted_linear_regression(x,y,w)
    % S is weighted objective function
    b = ((x.*w)'*x)\((x.*w)'*y);
    S = norm(sqrt(w) .* (y - b*x))^2;
end

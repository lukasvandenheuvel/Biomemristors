function [b,R2] = linear_regression(x,y,add_constant)
    arguments
        x double;
        y double;
        add_constant logical = false;
    end
    if add_constant
        x = [ones(size(x,1),1),x];
    end

    b  = x\y;
    yPred = x*b;
    R2 = 1 - sum((y - yPred).^2)/sum((y - mean(y)).^2);
end
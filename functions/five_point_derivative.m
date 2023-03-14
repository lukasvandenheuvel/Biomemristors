function dy = five_point_derivative(x)
    dy = (-x(5:end) + 8*x(4:end-1) - 8*x(2:end-3) + x(1:end-4)) / 12;
end
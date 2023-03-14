function plot_fill(x,y,error,fill_color)
    % Plot a line plot with errorbars
    if size(x,1)==1
        x = x';
    end
    if size(y,1)==1
        y = y';
    end
    if size(error,1)==1
        error = error';
    end
    % Check if the arrays are now the right size
    assert(size(x,1)~=1 & size(x,2)==1, 'Input array x must be of size (N,1)')
    assert(size(y,1)~=1 & size(y,2)==1, 'Input array y must be of size (N,1)')
    assert(size(error,1)~=1 & size(error,2)==1, 'Input array error must be of size (N,1)')
    
    xaxis = [x; flipud(x)];
    errobars = [y-error;... 
                flipud(y+error)];
    fill(xaxis,errobars,fill_color, 'linestyle', 'none')
end
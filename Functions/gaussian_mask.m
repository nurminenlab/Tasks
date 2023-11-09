function Gmask = gaussian_mask(xscreen,yscreen,xoff,yoff,xSigma,ySigma,theta)
  
% make the coordinates
[X,Y] = meshgrid(linspace(-xscreen/2,xscreen/2,xscreen),linspace(-yscreen/2,yscreen/2,yscreen));
% translate the coordinate

a = cos(theta)^2 / (2 * xSigma^2) + sin(theta)^2 / (2 * ySigma^2);
b = sin(2 * theta) / (4 * xSigma^2) - sin(2 * theta) / (4 * ySigma^2);
c = sin(theta)^2 / (2 * xSigma^2) + cos(theta)^2 / (2 * ySigma^2);

Gmask = exp(-(a * (X - xoff).^2 + 2 * b * (X - xoff) .* (Y - yoff) + c * (Y - yoff).^2));

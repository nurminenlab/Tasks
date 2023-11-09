function [cos_window] = raised_cosine(plat,edge,xscreen,yscreen,xoff,yoff,sw)
% [ cos_window ] = raised_cosine(plat,edge,xscreen,yscreen,xoff,yoff,sw)
% plat, plateau in pixels
% edge, full width of a single edge in pixels
% xscreen yscreen, display size in pixels
% xoff yoff, mask center offset from display center in pixels
% sw, character switch that makes pathces or donuts, 'R' = rod, 'H' =
% hole

% make the coordinats
[X,Y] = meshgrid(linspace(-xscreen/2,xscreen/2,xscreen), ...
                            linspace(-yscreen/2,yscreen/2,yscreen));
% translate the coordinate
X = X - xoff;
Y = Y + yoff;

[THETA,RHO] = cart2pol(X,Y);

cosini=(cos((RHO-plat/2).*(pi/edge))+1)/2;
colorado_plateau = double(RHO < plat/2);
salt_lake_desert = double(RHO > plat/2 + edge);
cosini(colorado_plateau==1) = colorado_plateau(colorado_plateau==1);
cosini(salt_lake_desert==1) = 0;

% case either R(od) or H(ole)
if sw == 'R'
    cos_window=cosini;
elseif sw == 'H'
    cos_window=cosini;        
    cos_window = 1-cos_window;
end

end
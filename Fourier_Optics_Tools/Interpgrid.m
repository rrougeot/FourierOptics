function [vi]=Interpgrid(vx,vy,I,xi,yi)

vx=vx(:);
vy=vy(:);


dx=vx(2)-vx(1);
dy=vy(2)-vy(1);

ix=uint64(floor((xi-vx(1))/dx)+1);
iy=uint64(floor((yi-vy(1))/dy)+1);

deltax=(xi-vx(ix))/dx;
deltay=(yi-vy(iy))/dy;

vi= I(sub2ind(size(I),ix,iy)).*(1-deltax).*(1-deltay)+...
    I(sub2ind(size(I),ix+1,iy)).*(deltax).*(1-deltay)+...
    I(sub2ind(size(I),ix,iy+1)).*(1-deltax).*deltay+...
    I(sub2ind(size(I),ix+1,iy+1)).*deltax.*deltay;
end
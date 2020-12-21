function Ry = ry(theta)
c = cos(theta);
s = sin(theta);

Ry = [c 0 s;
      0 1 0;
      -s 0 c];
end
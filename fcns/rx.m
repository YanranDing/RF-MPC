function Rx = rx(phi)
c = cos(phi);
s = sin(phi);

Rx = [1 0 0;
      0 c -s;
      0 s  c];
end
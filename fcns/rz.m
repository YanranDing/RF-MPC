function Rz = rz(psi)
c = cos(psi);
s = sin(psi);

Rz = [c -s 0;
      s c  0;
      0 0  1];
end
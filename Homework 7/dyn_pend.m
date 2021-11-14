function dz = dyn_pend(z, u, param)
  
  A = A_pend(z, param);
  b = b_pend(z, u, param);

  
  ddq = inv(A)*b;
  dim = length(z);
  dz = zeros(dim,1);
  dz(1:dim/2) = z(dim/2+1:end);
  dz(dim/2+1:end) = ddq;
  
end
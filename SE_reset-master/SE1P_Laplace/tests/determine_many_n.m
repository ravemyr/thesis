% Determine n for multiple values of M

clear

Mv = [12:4:40];

for M=Mv
  [nG, nK] = determine_single_n(M, false);
  fprintf('  %d    %d           %d\n', M, nG, nK);
end

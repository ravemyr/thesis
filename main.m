N = 1000;
avg = 0;
for i = 1:N
    run('first_ewald.m');
    avg = avg + ans;
end
res = avg/N
td = load('testdata5.txt');
N = 1000;
xi = 20;
L = 2;
tm = [];
te = [];
for i = 1:length(td(:,1))
    if(td(i,[1,2,3]) == [N,L,xi])
        tm = [tm,td(i,6)];
        te = [te,td(i,5)];
    end
end
semilogy(tm,te,'*-')
hold on
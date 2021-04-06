close all;
td = load('testdata7.txt');
[td_tols, idx] = sort(td(:,4));
td_errs = td(idx,5);
semilogx(flip(td_tols),1:length(td(:,1)),'ro')
hold on 
semilogx(flip(td_errs),1:length(td(:,1)),'b*')
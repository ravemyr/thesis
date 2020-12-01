function Q = rmse(x)
    N = max(size(x));
    QQ = 0;
    for i=1:N
       QQ = QQ + norm(x(i,:))^2; 
    end
    QQ = QQ/N;
    Q = sqrt(QQ);
end
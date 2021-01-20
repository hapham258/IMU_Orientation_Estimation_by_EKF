function zk = measurement_model(xk)
    g0 = 9.7867;
    m0 = [0.9876 0.0121 0.1567]';
    q1k = xk(1);
    q2k = xk(2);
    q3k = xk(3);
    q4k = xk(4);
    omegak = xk(5:7);
    bfk = xk(8:10);
    bmk = xk(11:13);
    bgk = xk(14:16);
    
    Cq = [ q1k^2-q2k^2-q3k^2+q4k^2  2*(q1k*q2k+q3k*q4k)      2*(q1k*q3k-q2k*q4k)
           2*(q1k*q2k-q3k*q4k)     -q1k^2+q2k^2-q3k^2+q4k^2  2*(q2k*q3k+q4k*q1k);
           2*(q1k*q3k+q2k*q4k)      2*(q2k*q3k-q4k*q1k)     -q1k^1-q2k^2+q3k^2+q4k^4];
    zk = zeros(9, 1);
    zk(1:3) = g0*Cq(:, 3) + bfk;
    zk(4:6) = Cq*m0 + bmk;
    zk(7:9) = omegak + bgk;
end
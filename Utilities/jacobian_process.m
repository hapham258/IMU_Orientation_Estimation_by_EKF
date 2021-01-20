function F = jacobian_process(xk, Ts)
    tauf = 200;
    taug = 300;
    ek = xk(1:3);
    q4k = xk(4);
    omegak = xk(5:7);
    
    S_omegak = skew_symmetric_matrix(omegak);
    S_ek = skew_symmetric_matrix(ek);
    tmp = [-S_omegak    omegak       S_ek+q4k*eye(3)  zeros(3, 3)    zeros(3, 3)   zeros(3, 3);
           -omegak'     0           -ek'              zeros(1, 3)    zeros(1, 3)   zeros(1, 3);
            zeros(3, 3) zeros(3, 1)  zeros(3, 3)      zeros(3, 3)    zeros(3, 3)   zeros(3, 3);
            zeros(3, 3) zeros(3, 1)  zeros(3, 3)     -2/tauf*eye(3)  zeros(3, 3)   zeros(3, 3);
            zeros(3, 3) zeros(3, 1)  zeros(3, 3)      zeros(3, 3)    zeros(3, 3)   zeros(3, 3);
            zeros(3, 3) zeros(3, 1)  zeros(3, 3)      zeros(3, 3)    zeros(3, 3)  -2/taug*eye(3)]; 
    F = eye(16) + Ts/2*tmp;
end
function H = jacobian_measurement(xk)
    g0 = 9.7867;
    m0 = [0.9876 0.0121 0.1567]';
    q1k = xk(1);
    q2k = xk(2);
    q3k = xk(3);
    q4k = xk(4);
    
    dc_3_dqk = 2*[ q3k -q4k  q1k -q2k;
                   q4k  q3k  q2k  q1k;
                  -q1k -q2k  q3k  q4k];
    dc1__dqk = 2*[ q1k -q2k -q3k  q4k;
                   q2k  q1k  q4k  q3k;
                   q3k -q4k  q1k -q2k];
    dc2__dqk = 2*[ q2k  q1k -q4k -q3k;
                  -q1k  q2k -q3k  q4k;
                   q4k  q3k  q2k  q1k];
    dc3__dqk = 2*[ q3k  q4k  q1k  q2k;
                  -q4k  q3k  q2k -q1k;
                  -q1k -q2k  q3k  q4k];
    dh1_dqk = g0*dc_3_dqk;
    dh2_dqk = [dc1__dqk'*m0 dc2__dqk'*m0 dc3__dqk'*m0]';
    H = [dh1_dqk     zeros(3, 3) eye(3)      zeros(3, 3) zeros(3, 3);
         dh2_dqk     zeros(3, 3) zeros(3, 3) eye(3)      zeros(3, 3);
         zeros(3, 4) eye(3)      zeros(3, 3) zeros(3, 3) eye(3)];
end
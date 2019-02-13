BEGIN{
    pi = atan2(1,1)*4.0
}{
    # Solve for x given a3*x^3 + a2*x^2 + a1*x +a0 = 0
    # Input is: a3 a2 a1 a0
    a3 = $1
    a2 = $2
    a1 = $3
    a0 = $4

    # Normalize equation so that a3 = 1
    a0 = a0/a3
    a1 = a1/a3
    a2 = a2/a3
    a3 = a3/a3

    Q = (3*a1-a2*a2)/9
    R = (9*a2*a1-27*a0-2*a2*a2*a2)/54

    D = Q*Q*Q + R*R

    # Depending on number of roots, calculate differently
    if (D<=0) {
        print "D<=0: ",D
        x = R/sqrt(-Q*Q*Q)
        theta = atan2(sqrt(1-x*x), x) 
        x1 = 2*sqrt(-Q)*cos(theta/3) - a2/3
        x2 = 2*sqrt(-Q)*cos((theta+2*pi)/3) - a2/3
        x3 = 2*sqrt(-Q)*cos((theta+4*pi)/3) - a2/3
        print x1,x2,x3
    } else {
        print "D POSITIVE: ",D
        S = (R+sqrt(D))^(1/3)
        T = (R-sqrt(D))^(1/3)
        x1 = S + T - a2/3
        print x1
    }
}END{

}

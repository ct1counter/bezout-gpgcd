# gpgcd-bezout-test 2nd variation
# GPGCD test routines for Bezout resultant
# Akira Terui, Boming Chi, October 2018

with (LinearAlgebra);
with (PolynomialTools);
with (ArrayTools);
with (VectorCalculus);

# The update function
bezout_update := proc (Btemplate, J, degH, F0, G0, V0, F1, G1, V1, alpha)
local i, j, dimF, dimG, dimV, dimVV1, F1poly, G1poly, B1, VV1, Q1, P1, R1, dim, BB, vlist, varlist, F0grad, F1grad, subslistB, subslistV, sublist, J1, Jrow, Jcol, A1, D1, DD1, DD2, VV2, F2, G2, V2, Normfg, st;
    # Jrow, Jcol, A1, R1, D1, f, k, B2, V2;

    # inputs:

    # B:            Bezout matrix (template)
    # J:             Jacobi matrix (template)
    # degH:     degree of the approximate GCD
    # F0, G0:    initial coefficient vector of F and G
    # F1, G1:    current coefficient vector of F and G
    # V0:          initial singular vector of the Bezout matrix
    # V1:          current singular vector of the Bezout matrix
    # alpha:      step width for an iteration

    # outputs:
    
    # F2, G2:    new (updated) coefficient vector of F and G
    # V2:          new (updated) singular vector of the Bezout matrix
    # Normfg:  the distance between F1, G1 and F2, G2

    # find dimension of initial coefficient vector( F0, G0 and V0)

    dimF := Dimension(F0);
    dimG := Dimension(G0);
    dimV := Dimension(V0);
    userinfo(2, bezout_update, `[dimF, dimG, dimV] =`,
             print([dimF, dimG, dimV]));

    # construct polynomials from coefficient vector( F1 and G1)

    st := time();
    F1poly := FromCoefficientVector(F1, x);
    G1poly := FromCoefficientVector(G1, x);
    st := time() - st;
    userinfo(1, bezout_update, `stpoly=`, print(st));
    userinfo(2, bezout_update, `F1poly=`, print(F1poly));
    userinfo(2, bezout_update, `G1poly=`, print(G1poly));

    # construct substitution table for B1 & J1

    st := time();
    sublist := [seq(f[j] = F1[j+1], j=0..(dimF-1)), seq(g[j] = G1[j+1], j=0..(dimG-1)), seq(v[j] = V1[j], j=1..dimV)];
    st := time() - st;
    userinfo(1, bezout_update, `stsublist=`, print(st));
    userinfo(2, bezout_update, `sublist = `, print(sublist, whattype(sublist)));

    # construct Bezout Matrix of F1 and G1 -> B1

    st := time();
    B1 := BezMatrix(F1poly, G1poly, x);
    userinfo(2, bezout_update, `B1 =`, print(B1));
    st := time() - st;
    userinfo(1, bezout_update, `stB1=`, print(st));

    # construct vector of objective function

    st := time();
    VV1 := <F1, G1, V1>;
    userinfo(2, bezout_update, `VV1 =`, print(VV1, whattype(VV1)));
    dimVV1 := dimF + dimG + dimV;
    st := time() - st;
    userinfo(1, bezout_update, `stVV1=`, print(st));

    # evaluate the constraint value -> Q1
    
    st := time();
    Q1 := B1 . V1[ 1 .. dimF-1 ];
    st := time() - st;
    userinfo(1, bezout_update, `stQ1=`, print(st));
    userinfo(2, bezout_update, `Q1 =`, print(Q1, whattype(Q1)));

    st := time();
    P1 := < (F1 - F0), (G1 - G0), Vector(dimV, datatype=double) >;
    st := time() - st;
    userinfo(1, bezout_update, `stP1=`, print(st));
    userinfo(2, bezout_update, `P1 =`, print(P1, whattype(P1)));

    # construct right-hand-side vector -> R1

    st := time();
    R1 := - <P1, Q1>;
    st := time() - st;
    userinfo(1, bezout_update, `stR1=`, print(st));
    userinfo(2, bezout_update, `R1 =`, print(R1, whattype(R1)));



    # construct Jacobian matrix -> J1

    st := time();
    J1 := subs(sublist, J);
    st := time() - st;
    userinfo(1, bezout_update, `stJ1=`, print(st));
    userinfo(2, bezout_update, `J1 = `, print(J1, whattype(J1)));

    # construct coefficient matrix for the linear system -> A1

    st := time();
    Jrow, Jcol := Dimension(J);
    A1 := Matrix(Jrow+Jcol);
    for i from 1 to Jcol do A1(i,i):=1 end do:
    for i from 1 to Jrow do for j from 1 to Jcol do A1(i+Jcol,j):=J1(i,j) end do end do:
    for i from 1 to Jcol do for j from 1 to Jrow do A1(i,j+Jcol):=-J1(j,i) end do end do:
    st := time() - st;
    userinfo(1, bezout_update, `stA1=`, print(st));
    userinfo(2, bezout_update, `A1 =`, print(A1));

    # solve the linear system: A1 . D1 == R1

    st := time();
    D1 := LinearSolve(A1, R1);
    st := time() - st;
    userinfo(1, bezout_update, `stD1=`, print(st));
    userinfo(2, bezout_update, `D1 =`, print(D1));

    DD1 := D1[1..dimVV1];
    userinfo(2, bezout_update, `DD1 =`, print(DD1));

    DD2 := D1[1..(dimF+dimG)];
    userinfo(2, bezout_update, `DD2 =`, print(DD2));
    
    st := time();
    VV2 := (VV1 + alpha * DD1)[..,1];
    st := time() - st;
    userinfo(1, bezout_update, `stVV2=`, print(st));
    userinfo(2, bezout_update, `VV2 =`, print(VV2, whattype(VV2)));

    F2 := VV2[1..dimF];
    G2 := VV2[(dimF + 1)..(dimF + dimG)];
    V2 := VV2[(dimF + dimG + 1)..(dimVV1)]; 
    Normfg := LinearAlgebra:-Norm(DD2,2);
    userinfo(2, bezout_update, `Normfg =`, print(Normfg));
    
    return F2, G2, V2, Normfg;
    
end proc:

# The initial function
bezout_gpgcd := proc(F, G, x, degH, stopcriterion, alpha)
local degF, degG, hh, Ftemplate, Gtemplate, Btemplate, V, BV, i, j, m, n, J, Fcoef, Gcoef, B, V0, VV0, v0, st, st0, st1, st2, ff, gg, DF, DG, DD, numofiteration, sublist, B0;

    # inputs: 

    # F: the initial polynomial F
    # G: the initial polynomial G
    # degH: the degree of the approximate GCD
    # alpha: step width for an iteration

    # outputs:

    # ff, gg: the polynomials finally we get
    # hh: the approximate GCD finally we get
    # DF, DG: the perturbation of ff and gg
    # DD: the whole perturbation
    # numofiteration: number of iterations
    # B0: the Bezout Matrix of ff and gg

    # Construct the Bezout matrix of the initial polynomials

    st0 := time();
    degF := degree(F,x);
    degG := degree(G,x);
    # if degF < degG then degF, degG := degG, degF; F, G := G, F end if;
    userinfo(2, bezout_update, `degG =`, print(degG));
    Ftemplate := add(f[i]*x^i, i=0..degF);
    Gtemplate := add(g[i]*x^i, i=0..degG);
    userinfo(2, bezout_update, `Gtemplate =`, print(Gtemplate));    st := time();
    Btemplate := BezMatrix(Ftemplate, Gtemplate, x);

    # Construct the Jacobian matrix

    st := time() - st;
    userinfo(1, bezout_update, `stBtemplate=`, print(st));
    V := Vector(degF,symbol=v);
    BV := Vector(degF);
    for i from 1 to degF do BV[i] := add(Btemplate[i,j]*V[j], j=1..degF) end do;
    userinfo(2, bezout_update, `BV =`, print(BV));
    st := time();
    J := bezout_jacobian(degF, degG, Btemplate);
    st := time() - st;
    userinfo(1, bezout_update, `stJtemplate=`, print(st));
    userinfo(1, bezout_update, `J =`, print(J));
    
    # Construct coefficient vectors

    Fcoef := PolynomialTools:-CoefficientVector(F, x);
    Gcoef := PolynomialTools:-CoefficientVector(G, x);
    
    # Construct singular vector of the d-th from the right side

    B := BezMatrix(F, G, x);
    V0 := LinearAlgebra:-Transpose(LinearAlgebra:-SingularValues(B, output='Vt'));
    v0 := V0[..,-degH];
    userinfo(2, bezout_update, `v0 =`, print(v0));
    st1 := time() - st0;
    userinfo(1, bezout_update, `stbeforeiteration=`, print(st1));
    ff, gg, hh, DF, DG, DD, numofiteration, B0 := bezout_iteration(Btemplate, J, degF, degG, degH, Fcoef, Gcoef, v0, Fcoef, Gcoef, v0, stopcriterion, 1.0, 0);
    st2 := time() - st0;
    userinfo(1, bezout_update, `stwholeiteration=`, print(st2));

    return ff, gg, hh, DF, DG, DD, numofiteration, B0;
    # return B, J, V0, v0
end proc:

# The Jacobian Matrix
bezout_jacobian := proc (m, n, Btemplate)
local J, i, j;

    # inputs: 
    # degF, degG： the degree of polynomials F and G

    J := Matrix(m, 2*m + n + 2);

    for i from 1 to m do for j from i+1 to m+1 do J(i,j):= -g[i-1]*v[j-1] end do end do;
    for i from 2 to m do for j from i+1 to m+1 do J(i,j):= J(i,j)+J(i-1,j-1) end do end do;

    for j from 1 to n do for i from j to n do J(i,j):= g[i]*v[j] end do end do;
    for j from n-1 by -1 to 1 do for i from n-1 by -1 to j do J(i,j):= J(i,j)+J(i+1,j+1) end do end do;
    
    for i from 1 to n+1 do for j from i+m+2 to m+n+2 do J(i,j):= f[i-1]*v[j-m-2] end do end do;
    for i from 2 to n+1 do for j from i+m+2 to m+n+2 do J(i,j):= J(i,j)+J(i-1,j-1) end do end do;

    for j from m+2 to 2*m+1 do for i from j-m-1 to m do J(i,j):= -f[i]*v[j-m-1] end do end do;
    for j from 2*m by -1 to m+2 do for i from m-1 by -1 to j-m-1 do J(i,j):= J(i,j)+J(i+1,j+1) end do end do;

    for j from m+n+3 to 2*m+n+2 do for i from 1 to m do J(i,j):= Btemplate(i,j-m-n-2) end do end do;
    return J
end proc:

# The iteration function
bezout_iteration := proc (Btemplate, J, degF, degG, degH, F0, G0, V0, F1, G1, V1, stopcriterion, alpha, numofiteration)
local FF, GG, VV, Normfg, B, ff, gg, f0, g0, DF, DG, DD, hh, st, numofiteration0, sublist;

    # inputs: 

    # Btemplate:       the template of the Bezout Matrix
    # J:                      the template of the Jacobian Matrix
    # degF, degG:     the degree of the initial polynomials
    # degH:               the degree of the approximate GCD
    # F0, G0:              initial coefficient vector of F and G
    # F1, G1:              current coefficient vector of F and G
    # V0:                    initial singular vector of the Bezout matrix
    # V1:                    current singular vector of the Bezout matrix
    # alpha:                step width for an iteration
    # stopcriterion:     the stop criterion
    # numofiteration: number of current iteration( initial number of iteration: 0)

    # outputs:

    # ff, gg: the polynomials finally we get
    # hh: the approximate GCD finally we get
    # DF, DG: the perturbation of ff and gg
    # DD: the whole perturbation
    # numofiteration: number of current iteration
    # B0: the Bezout Matrix of ff and gg

    FF, GG, VV, Normfg := bezout_update(Btemplate, J, degH, F0, G0, V0, F1, G1, V1, alpha);
    if Normfg < stopcriterion then
    userinfo(2, bezout_update, `FF=`, print(FF));  
    st := time();
    f0 := FromCoefficientVector(F0,x);
    g0 := FromCoefficientVector(G0,x);
    ff := FromCoefficientVector(FF,x);
    gg := FromCoefficientVector(GG,x);
    DF := Norm(PolynomialTools:-CoefficientVector(ff-f0,x),2);
    DG := Norm(PolynomialTools:-CoefficientVector(gg-g0,x),2);
    DD := sqrt(DF^2+DG^2);
    sublist := [seq(f[j] = FF[j+1], j=0..degF), seq(g[j] = GG[j+1], j=0..degG)];
    B := BezMatrix(ff, gg, x);
    userinfo(2, bezout_update, `B=`, print(B));  
    hh := bezout_gcd(B, max(degF, degG), degH);
    st := time() - st;
    userinfo(1, bezout_update, `stafteriteration=`, print(st));  
        return ff, gg, hh, DF, DG, DD, numofiteration, B;  
    else
    numofiteration0 := numofiteration + 1;
    userinfo(2, bezout_update, `numofiteration0 =`, print(numofiteration0));
    return bezout_iteration(Btemplate, J, degF, degG, degH, F0, G0, V0, FF, GG, VV, stopcriterion, alpha, numofiteration0); end if
end proc:

    # Find Bezout Matrix from coefficient vector

BezMatrix := proc(F, G, x)
local degF, degG, Bez, Bez1, i, j, fv, gv, gv0;

    # inputs: 
    # F:    the polynomial F
    # G:    the polynomial G

    # outputs:

    # Bez:    the Bezout Matrix of F and G

    degF := degree(F, x);
    degG := degree(G, x);  
    fv := CoefficientVector(F, x);
    gv0 := CoefficientVector(G, x);
    gv := Vector(degF+1);
    for i from 1 to degG+1 do gv(i):=gv0(i) end do;

    # let the degree of F be the larger one

    #if degF < degG then F, G := G, F ; degF, degG := degG, degF end if;

    # construct the Bezout Matrix

    Bez := Matrix(degF);

    for j from 1 to degF do 
       for i from 1 to j do
          Bez(i,j) := fv[i]*gv[j+1]-gv[i]*fv[j+1]
       end do
    end do;

    for i from 2 to degF-1 do
       for j from i to degF-1 do
          Bez(i,j) := Bez(i,j)+Bez(i-1,j+1)
        end do
    end do;

    for j from 1 to degF-1 do 
       for i from j+1 to degF do
          Bez(i,j) := Bez(j,i)
       end do
    end do;

    return Bez
end proc:

# Find norm between two polynomials

normpoly := proc(f, g, x)
local n;

    n := Norm(CoefficientVector(f-g, x),2);
    return n
end proc:

# Find GCD of polynomials from their Bezout Matrix

bezout_gcd := proc (B, n, h)
local p, l, u, b, i, j, X, k, y, CoeffVector, gcdB;

    # inputs: 
    # B: the Bezout Matrix
    # n: the degree of the Bezout Matrix
    # k: the degree of GCD

    # outputs:

    # gcdB: the GCD calculated by the Bezout Matrix

    p,l,u:=LUDecomposition(B[..,h+1..n]);
    CoeffVector := Vector(h+1);
    CoeffVector(h+1) := 1;
    for k from 1 to h do
        b:=Transpose(p).B[..,k];
        y:=Vector(1..n);
        for i from 1 to n do y(i):=b(i)-add(l[i,j]*y[j],j=1..i-1) end do:
        X:=Vector(1..n-h);
        for i from 1 to n-h do X(n-h+1-i):=(y(n-h+1-i)-add(u[n-h+1-i,j]*X[j],j=n-h+2-i..n-h))/u(n-h+1-i,n-h+1-i) end do:
        CoeffVector(k) := X(1);
    end do:
    #p,l,u:=LUDecomposition(B[1..n-h,h+1..n]);
    #CoeffVector := Vector(h+1);
    #CoeffVector(h+1) := 1;
    #for k from 1 to h do
    #    LinearSolve([p,l,u],B[..,k],method=LU);
    #end do:
    gcdB := FromCoefficientVector(CoeffVector,x);
    return gcdB
end proc:


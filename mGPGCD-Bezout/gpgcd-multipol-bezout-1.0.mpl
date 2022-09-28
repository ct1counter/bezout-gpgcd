# gpgcd-bezout-test with more than 2 univariate polynomials
# GPGCD test routines for Bezout resultant
# Akira Terui, Boming Chi, May 2021

# Update with the initial values

with(LinearAlgebra):
with(PolynomialTools):
with(ArrayTools):
with(VectorCalculus):
with(StringTools):
with(ExcelTools):

# Input the initial values
bezout_gpgcd_mul := proc(F0, u, numpol0, var, d0, stopcriterion, alpha, num_bound, exname)
local 
    numpol,
    pol, 
    sigma, # weight of the singular values of R1 in bezout_update (the modified Newton method)
    d, 
    hh, 
    F, 
    G, 
    Btemplate, 
    Atemplate,
    btemplate,
    V, 
    BV, 
    i, 
    j, 
    m, 
    n, 
    J, 
    pol0,
    polvect0,
    B, 
    A,
    b,
    V0, 
    X0,
    VV0, 
    v0, 
    st, 
    st0, 
    st1, 
    st2, 
    f2, 
    data,
    ff, 
    gg, 
    DF, 
    DG, 
    DD, 
    numofiteration, 
    sublist, 
    polvect2, V2, Normpol,
    B0, f1, norm_perturbation, norm_remainder, h0;

global excelname;

    # inputs: 

    # pol: the initial polynomials
    # var: the variable
    # d0: the degree of the approximate GCD
    # stopcriterion: the stop criterion with the modified Newton method
    # alpha: step width for an iteration    
    # num_bound: the boundary of number of iterations

    # outputs:

    # ff, gg: the polynomials finally we get
    # hh: the approximate GCD finally we get
    # DF, DG: the perturbation of ff and gg
    # DD: the whole perturbation
    # numofiteration: number of iterations
    # B0: the Bezout Matrix of ff and gg

    excelname := exname;

    # Construct the Bezout matrix of the initial polynomials

    st0 := time();

    # Count the number of the polynomials

    pol := convert(F0[u, 1..numpol0], list);
    numpol := nops(pol);
    userinfo(2, bezout_gpgcd, `numpol =`, print(numpol));

    # Show the weight of the first polynomial

    sigma := norm(pol[1], 2);
    userinfo(2, bezout_gpgcd, `sigma =`, print(sigma));

    # Find the degrees of the polynomials

    d := Array(1 .. numpol);
    for i to numpol do
        d[i] := degree(pol[i], var);
    end do;
    userinfo(2, bezout_gpgcd, `d =`, print(d));

    m := d[1]; # The largest degree of the polynomials
    userinfo(2, bezout_gpgcd, `m =`, print(m));

    F:=Matrix(numpol,m+1,symbol=f);
    G:=Array(1..numpol);
    for i to numpol do G[i]:=Array(F[i,..]) end do;

    # Find the template of the Bezout matrix
    
    Btemplate := BezMatrix_2(G[1],G[2]);
    for i from 3 to numpol do Btemplate:=<Btemplate, BezMatrix_2(G[1],G[i])> end do;
    Atemplate := Btemplate[..,1..m-d0];
    btemplate := Btemplate[..,m-d0+1];

    # Construct the Jacobian matrix

    st := time();
    J := bezout_jacobian(d, d0, numpol, Atemplate);
    st := time() - st;
    userinfo(2, bezout_gpgcd, `stJtemplate=`, print(st));
    userinfo(2, bezout_gpgcd, `J =`, print(J));
    
    sublist := [seq(seq(f[i,j]=CoefficientVector(pol[i], x)[j], j=1..d[i]+1),i=1..numpol), seq(seq(f[i,j]=0, j=d[i]+2..d[1]+1),i=1..numpol)];
    pol0 := subs(sublist, G);
    
    # Construct the Bezout matrix

    B := BezMatrix_mul(pol, var);
    
    # Find the first m-d0 column vectors of B

    A := B[..,1..m-d0];

    # Find the m-d0+1-th column vector of B

    b := B[..,m-d0+1];

    # Find the initial X with AX=b with least squares

    X0 := LeastSquares(A, b, optimize);
    polvect0 := convert(convert(pol0, listlist), Vector);
    userinfo(2, bezout_gpgcd, `A =`, print(A));
    userinfo(2, bezout_gpgcd, `b =`, print(b));
    userinfo(2, bezout_gpgcd, `X0 =`, print(X0));
    st1 := time() - st0;
    userinfo(1, bezout_gpgcd, `stbeforeiteration=`, print(st1));
    polvect2, V2, Normpol, numofiteration := bezout_iteration(Atemplate, btemplate, J, d, d0, numpol, sigma, polvect0, X0, polvect0, X0, stopcriterion, alpha, 0, num_bound, u, var, pol);

    st := time();
    f1 := Array(1..numpol):
    for i to numpol do f1[i] := polvect2[(d[1]+1)*(i-1)+1..(d[1]+1)*i] end do:
    B0 := BezMatrix_2(f1[1],f1[2]):
    for i from 3 to numpol do B0 := <B0, BezMatrix_2(f1[1],f1[i])> end do:
    h0 := bezout_gcd(B0, d[1], d0);
    norm_perturbation := Vector(numpol);
    for i to numpol do norm_perturbation[i] := norm(FromCoefficientVector(polvect2[(d[1]+1)*(i-1)+1..(d[1]+1)*i], var)-pol[i],2) end do;
    
    norm_remainder := Vector(numpol);
    for i to numpol do norm_remainder[i] := norm(SNAP:-Remainder(FromCoefficientVector(f1[i], var), h0, var),2) end do;
    
    st := time() - st;
    userinfo(1, bezout_iteration, `stafteriteration=`, print(st));  
    st2 := time() - st0;
    userinfo(1, bezout_gpgcd, `stwholeiteration=`, print(st2));
    
    f1 := Array(1..numpol):
    for i to numpol do f1[i] := polvect2[(d[1]+1)*(i-1)+1..(d[1]+1)*i] end do:
    B0 := BezMatrix_2(f1[1],f1[2]):
    for i from 3 to numpol do B0 := <B0, BezMatrix_2(f1[1],f1[i])> end do:
    h0 := bezout_gcd(B0, d[1], d0);
    f2 := Array(1..numpol):
    for i to numpol do f2[i] := expand(h0*remainder_poly(pol[i], h0, var)[3]) end do:
    norm_perturbation := Vector(numpol);
    for i to numpol do norm_perturbation[i] := norm(f2[i]-pol[i],2) end do;
    userinfo(2, bezout_iteration, `f1 =`, print(f1));
    userinfo(2, bezout_iteration, `h0 =`, print(h0));
    norm_remainder := Vector(numpol);
    for i to numpol do norm_remainder[i] := norm(SNAP:-Remainder(f2[i], h0, var),2) end do;
    userinfo(1, bezout_iteration, `norm_perturbation =`, print(norm(norm_perturbation,2))); userinfo(1, bezout_iteration, `norm_remainder =`, print(norm(norm_remainder,2)));

    # data := Transpose(<u, st2, numofiteration, Normpol, norm(norm_perturbation,2), norm(norm_remainder,2)>);
    data := Transpose(<u, st2, numofiteration, norm(norm_perturbation,2)>);
    Export(data, Join(["test-", convert(d[1], string), "-", convert(d0, string), "-", excelname, ".xlsx"], ""), "Bezout", Join(["A", convert(u+1, string)], ""));
    
    return f2, V2, Normpol, numofiteration, st2, norm_perturbation, norm_remainder, B0, h0;
    # return B, J, V0, v0
end proc:

# The Jacobian Matrix with the modified Newton method
bezout_jacobian := proc (d, d0, numpol, Atemplate)
local m, J, i, j, k, l;

    # inputs: 
    # d： the degrees of the input polynomials
    # numpol: the number of the input polynomials
    # Btemplate: the template of the Bezout matrix with the input polynomials

    # output:
    # J: the Jacobian matrix with the modified Newton method





    m := d[1];
    J := Matrix((numpol-1)*m, add(d[i]+1,i=1..numpol)+m-d0);

    # the left part

    for k from 2 to numpol do
    for i from 1 to m do for j from i+1 to m-d0+1 do J((k-2)*m+i,j):= f[k,i]*v[j-1] end do end do;
    for i from 2 to m do for j from i+1 to m+1 do J((k-2)*m+i,j):= J((k-2)*m+i,j)+J((k-2)*m+i-1,j-1) end do end do;

    for j from 1 to m-d0 do for i from j to m do J((k-2)*m+i,j):= -f[k,i+1]*v[j] end do end do;
    for j from m-1 by -1 to 1 do for i from m-1 by -1 to j do J((k-2)*m+i,j):= J((k-2)*m+i,j)+J((k-2)*m+i+1,j+1) end do end do;

    for i from 1 to d0 do for j from 1 to m+1-d0 do J((k-2)*m+i+j-1,j) := J((k-2)*m+i+j-1,j)+f[k,i+m+1-d0] end do end do;
    for i from 1 to m+1-d0 do for j from 1 to d0 do J((k-2)*m+i+j-1,j+m+1-d0) := J((k-2)*m+i+j-1,j+m+1-d0)-f[k,i] end do end do;

    # the middle part

    for i from 1 to d[k]+1 do for j from i to m-d0 do J((k-2)*m+i,j+add(d[l]+1,l=1..k-1)+1):= -f[1,i]*v[j] end do end do;
    for i from 2 to d[k]+1 do for j from i to d[k] do J((k-2)*m+i,j+add(d[l]+1,l=1..k-1)+1):= J((k-2)*m+i,j+add(d[l]+1,l=1..k-1)+1)+J((k-2)*m+i-1,j+add(d[l]+1,l=1..k-1)) end do end do;


    for j from 1 to m-d0 do for i from j to m do J((k-2)*m+i,j+add(d[l]+1,l=1..k-1)):= f[1,i+1]*v[j] end do end do;
    for j from m-1 by -1 to 1 do for i from m-1 by -1 to j do J((k-2)*m+i,j+add(d[l]+1,l=1..k-1)):= J((k-2)*m+i,j+add(d[l]+1,l=1..k-1))+J((k-2)*m+i+1,j+add(d[l]+1,l=1..k-1)+1) end do end do;

    for i from 1 to d0 do for j from 1 to m+1-d0 do J((k-2)*m+i+j-1,j+add(d[l]+1,l=1..k-1)):= J((k-2)*m+i+j-1,j+add(d[l]+1,l=1..k-1))-f[1,i+m+1-d0] end do end do;
    for i from 1 to m+1-d0 do for j from 1 to d[k]-m+d0 do J((k-2)*m+i+j-1,j+add(d[l]+1,l=1..k-1)+m+1-d0):= J((k-2)*m+i+j-1,j+add(d[l]+1,l=1..k-1)+m+1-d0)+f[1,i] end do end do;

    for i from 1 to m do for j from 1 to d[k] do J((k-2)*m+i,min(j+add(d[l]+1,l=1..k),add(d[i]+1,i=1..numpol)+m-d0)) := 0 end do end do;

    end do;

    # the right part

    for i to (numpol-1)*m do for j to m-d0 do J(i,-j) := Atemplate(i,-j) end do end do:


    return J
end proc:

# The iteration function
bezout_iteration := proc (Atemplate, btemplate, J, d, d0, numpol, sigma, polvect0, X0, polvect1, X1, stopcriterion, alpha, numofiteration, num_bound, u, var, pol)
local X2, Normpol, polvect2, st, numofiteration0, sublist, f1, f2, i, B0, h0, norm_perturbation, norm_remainder, data;

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

    numofiteration0 := numofiteration + 1;
    polvect2, X2, Normpol := bezout_update(Atemplate, btemplate, J, d, d0, numpol, sigma, polvect0, X0, polvect1, X1, alpha);
    #if type(numofiteration0/20,integer) = true then userinfo(1, bezout_iteration, `numofiteration0 =`, print(numofiteration0)); userinfo(1, bezout_iteration, `Normpol =`, print(Normpol)) end if;
    userinfo(1, bezout_iteration, `numofiteration0 =`, print(numofiteration0)); userinfo(1, bezout_iteration, `Normpol =`, print(Normpol));


    if Normpol < stopcriterion or numofiteration0 > num_bound-1 then
    userinfo(1, bezout_iteration, `numofiteration0=`, print(numofiteration0));  
    userinfo(1, bezout_iteration, `Normpol =`, print(Normpol));
        return polvect2, X2, Normpol, numofiteration0;  
    else
    userinfo(2, bezout_iteration, `numofiteration0 =`, print(numofiteration0));
    return bezout_iteration(Atemplate, btemplate, J, d, d0, numpol, sigma, polvect0, X0, polvect2, X2, stopcriterion, alpha, numofiteration0, num_bound, u, var, pol); end if
end proc:

    # Find Bezout Matrix from coefficient vector

# The update function
bezout_update := proc (Atemplate, btemplate, J, d, d0, numpol, sigma, polvect0, X0, polvect1, X1, alpha)
local i, 
    j, 
    k, 
    degVV1, 
    F1poly, G1poly, 
    A1, b1,
    VV1, 
    X2,
    Q1, P1, R1, 
    dim, 
    BB, 
    vlist, 
    varlist, 
    F0grad, 
    F1grad, 
    sublist, 
    J1, Jrow, Jcol, 
    A, 
    polvect2, 
    D1, DD1, DD2, 
    VV2, V2, 
    Normpol, 
    st;
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
    # Normpol:  the distance between F1, G1 and F2, G2



    # construct substitution table for B1 & J1

    st := time();
    sublist := [seq(seq(f[i,j]=polvect1[add(d[k]+1, k=1..i-1)+j], j=1..d[i]+1),i=1..numpol), seq(seq(f[i,j]=0, j=d[i]+2..d[1]+1),i=1..numpol), seq(v[i] = X1[i], i=1..d[1]-d0)];
    st := time() - st;
    userinfo(1, bezout_update, `stsublist=`, print(st));
    userinfo(2, bezout_update, `sublist = `, print(sublist, whattype(sublist)));

    # construct Bezout Matrix of the input polynomials -> B1

    st := time();
    A1 := subs(sublist, Atemplate);
    b1 := subs(sublist, btemplate);
    userinfo(2, bezout_update, `A1 =`, print(A1));
    userinfo(2, bezout_update, `b1 =`, print(b1));
    st := time() - st;
    userinfo(1, bezout_update, `stA1b1=`, print(st));

    # construct vector of objective function

    st := time();
    VV1 := <polvect1, X1>;
    VV1 := convert(VV1, Vector);
    userinfo(2, bezout_update, `VV1 =`, print(VV1, whattype(VV1)));
    degVV1 := Dimension(VV1);
    st := time() - st;
    userinfo(1, bezout_update, `stVV1=`, print(st));

    # evaluate the constraint value -> Q1
    
    st := time();
    Q1 := A1 . X1 - b1 ;
    st := time() - st;
    userinfo(1, bezout_update, `stQ1=`, print(st));
    userinfo(2, bezout_update, `Q1 =`, print(Q1));

    st := time();
    P1 := < (polvect1 - polvect0), Vector(d[1]-d0) >;
    st := time() - st;
    userinfo(1, bezout_update, `stP1=`, print(st));
    userinfo(2, bezout_update, `P1 =`, print(P1));

    # construct right-hand-side vector -> R1

    st := time();
    R1 := - <P1, Q1>;
    st := time() - st;
    userinfo(1, bezout_update, `stR1=`, print(st));
    userinfo(2, bezout_update, `R1 =`, print(R1));



    # construct Jacobian matrix -> J1

    st := time();
    J1 := subs(sublist, J);
    st := time() - st;
    userinfo(1, bezout_update, `stJ1=`, print(st));
    userinfo(2, bezout_update, `J1 = `, print(J1));

    # construct coefficient matrix for the linear system -> A1

    st := time();
    Jrow, Jcol := Dimension(J);
    A := Matrix(Jrow+Jcol);
    for i from 1 to Jcol do A(i,i):=1 end do:
    for i from 1 to Jrow do for j from 1 to Jcol do A(i+Jcol,j):=J1(i,j) end do end do:
    for i from 1 to Jcol do for j from 1 to Jrow do A(i,j+Jcol):=-J1(j,i) end do end do:
    st := time() - st;
    userinfo(1, bezout_update, `stA=`, print(st));
    userinfo(2, bezout_update, `A =`, print(A));

    # solve the linear system: A1 . D1 == R1

    st := time();
    D1 := LinearSolve(A, R1);
    st := time() - st;
    userinfo(1, bezout_update, `stD1=`, print(st));
    userinfo(2, bezout_update, `D1 =`, print(D1));

    DD1 := D1[1..degVV1];
    userinfo(2, bezout_update, `DD1 =`, print(DD1));

    DD2 := D1[1..add(d[i]+1, i=1..numpol)];
    userinfo(2, bezout_update, `DD2 =`, print(DD2));
    
    st := time();
    VV2 := (VV1 + convert(alpha * DD1, Vector));
    st := time() - st;
    userinfo(1, bezout_update, `stVV2=`, print(st));
    userinfo(2, bezout_update, `VV2 =`, print(VV2, whattype(VV2)));

    polvect2 := VV2[1..add(d[i]+1, i=1..numpol)];
    X2 := VV2[-(d[1]-d0)..-1]; 
    Normpol := LinearAlgebra:-Norm(DD2,2);
    userinfo(2, bezout_update, `Normpol =`, print(Normpol));
    
    return polvect2, X2, Normpol;
    
end proc:

# The Bezout matrix assosiated with the polynomials

BezMatrix := proc(F, G, var)
local degF, degG, Bez, Bez1, i, j, fv, gv, gv0;

    # inputs: 

    # F:    the polynomial F
    # G:    the polynomial G

    # outputs:

    # Bez:    the Bezout Matrix of F and G

    degF := degree(F, var);
    degG := degree(G, var);  
    fv := CoefficientVector(F, var);
    gv0 := CoefficientVector(G, var);
    gv := Vector(degF+1);
    for i from 1 to degG+1 do gv(i):=gv0(i) end do;

    # let the degree of F be the larger one

    #if degF < degG then F, G := G, F ; degF, degG := degG, degF end if;

    # construct the Bezout Matrix

    Bez := Matrix(degF);

    for j from 1 to degF do 
       for i from 1 to j do
          Bez(i,j) := gv[i]*fv[j+1]-fv[i]*gv[j+1]
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


# The Bezout matrix assosiated with multiple polynomials

BezMatrix_mul := proc(pol, var)
local B, numpol, i;
    numpol :=nops (pol);
    B := BezMatrix(pol[1], pol[2], var);
    for i from 3 to numpol do B :=  <B, BezMatrix(pol[1],pol[i],var) > end do;
    return B
end proc:

# The Bezout matrix assosiated with the coefficients of polynomials

BezMatrix_2 := proc(fv, gv)
local degF, Bez, Bez1, i, j;

    # inputs: 

    # fv:    the array of the coefficients of the polynomial F
    # gv:    the array of the coefficients of the polynomial G

    # outputs:

    # Bez:    the Bezout Matrix of F and G

    degF := upperbound(fv)-1;

    # let the degree of F be the larger one

    #if degF < degG then F, G := G, F ; degF, degG := degG, degF end if;

    # construct the Bezout Matrix

    Bez := Matrix(degF);

    for j from 1 to degF do 
       for i from 1 to j do
          Bez(i,j) := gv[i]*fv[j+1]-fv[i]*gv[j+1]
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
local B0, Brow, Bcol, p, l, u, b, i, j, X, k, y, CoeffVector, gcdB;

    # inputs: 
    # B: the Bezout Matrix
    # n: the degree of the Bezout Matrix
    # h: the degree of GCD

    # outputs:

    # gcdB: the GCD calculated by the Bezout Matrix

    # extend B0 to square

    B0 := B[..,h+1..n];
    Brow, Bcol := Dimension(B0);
    B0 := <B0 | Matrix(Brow,Brow-Bcol)>;

    # LU Decompositon with B0

    p,l,u := LUDecomposition(B0);
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

remainder_poly := proc(F, G, var)
local degF, degG, VecF, A, QR, tor, VecF3, VecF4, F3, F4;

# inputs:

# F: the polynomial with higher degree
# G: the polynomial with smaller degree
# var: the variable

degF := degree(F);
degG := degree(G);
VecF := Vector(CoefficientVector(F, var, 'termorder' = 'reverse'), datatype=float);
A := Matrix(Transpose(SylvesterMatrix(G, var^(degF-degG+1), var))[..,1..degF-degG+1], datatype=float);
QR, tor := QRDecomposition(A, output='NAG');
VecF3 := LeastSquares([QR, tor], VecF);
VecF4 := VecF-convert(A.VecF3, Vector);
F3 := FromCoefficientVector(VecF3, var, 'termorder' = 'reverse');
F4 := FromCoefficientVector(VecF4, var, 'termorder' = 'reverse');
return VecF3, VecF4, F3, F4
end proc:


#remainder_vec := proc(F, G, var)
#local degF, degG, VecF, A, QR, tor, VecF3, VecF4, F3, F4;

# inputs:

# F: the coefficient vector of polynomial with higher degree
# G: the coefficient vector of polynomial with smaller degree
# var: the variable

#A := Matrix(Transpose(SylvesterMatrix(G, var^(degF-degG+1), var))[..,1..degF-degG+1], datatype=float);
#QR, tor := QRDecomposition(A, output='NAG');
#VecF3 := LeastSquares([QR, tor], VecF);
#VecF4 := F-convert(A.VecF3, Vector);
#F3 := FromCoefficientVector(Les, var, 'termorder' = 'reverse');
#F4 := FromCoefficientVector(VecF4, var, 'termorder' = 'reverse');
#return VecF3, VecF4, F3, F4
#end proc:

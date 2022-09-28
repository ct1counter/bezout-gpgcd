# Akira Terui, Boming Chi, June 2021
# The src for generating test polynomials of approximate GCD with multipol polynomials


gen_testpol_multi := proc(numpol, d1, d0, numexp, var, perturbation, Seed, maximum)
local d, i, mulcoeff, dmulcoeff, F, dF, j, k;


    # inputs:

    # numpol: the number of polynomials for each test
    # d1: the degree of the test polynomials
    # d0: the degree of approximate GCD
    # numexp: the number of the tests for each group
    # var: the variable
    # perturbation: the initial perturbation 
    # seed: the seed of the randomization
    # maximum: the maximum absolute values of the coefficients of the cofactors and approximate GCD

    with(RandomTools[MersenneTwister]):
    with(PolynomialTools):

    # Seed of the randomization

    #Seed := randomize(seed);

    # Degree of the test polynomials

    d := Vector(numpol);

    for i to numpol do d[i] := d1 end do:
    mulcoeff := seq(10*(GenerateFloat()-0.5),i=1..(add(d[j]-d0+1, j=1..numpol)+d0+1)*numexp);
    F := Array(1..numexp, 1..numpol);
    for i to numexp do 
        for k to numpol do 
            F[i,k] := simplify(FromCoefficientVector(convert([mulcoeff[add(d[j]-d0+1, j=1..k-1)+1+(i-1)*(add(d[j]-d0+1, j=1..numpol)+d0+1)..add(d[j]-d0+1, j=1..k)+(i-1)*(add(d[j]-d0+1, j=1..numpol)+d0+1)]],Vector),var)*FromCoefficientVector(convert([mulcoeff[add(d[j]-d0+1, j=1..numpol)+1+(i-1)*(add(d[j]-d0+1, j=1..numpol)+d0+1)..add(d[j]-d0+1, j=1..numpol)+d0+1+(i-1)*(add(d[j]-d0+1, j=1..numpol)+d0+1)]],Vector),var)); 
        end do; 
    end do;



    dmulcoeff := seq(GenerateFloat()-0.5,i=1..(add(d[j]+1, j=1..numpol))*numexp);
    dF := Array(1..numexp, 1..numpol);
    for i to numexp do 
        for k to numpol do 
            dF[i,k] := simplify(FromCoefficientVector(convert([dmulcoeff[add(d[j]+1, j=1..k-1)+1+(i-1)*(add(d[j]+1, j=1..numpol))..add(d[j]+1, j=1..k)+(i-1)*(add(d[j]+1, j=1..numpol))]],Vector),var)); dF[i][k] := dF[i][k]/norm(dF[i][k],2)*perturbation 
        end do; 
    end do;

    return simplify(F+dF), F, dF
    
end proc:

#save_gen_testpol_multi := proc(numpol, d1, d0, numexp, var, perturbation, Seed, maximum)
#
#    gen_testpol_multi(numpol, d1, d0, numexp, var, perturbation, Seed, maximum);

#    save Join([Join(["testpol", "mul", convert(d1, string), convert(d0, string), convert(perturbation, string)],"-")), ".mpl"], "");


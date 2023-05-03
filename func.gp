/* Various routines with rational functions,
   mostly with Laurent polynomials */

\\ We will do Taylor expansion around w0=1.
Taylor(X,i,w,w0=0) = polcoef(subst(X,w,w0+t),i,t);

\\ Reverse the Laurent polynomial, or other function of z
rev(F,z=z) = subst(F,z,1/z);

\\ Check that function is a Laurent polynomial
islaurent(F,z=z) =
{
nF = numerator(F,z);
dnF = poldegree(nF,z);
poldegree(rev(nF)*z^dnF,z)==0;
}

\\ ldeg returns -1 if F is NOT Laurent, otherwise returns its degree ≥0.
ldeg(F,z=z) = if(!islaurent(F,z),-1,max(poldegree(F,z),poldegree(rev(F),z)));

\\ grading(F,v) returns [integer d,boolean b] where b = ishomogeneous(F), d = grading(F)
grading(F,v=[a,b]) = [dd=poldegree(F=substvec(F,v,hh*v),hh),poldegree(hh^dd*rev(F,hh),hh)==0];
/* unused */ ishomogeneous(F,v=[a,b]) = grading(F,v)[2];

T0(F,z=z) = polcoef(F,0,z);
T(N,F,z=z,de=ldeg(F)) =  sum(k=ceil(-de/N),de/N,polcoef(F,N*k,z)*z^k);

\\ Create an empty list to put all the respective invariants
functions = List();
names = List();

/*
    Function verify(F) verifies that F
    is a Laurent polynomial, gives its degree,
    then specifies for which N ≤ degree
    T_N(F) equals to T_0(F).
*/
verify(F,name="F",z=z) = 
{
listput(~functions,F);
listput(~names,name);
if((de=ldeg(F,z))<0,print(name," is NOT a Laurent polynomial");0, print("(T_N-T_0) ",name," = 0 for N in ",vecsort(vector(de,n,n*(T(n,F,z,de)==T0(F,z))),,8)," and N>",de,".");1);
}

\\ myfac(G) prints factorization of G
myfac(G,name="F",z=z) = print(name," = ",G,", and it factorizes as ",Fac=factor(G)," * ",G/factorback(Fac));
\\," where w=exp(2πi k/N).");

check(F,name="F",z=z) = if(verify(F,name,z),myfac(T0(F),strjoin(["T_0 ",name]),z));

/* 
   func.gp : routines with functions
   geom.gp : basic geometry - parametrization and basic functions (area, scalar product, cot, etc)
*/
\r func.gp
\r geom.gp

check(1,"1 (constant)");

\\ 4.1 areas
check(area(S*P),"{5} A");
check(area(S*(P-Prev(P))), "Area(P(i-1)P(i)P(i+1)");

\\ distance to focus

check(xy2z(a+cc*x/a), "{13} d1");

\\ 4.2 scalar products [4] 
print("(8:C2), that implies (6:L1)");
check( Cor2 = scalar(S*(P-[x0,y0]~)), "{8} |P-[x0,y0]|^2" );
print("(7:P4)");
check( scalar(S*(P-Next(P))), "{7} |P(i+1)-P(i)|^2" );

\\ 4.3 cotangent [5]
print("(9:C3T1)");
check( The1 =  -Cot(S*(P-Prev(P))),"{9} Cot θ" );
print("(10:O1a)");
check( Obs1 = Cot(S*(P-Prev(P)))^2,"{10} Cot(θ)^2");
eps = (w-1)/I; \\ eps = φ + O(φ^2) for φ small. If φ=2π/N, then eps = 2π/N + O(1/N^2)
print("{11a:O1b} The limit of Cot*(2π/N) is ",subst(T0(The1)*eps,w,1));
print("{11b:01c} The limit of (Cot*(2π/N))^2 is ",subst(T0(Obs1)*eps^2,w,1));

if(T0(Cot(S*(P-Prev(P))) - Cot(S*P)) == 0, print("{NEW} Observation 1+ε: T_0 Cot θ = T_0 Cot φ"));
if(T0(Cot(S*(P-Prev(P)))^2 - Cot(S*P)^2) == 0, print("{NEW} Observation 1+2ε: T_0 (Cot θ)^2 = T_0 (Cot φ)^2"));

\\ 4.4 curvature [6]
if(prod(n=0,2,polcoef(QQ,n,t)==0),print("Osculation verified."),print("ERROR (osculation)"));

if(ispower(R2/a/b,3,&kapa),print("(14:C4)");check(kapa,"{14} (ab/κ^2)^(1/3)"));

print("(ADD:C4+ε): z coordinate of the fourth intersection point is ",zFourthPoint);
print("(ADD:C4+2ε): the coordinates of the center of the osculating circle are:"); 
myfac(Xc,"X_c");
myfac(Yc,"Y_c");

\\ 4.5 evolute [7]

print("(15:P6)");
verify(As ,"{15} As");
T0As=T0(As);

print("(C3:C5) The degree of A_s with respect to s equals ",zAs = poldegree(As,s),", thus for every w there are ",zAs," values of s such that A_s vanish");

print("(O1:O2)");
factor(poldisc(T0As,s))~
\\ then see that the last term is a square iff
\\ w is 3rd or 6th root of 1, or for 4 real values of w (not on unit circle).
factor(substpol(substpol(T0As,w^3,1),w^2,-1-w))~

/*
print("Observation 3  (cannot verify)");
CotPs = Cot(Ps);
verify(CotPs,"Cot(Ps)");
*/

print("(16:O4): lim As/2π = ",subst(T0As/eps,w,1));

print("(17:P7)");
verify( ell2 = scalar(Ps-Next(Ps)), "{17} ell^2" );
print("{17} T_0 ell2 = ",T0(ell2));

names = Vec(names);
functions = Vec(functions);

print("Now we have collected ",#names," polynomial quantities, namely ",strjoin(names,", "));
print("Their degrees are ", degs = apply(poldegree,functions));
grau = vecmax(degs);
print("Moreover, they are linearly dependent. For example, we have linear dependence ",lindep(Vec(functions)));

FF = matrix(#functions, 2*grau+1, i, j, polcoef(functions[i],j-1-grau,z));
print("We have the following linear dependencies: ",matker(FF~)~);

print("In particular, we see a linear dependence between 1, cot(θ) and κ^-(2/3) ");

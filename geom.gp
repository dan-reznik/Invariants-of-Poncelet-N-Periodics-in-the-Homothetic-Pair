/* The parametrization of points P = [x(z),y(z)] on the circle by complex numbers z */
P = [(z+1/z)/2,(z-1/z)/2/I]~;
xy2z(q,V=[x,y]) = substvec(q, V, P);
/* unused */ z2xy(q,z=z,V=[x,y]~) = subst(q,z,[1,I]*V);
/* gives an isomorphism of the circle rotations group and the multiplicative group. */
Next(q) = subst(q,z,z*w);
Prev(q) = subst(q,z,z/w);
/* The affine transformation                                                        */
S = [a,0;0,b];
Scale(Q,T=S,V=[x,y]~) = substvec(Q,V,T*V);
/* Wedge product - oriented area of the triangle <0, V, V2>, default V2 is next V   */
area(V,V2=Next(V)) = matdet(Mat([V,V2]))/2;
/* Scalar product of V and V2 (by default of V with itself)                         */
scalar(V,V2=V) = V2~ * V;
/* cotangent of an angle                                                            */
Cot(V,V2=Next(V)) = scalar(V, V2) / (2*area(V,V2)) ;
/* Below [Xc,Yc] is the center of the osculating circle, and R2 is its radius^2     */ 
/* (i.e. 1/curvature^2), functions of z.                                            */
/* Circle is the equation of the circle with center [xc,yc] and R^2 = xc^2+yc^2-zc  */
/* Q is S^-1(Circle) restricted to UnitCircle, in z shifted by t.                   */
Circle = x*(x-2*xc)+y*(y-2*yc)+zc;
Q = subst(numerator(xy2z(Scale(Circle))),z,z+t);
/* Next 4 lines compute [Xc,Yc] and R2                                              */
M =   matrix(3,3,i,j,  polcoef(polcoef(Q,i-1,t),1,[xc,yc,zc][j]));
B = -vectorv(3,  i,   substvec(polcoef(Q,i-1,t),  [xc,yc,zc],[0,0,0]));
[Xc,Yc,Zc] = matsolve(M,B);
R2 = Xc^2 + Yc^2 - Zc;
/* Let us verify the tangency and see the fourth intersection point                 */
QQ = substvec(Q,[xc,yc,zc],[Xc,Yc,Zc]);
zFourthPoint = z - polcoef(QQ,3,t)/polcoef(QQ,4,t);
FourthPoint = subst(P,z,zFourthPoint);
/* A linear combination of the center [C,D] and the original point P, and its area. */
As = area(Ps = (1-s)*S*P + s*[Xc,Yc]~);

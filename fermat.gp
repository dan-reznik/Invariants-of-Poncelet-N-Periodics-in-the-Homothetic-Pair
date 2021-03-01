/* 
https://faculty.evansville.edu/ck6/encyclopedia/ETC.html
X(13) = 1st isogonic center (Fermat point, Torricelli point).
It minimizes the sum of distances to the vertices of given a triangle ABC,
|AX|+|BX|+|CX| is minimum when X=X13
*/
X13(A,B,C) = {
/* Let a,b,c be squares of triangle lenght. a = |BC|ˆ2 &c 
       a,b,c - symbolic, aa,bb,cc - numeric.
       Aa,Bb,Ca - the 3 respective vectors
       I use formula for oriented area, 
       so just choose the orientation of the triangle properly
*/ 
\\\ Det = matdet(...)		\\ area = det/2
\\ Q = 2*sqrt(3);
k = Q*Det;
w = vector(3);
w[1] = a^2 - 2 * (b - c)^2 + a * (b + c + k);
w[2] = substvec(w[1],[a,b,c],[b,c,a]);
w[3] = substvec(w[1],[a,b,c],[c,a,b]);
Aa = C-B; aa = Aa * Aa~;
Bb = A-C; bb = Bb * Bb~;
Cc = B-A; cc = Cc * Cc~;
w = substvec(w,[a,b,c],[aa,bb,cc]);
P = w*[A,B,C]~ / vecsum(w);
return(P);
}

-- checking zero divisors in M2

A = ZZ[a,b,c];
I = ideal((a+b)*(b+c)); 
h = a^2
matrix{{h}}**(A/I)
(ker(matrix{{h}}**(A/I))) 

h = a+b
matrix{{h}}**(A/I)
(ker(matrix{{h}}**(A/I))) 
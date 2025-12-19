# matrixEquationSolver
This program is designed to solve the linear equation Ax = b.  
Let me introduce the algorithm: we compute the reduced row echelon form of augmented matrix [A|b], 
then it is obvious to see whether the equation can be solved or not by comparing the rank of A and [A|b] (counting pivots to solve the rank).  
If it can, we let all free variables to be 0, then we have the Xp(X particular);  
Then we let a free variable to be 1 while others are 0, thus construct the null space of A;  
the column vectors of N(A) are the X(n)'s  
The solution to Ax = b is X = Xp + ci * Xni

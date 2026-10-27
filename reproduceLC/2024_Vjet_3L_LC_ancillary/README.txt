**************************************************************************************************************************
*  Planar three-loop QCD helicity amplitudes for $V$+jet production at hadron colliders                                  *
*  Thomas Gehrmann, Petr Jakubčík, Cesare Carlo Mella, Nikolaos Syrrakos, Lorenzo Tancredi                               *
**************************************************************************************************************************

The ancillary files to this article include:
 - projectors onto the helicity amplitude coefficients
 - helicity amplitude coefficients after renormalization but before IR subtraction in the decay region
 - finite remainders of the helicity amplitudes in decay kinematics and continued to all production regions

The files are in Mathematica-readable format. 

1) Projectors
The file Projectors.m contains a matrix where the coefficient of the j-th tensor structure in the projector onto the i-th helicity 
amplitude coefficient is given by the element in the i-th row and j-th column. Helicity amplitude coefficients are in this order:
(alpha_1, alpha_2, alpha_3, gamma_1, gamma_2, gamma_3)

2) Renormalized helicity amplitude coefficients
The file RenormalizedAmpCoeffs.m contains the six helicity amplitude coefficients after UV renormalization but before
IR subtraction in Minkowskian decay kinematics. The file is a list of replacement rules in the format
        coefficient[index, loops] -> expression,
eg.     alpha[1, 3] -> 0
is the three-loop expression for the helicity coefficient alpha_1. 

3) Finite remainders and their analytic continuation
In the directory FiniteRemainders, we provide the finite remaiders of the independent helicity amplitudes in the decay region, 
and continued to all relevant production regions as explained in Section 6.3 of 2306.10170. The structure is as follows:

        - decay.m             V_L -> qb_R(p1) + q_L(p2) + g_+(p3)  AND   V_L -> qb_R(p1) + q_L(p2) + g_-(p3) 
        - lmL.m               qb_L(p1) + g_-(p3) -> V_R + qb_L(p2)     
        - lpL.m               qb_L(p1) + g_+(p3) -> V_R + qb_L(p2) 
        - rmR.m               q_R(p1) + g_-(p3)  -> V_R + q_R(p2) 
        - rpR.m               q_R(p1) + g_+(p3)  -> V_R + q_R(p2) 
        - lrP.m               qb_L(p1) + q_R(p3) -> V_R + g_+(p2) 
        - lrM.m               qb_L(p1) + q_R(p3) -> V_R + g_-(p2) 

Each file is a list of replacement rules in the format
        coefficient[index, loops] -> expression,
eg.     alpha[1, 3] -> 0
is the three-loop expression for the helicity coefficient alpha_1. 
Some files contain only the coefficients alpha, others only gamma, based on which of the two independent helicity amplitudes
from the decay region was continued. Expressions in the production regions are in the variables v = 1/y and u = -z/y,
where y = s13, z = s23. We denote M_V^2 with q2.
   
   
   
   
   
   
   
   
  




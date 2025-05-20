# Parallel-Sequence-Alignment

Project by Ángela Garcinuño Feliciano and Lucas Massot

Interesting links:
  - https://en.wikipedia.org/wiki/Edit_distance
  - https://dl.acm.org/doi/pdf/10.1145/2893488?casa_token=V8yiw_twDMwAAAAA:EVi4esMB8pIzQawnjLJVF9HPxXBe1FIbJQLbtktXFVJpZQVV6nQVVbA6NcLLQRdpRsfe1gOf279X7Rs


Deadlines:
  - 4th May: Read paper and choose algorithm (Unofficial)
  - 21th May: Initial report (Official)
  - 25th May: Finish Code (Official)
  - 1st June: Finish Report (Unofficial) (Deadline to submit project for review)
  - 7th June: Project Deadline (Official)
  - 12th June: Project Defense (Offficial)



## Algorithm Description
1. Sequence Alignment Problem
   - Input: A, B sequences of lengths n and m respectively
   - Output: Optimal alignment and score
  
   - Definitions:
       - Extended_Problem(A, B, start_type, end_type): takes two subsequences A and B as well as the start and end alignment types (i.e. starts/ends with alignment to gap or sequence).
       - Partial Balanced Partition of A and B: a list $P = ((i_0, j_0, t_0), (i_1, j_1, t_1), ..., (i_{|P|-1}, j_{|P|-1}, t_{|P|-1})$. It is the list of positions where an alignment can be split into subproblems. (Propositions and properties found in article). A partial balanced partition contains enough information to solve the problem. Algorithm to find such a partition given in section 5 of the article.

   - Problem: Align(A,B). This extended_Problem allows us to divide the original problem into subproblems and then combine them to find the optimal solution for the original problem. 
       - Base case: solve alignment problem for small enough sequence.
       - Determine split point: we need to split the algorithm at a point which we know will be in the optimal alignment.
       - Recursively solve each subproblem.
       - Join solutions
    
### Find Paritial Balanced Partition:
We use 3 dynamic programming tables: T1, T2 and T3. Let $A = a_{i'}, a_{i'+1}, ..., a_{i''}$ and $B = b_{j'}, b_{j'+1}, ..., b_{j''}$. We want to solve the subproblem (A,B,s,e). We will only use cells [i,j] such that $i'-1 \leq i \leq i''$ and $j'-1 \leq j \leq j''$ in each table. The value of entry [i,j] is the score for optimally aligning $a_{i'}, a_{i'+1}, ..., a_{i}$ with $b_{j'}, b_{j'+1}, ..., b_{j}$ with the following conditions:
 - In T1, $a_i$ must be matched with $b_j$.
 - In T2, $b_j$ is matched to a gap.
 - In T3, $a_i$ is matched to a gap.
Then, $T_k[i,j]$ contains the score of an optimal alignment to the problem ($a_{i'}, a_{i'+1}, ..., a_{i}$, $b_{j'}, b_{j'+1}, ..., b_{j}$, s, k)

#### Algorithm:
##### Compute T1, T2 and T3
Initialisation:
  - [i' -1, j' -1] =
    - If start_type <= 1: T|start_type = 0
    - otherwise: - \infty
  - if start_type > 0:
      - start_type == 2: T2[i'-1, j'] = -(g+h)
      - start_type != 2: T2[i'-1, j'] = - \infty
      - start_type == 3: T3[i'-1, j'] = -(g+h)
      - start_type != 3: T3[i'-1, j'] = - \infty

Compute remaining cells row by row:
- T1[i,j] = f(ai, bj) + max(T1[i-1, j-1], T2[i-1, j-1], T3[i-1, J-1])
- T2[i,j] = max(T1[i, j-1] - (g+h), T2[i, j-1] - g, T3[i, J-1] - (g+h))
- T3[i,j] = max(T1[i-1, j] - (g+h), T2[i-1, j] - (g+h), T3[i-1, J] - g)
- Top row and leftmost column have predetermined values (see article page 7).

Once tables are filled:
- h'(k) = h if k = end_type and end_type in {-2, -3} ; 0 otherwise.
- If end_type > 0: optimal score = T|end_type[i'',j'']
- If end_type < 0: optimal score = max(T1[i'',j''], T2[i'',j''] + h'(-2), T3[i'',j''] +h'(3)]
- Alignment extracted using an origin traceback procefure starting from the entry said to contain the optimal score.

Note: value of each cell depends only on left, upper-left and upper neighbours.

We define the *orign* of a cell [i,j]k to be the cell from which Tk[i,j] was calculated. [i'-1, j'-1] cells and cells whose values are - \infty are without origin.

##### Find elements of a partial balanced partition
We use tables T1, T2 and T3 to find cells where we can perform a recursive decomposition of the problem.

We want to find cells thar lie on the solution.

Let $rev(A) = a_{i''}, a_{i''-1}, ..., a_{i'}$ and $rev(B) = b_{j''}, b_{j''-1}, ..., b_{j'}$. Let $T_k^R[i,j]$ denote the score of an optimal alignment of $a_{i''}, a_{i''-1}, ..., a_{i+1}$ and $b_{j''}, b_{j''-1}, ..., b_{j+1}$ under the extended alignment problem with start_type = e and end_type = k. The process of computing the tables $T_1^R$,  $T_2^R$ and $T_3^R$ is similar to computing T1, T2, T3 but the seed cell is [i'',j''], the optimal score is in [i'-1, j'-1] and the computation and it's rules are reversed.

Consider any cell [i,j] of tables T1, T2 and T3. Let $A_1 = a_{i'}, a_{i'+1}, ..., a_{i}$, $B_1 = b_{j'}, b_{j'+1}, ..., b_{j}$, $A_2 = a_{i+1}, a_{i+2}, ..., a_{i''}$, $B_2 = b_{j+1}, b_{j+2}, ..., b_{j''}$. The best possible alignment of A and B that passes through [i,j] has score:

$$opt(i,j) = max(T_{max}[i,j] + T^R_{max'}[i,j], T_2[i,j] + T^R_2[i,j] + h, T_3[i,j] + T^R_3[i,j] + h)$$

Where $T_{max} = max(T_1[i,j], T_2[i,j], T3[i,j])$ (same for $T^R_{max'}$). A solution passes through a cell if opt(i,j) is equal to the score of an optimal alignment.

#### Alg
We assume we are given p processors, with the sequences A' and B' distributed such that $P_k$ is given $b_{k\frac{n}{p} +1},...,b_{(k+1)\frac{n}{p}}$ and $a_{k\frac{m}{p} +1},...,a_{(k+1)\frac{m}{p}}$.

[Stopped here - end of page 7]


Computation of Recombination Probabilities
------------------------------------------

As mentioned earlier, probabilities of recombination events, denoted
$g_{{{{{\epsilon_i}}}}}$, play a key role in
linkage analysis (e.g., AHGL pp 117-120). Here we describe the
relationship between these recombination probabilities and map distances
between loci. Using this relationship, all recombination probabilities
can be computed given the map distances between loci.

To establish this relationship for $k$ loci, we let
$W_{{{{{\delta_i}}}}}$ denote a region of the
chromosome composed of inter-locus segments. Element $j$ of the
$k-1\times1$ vector  is $1$ if the segment between loci $j$ and $j+1$ is
included in $W_{{{{{\delta_i}}}}}$ and is $0$
otherwise. The length of the region
$W_{{{{{\delta_i}}}}}$ is

$$x({{{{\delta_i}}}})=\sum{{{{\delta_i}}}}_{j}x_{j}$$

where $x_{j}$ is the map distance between loci $j$ and $j+1$.

The probability of an odd number of crossovers occurring in
$W_{{{{{\delta_i}}}}}$ is denoted
$R({{{{\delta_i}}}})$ and is called the
recombination value for . Given a map function $r_{x}=M(x)$, the
recombination value for  can be computed as

$$R({\delta_i})=M[x({\delta_i})]$$

The recombination value for  can also be computed as the sum of those
$g_{\epsilon_i}$ for which there is an
odd number of recombinations in
$W_{\delta_i}$ (This rule works because
the sum of an odd number of odd numbers is odd; the sum of an even
number of odd numbers is even; and the sum of even numbers is always
even.). For example, if $k=4$ and

${{{{\delta_i}}}}=[1,0,1]'$,
$R({{{{\delta_i}}}})=g_{001}+g_{011}+g_{100}+g_{110}$.

The number $s_{ij}$ of recombinations in the region
$W_{{{{{\delta_i}}}}}$ given recombination
event  is

$$s_{ij}={{{{\delta_i}}}}'{{{{\epsilon_j}}}}$$

So, recombination value for  can be written as

$$\begin{split}
R({\delta_i})
& =\sum_{j\,{for {s_{ij}}odd}}g_{{{{{\epsilon_j}}}}}\\
& ={{\frac{1}{2}}}[1-\sum_{j=1}^{2^{k-1}}(-1)^{s_{ij}}g_{{{{{\epsilon_j}}}}}]
\end{split}
$$

In matrix notation, the above relationship between the
$R({\delta_i})$’s and the
$g_{ \epsilon_j}$’s can be written as

$$r={\frac{1}{2}}[1-{Ag}]$$

where is a $2^{k-1}\times1$ vector of recombination values, is a
$2^{k-1}\times1$ vector of 1’s, the matrix
${A}=\{(-1)^{s_{ij}}\}$, and is a
$2^{k-1}\times1$ vector of recombination probabilities. Rearranging
([Matrix-R~d~elta]) gives
$${Ag}={1}-2{r}$$
The following properties can be shown to be true for the matrix :
${{{A}}}={{{A}}}'$,
${{{a}}}_{i}'{{{a}}}_{i}=2^{n-1}$
for $i=1,\ldots,2^{k-1}$, and
${{{ a}}}_{i}'{{{a}}}_{j}=0$
for $i\neq j$. So,
${{{AA}}}={{{I}}}2^{k-1}$
and
${A}^{-1}={{{A}}}\frac{1}{2^{k-1}}$.
Now, $g$ can be written in terms of $r$ as

$$\begin{split}
g & = A^{-1}(1-2r)\\
& =\frac{A(1-2r)}{2^{k-1}}
\end{split}
$$

In scalar notation, the recombination probability for  can be written as

$$g_{\epsilon_i}=\sum_{j}^{2^{k-1}}(-1)^{s_{ij}}\frac{1-2R(\delta_{j})}{2^{k-1}}
$$

Equations ([R~d~elta1]) and ([g~e~psilon]) establish the relationship
between map distances and recombination probabilities.

Consider the numerical example in AHGL (pp 125–126). Here, $k=4$ and the
recombination rates in the three intervals and the corresponding map
distances for the Haldane and the binomial ($N=2$) map functions are
given in table ([lengths-table]).

[lengths-table]

|Segment| $r_{j}$ | $x_{j}$|$x_{j}$|
|---------|---------|--------|
|$j$||Haldane|Binomial|
|1 | 0.1 | 0.1116 | 0.1056|
|2 | 0.05| 0.0527 | 0.0513|
|3 | 0.2 | 0.2554 | 0.2254|


The set of vectors $\epsilon_i$ and $\delta_{j}$ are
given by the rows of the matrix :
<br><br>

$${U}=\begin{bmatrix}0 & 0 & 0\\
0 & 0 & 1\\
0 & 1 & 0\\
0 & 1 & 1\\
1 & 0 & 0\\
1 & 0 & 1\\
1 & 1 & 0\\
1 & 1 & 1
\end{bmatrix}$$

Thus, the matrix of $s_{ij}$’s is
$${S}=UU'=
\begin{bmatrix}0 & 0 & 0 & 0 & 0 & 0 & 0 & 0\\
0 & 1 & 0 & 1 & 0 & 1 & 0 & 1\\
0 & 0 & 1 & 1 & 0 & 0 & 1 & 1\\
0 & 1 & 1 & 2 & 0 & 1 & 1 & 2\\
0 & 0 & 0 & 0 & 1 & 1 & 1 & 1\\
0 & 1 & 0 & 1 & 1 & 2 & 1 & 2\\
0 & 0 & 1 & 1 & 1 & 1 & 2 & 2\\
0 & 1 & 1 & 2 & 1 & 2 & 2 & 3
\end{bmatrix}$$

and the matrix is

$${A}=(-1)^{s_{ij}}=
\begin{bmatrix}
1 & 1 & 1 & 1 & 1 & 1 & 1 & 1\\
1 & -1 & 1 & -1 & 1 & -1 & 1 & -1\\
1 & 1 & -1 & -1 & 1 & 1 & -1 & -1\\
1 & -1 & -1 & 1 & 1 & -1 & -1 & 1\\
1 & 1 & 1 & 1 & -1 & -1 & -1 & -1\\
1 & -1 & 1 & -1 & -1 & 1 & -1 & 1\\
1 & 1 & -1 & -1 & -1 & -1 & 1 & 1\\
1 & -1 & -1 & 1 & -1 & 1 & 1 & -1
\end{bmatrix}$$

The map lengths for $W_{\delta_i}$ and
corresponding recombination values computed from Haldane and binomial
($N=2$) map functions are in table ([r-probs-table]). To obtain the
recombination probabilities using the Haldane and binomial map
functions, vectors from the fifth and seventh columns of table
([r-probs-table]) were used with the matrix from ([A-matrix]) in
equation ([Matrix-g~e~psilon]). These probabilities are given in table
([g~e~psilon-table]).

The recombination probabilities under the Haldane map function can be
computed much more simply due to the lack of interference. For example,
$$g_{110}=r_{1}r_{2}(1-r_{3})=(0.1)(0.05)(1-0.2)=0.004$$.

However, this approach cannot be used when interference is present. Note that
the Kosambi map function cannot be used for mapping more than 3 loci. When the
Kosambi map function was used for this example, $g_{111}$ was negative (AHGL p
126).

Table 2: Map length $x(\delta_i)$ and recombination value $R(\delta_i)$ for each
$\delta_i$, computed from Haldane and binomial (N = 2) map functions.

|$\delta_i$| | | Haldane|       |Binomial|        |
|:----:|:----:|:----:|:--------:|:-------:|:--------:|:--------:|
|    |    |    |$x(\delta_i)$|$R(\delta_i)$|$x(\delta_i)$|$R(\delta_i)$|
| 0  | 0  | 0  |  0     | 0     | 0      | 0      |
| 0  | 0  | 1  | 0.2554 | 0.2   | 0.2254 | 0.2    |
| 0  | 1  | 0  | 0.0527 | 0.05  | 0.0513 | 0.05   |
| 0  | 1  | 1  | 0.3081 | 0.23  | 0.2767 | 0.2384 |
| 1  | 0  | 0  | 0.1116 | 0.1   | 0.1056 | 0.1    |
| 1  | 0  | 1  | 0.3670 | 0.26  | 0.3310 | 0.2762 |
| 1  | 1  | 0  | 0.1642 | 0.14  | 0.1569 | 0.1446 |
| 1  | 1  | 1  | 0.4197 | 0.284 | 0.3823 | 0.3092 |

Table 3: Probabilities of recombination events (g\epsilon_i)  computed from
Haldane Binomial
and binomial (N = 2) map functions.

|$\epsilon_i$| | |g\epsilon_i||
|:----:|:----:|:----:|:--------:|:-------:|:--------:|:--------:|
|  |   |   |Haldane| Binomial|
|0 | 0 | 0 | 0.684 | 0.6704|
|0 | 0 | 1 | 0.171 | 0.1823|
|0 | 1 | 0 | 0.036 | 0.0415|
|0 | 1 | 1 | 0.009 | 0.0058|
|1 | 0 | 0 | 0.076 | 0.0854|
|1 | 0 | 1 | 0.019 | 0.0119|
|1 | 1 | 0 | 0.004 | 0.0027|
|1 | 1 | 1 | 0.001 | 0     |


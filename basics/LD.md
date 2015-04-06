
#1. Linkage Disequilibrium in an Infinite Population


Gametic disequilibrium, which is more commonly referred to as linkage
disequilirium (LD), is the statistical dependence between alleles in a
haplotype. Under gametic equilibrium, alleles in a haplotype are
independent. This is also called linkage equilibrium. Note that it is
possible for two loci that are linked to be in gametic equilibrium;
also, loci that are unlinked can be in gametic disequilibrium.

Suppose that starting from generations 0, all individuals are produced
by random mating. Then, the probability of haplotype ($A_{i},B_{j}$) in
generation 1 can be written as

$${\textstyle
\Pr_{1}(A_{i},B_{j})=(1-r)\Pr_{0}(A_{i},B_{j})+r\Pr(A_{i})\Pr(B_{j})}$$

where $r$ is the recombination rate between loci $A$ and $B$, and
$\Pr_{0}(A_{i},B_{j})$ is the probability of ($A_{i},B_{j}$) in
generation 0. The disequilibrium in generation 1 is

$$\begin{split}
\Delta_{1} & ={\Pr}_{1}(A_{i},B_{j})-\Pr(A_{i})\Pr(B_{j})\\
 & =(1-r){Pr}_{0}(A_{i},B_{j})+r\Pr(A_{i})\Pr(B_{j})-\Pr(A_{i})\Pr(B_{j})\\
 & =(1-r){Pr}_{0}(A_{i},B_{j})-(1-r)\Pr(A_{i})\Pr(B_{j})\\
 & =(1-r)\Delta_{0}
\end{split}$$

where $\Delta_{0}$ is the disequilibrium in generation 0.
Similarly, the probability of haplotype ($A_{i},B_{j}$) in generation 2
is

$${\textstyle
\Pr_{2}(A_{i},B_{j})=(1-r)\Pr_{1}(A_{i},B_{j})+r\Pr(A_{i})\Pr(B_{j})}$$

and the disequilibrium in generation 2 is

$$\begin{split}
\Delta_{2} & ={\Pr}_{2}(A_{i},B_{j})-\Pr(A_{i})\Pr(B_{j})\\
 &
=(1-r){\mbox{Pr}}_{1}(A_{i},B_{j})+r\Pr(A_{i})\Pr(B_{j})-\Pr(A_{i})\Pr(B_{j})\\
 & =(1-r){\mbox{Pr}}_{1}(A_{i},B_{j})-(1-r)\Pr(A_{i})\Pr(B_{j})\\
 & =(1-r)\Delta_{1}\\
 & =(1-r)^{2}\Delta_{0}
\end{split}$$

It follows that in generation $n$, the disequilibrium is

<br><br>
$$\Delta_{n}=(1-r)^{n}\Delta_{0}$$

Thus, with each generation of random
mating the haplotype distribution moves closer to equilibrium
(statistical independence). For loci that are unlinked, $(1-r)=1/2$ and
equilibrium is reached quickly; for example, $(1/2)^{10}=1/1024$. On the
other hand, loci that are tightly linked will take much longer to reach
equilibrium; for example, $(1-r)^{10}>1/3$ for $r=0.1$. In the limit,
however, equilibrium is reached in an infinite population.

#2. Linkage Disequilibrium in a Finite Population


In a closed finite population, in the absence of mutation, after a
sufficient number of generations of random mating all alleles will
become identical by descent; thus all alleles also will become identical
in state. In other words, all loci will become “fixed”. In such a
population genetic variability is absent and LD is not defined.

Loci that are moving toward fixation will pass through a phase where
alleles that are identical by state (IBS) will also be identical by
descent (IBD). In other words, all alleles that are IBS will trace back
to a common ancestral mutant allele. Thus, at a biallelic locus A with
alleles $A_{1}$ and $A_{2}$, all the $A_{1}$ alleles will have a common
ancestor and all the $A_{2}$ alleles will have a different common
ancestor. In such loci, Sved @sved1971linkage has shown algebraically
that the expected value of LD between loci A and B as measured by the
squared correlation ($\rho^{2}$) between allele states is related to the
probability of joint identity-by-descent at loci A and B for two
randomly sampled gametes. At first, this relationship may not be
obvious. To get an intuitive feel for this relationship, consider a
population where alleles $A_{1}$ and $A_{2}$ are segregating at the A
locus and alleles $B_{1}$ and $B_{2}$ at the B locus. Suppose all
haplotypes with $A_{1}$ also have $B_{1}$ and those with $A_{2}$ have
$B_{2}$. Then, $\rho$ would be 1. If the allelic associations were
reversed, i.e., $A_{1}$ goes with $B_{2}$ and $A_{2}$ goes with $B_{1}$,
then $\rho$ would be -1. In both cases $\rho^{2}$ is 1. In such a
population for two randomly sampled gametes, the conditional probability
would be one that alleles at the B locus are identical-by-state (IBS)
given they are IBS at the A locus. On the other hand, if there are a few
haplotypes where the association between alleles is different from the
other haplotypes, then $\rho^{2}$ would be less than 1. Further, in this
population for two randomly sampled gametes, the conditional probability
would be less than 1 that alleles at the B locus are identical-by-state
(IBS) given they are IBS at the A locus. Hopefully, this discussion
helps to see that LD as measured by $\rho^{2}$ should be related to the
probability of joint IBS. Recall, however, that we assumed that all
$A_{1}$ alleles have descended from a common ancestor and similarly all
$A_{2}$ alleles have descended from a different common ancestor. Thus,
given alleles at the A locus are IBS, they have descended from a common
ancestor (IBD), and provided that no recombination between the two loci
has happened in the two paths descending from the common ancestor to the
two randomly sampled haplotypes, the B alleles will also be IBD. Sved
@sved1971linkage denoted this conditional probability by $Q$, and
reasoned that in the pool of gametes where alleles at the B locus are
IBD given they are IBD at the A locus, LD as measured by the squared
correlation ($\rho^{2}$) between allele states will be 1.0
@sved2009linkage. On the other hand, in the pool of gametes where
recombination has taken place between loci A and B, given random mating,
$\rho^{2}$ is expected to be null @sved2009linkage.

Let $C=1$ denote the condition that for a randomly sampled pair of
gametes the alleles at locus B are IBD given that alleles are IBD at
locus A, and $C=0$ denote that this condition is not met. Then, the
expected value of $\rho^{2}$ can be written as

(1)
$$\begin{aligned}
E(\rho^{2})
& = \underset{C}{E}[E(\rho^{2}|C)]\\
& = E(\rho^{2}|C=1)\Pr(C=1)+E(\rho^{2}|C=0)\Pr(C=0)\\
& = 1Q+0(1-Q)\\
& = Q.
\end{aligned}$$

Let $Q_{t}$ denote the probability of $C=1$ in generation $t$. This
probability can be recursively written as

(3)
$$Q_{t}=[\frac{1}{2N}+(1-\frac{1}{2N})Q_{t-1}](1-r)^{2}$$

where $N$ is the effective population size, $\frac{1}{2N}$ is the
probability that two randomly sampled gametes in generation $t$ are both
inherited from the same gamete in the previous generation,
(1-$\frac{1}{2N}$) is the probability that they are inherited from
different gametes in the previous generation, and $(1-r)^{2}$ and the
probability that loci A and B do not recombine in these gametes in the
last generation. At equilibrium, $Q_{t}=Q_{t-1}=Q_{E}$. So, setting
$Q_{t}$ and $Q_{t-1}$ in to $Q_{E}$ gives

(3)
$$\begin{split}
Q_{E} & =\frac{(1-r)^{2}}{2N-(2N-1)(1-r)^{2}}\\
 & \approx\frac{1}{4Nr+1}
\end{split}$$

In the derivation of (3) we only considered pairs of loci where
alleles are segregating at each locus. Further, it was assumed that all
haplotypes in the current generation descend from two ancestral
haplotypes. So, if we code the four possible haplotypes at a pair of
biallelic loci as: 00, 01, 10, an 11, the ancestral pair of haplotypes
must be either (00,11) or (01,10). Any other pair would lead to one
locus being fixed. For example, the pair (00,01) has locus one fixed for
the allele coded as 0. So, if $4Nr$ is close to zero, most haplotypes in
the current generations will be non-recombinants of the ancestral type
and $\rho^{2}$ is close to 1. The consequences of a few recombinants are
examined below using the following Julia function.


    function RSqr(nij) 
        N = sum(nij)
        Exy = nij[4]/N
        Ex = (nij[3] + nij[4])/N
        Ey = (nij[2] + nij[4])/N
        Vx = Ex * (1 - Ex)
        Vy = Ey * (1 - Ey)
        Cxy = Exy - Ex * Ey
        res = Cxy^2/(Vx * Vy)
        return(res)
    end




    RSqr (generic function with 1 method)



Here we only have non-recombinants:


    nij = [80,  # ancestral haplotype 00
           0,   # recombinant         01
           0,   # recombinant         10
           20]  # ancestral haplotype 11
    @printf "r^2 = %5.3f \n" RSqr(nij)

    r^2 = 1.000 


Now we introduce two recombinants:


    nij = [80,  # ancestral haplotype 00
           1,   # recombinant         01
           1,   # recombinant         10
           18]  # ancestral haplotype 11
    @printf "r^2 = %5.3f \n" RSqr(nij)

    r^2 = 0.874 


Here is another example, where one of the ancestral haplotypes has a
very low frequency:


    nij = [98,  # ancestral haplotype 00
           1,   # recombinant         01
           0,   # recombinant         10
           1]   # ancestral haplotype 11
    @printf "r^2 = %5.3f \n" RSqr(nij)

    r^2 = 0.495 


In a finite population, however, most loci are fixed. Then, LD is not
defined for these loci. When mutation introduces variability into such a
locus, LD is defined but will be low. This is demonstrated in the
following example:


    nij = [80,  # ancestral haplotype 00
           19,  # recombinant         01
           1,   # recombinant         10
           0]   # ancestral haplotype 11
    @printf "r^2 = %5.3f \n" RSqr(nij)

    r^2 = 0.002 


At mutation-drift equilibrium, most loci will be of this type, where
mutation has recently introduced variability. Thus, $E(\rho^{2})$ can be
much lower in a population that has reached mutation-drift equilibrium
than indicated by (3). To examine this further, the exact
distribution of is $\rho^{2}$ is recursively computed next.

#3. Distribution of $\rho^{2}$ in the Presence of Mutation


Computing the distribution of $\rho^{2}$ involves computing the joint
distribution for allele frequencies at two loci. Thus, we will review
first how to compute the distribution for allele frequency at a single
locus.

### Computing the distribution of allele frequency at a single locus

Consider a population of $2N$ gametes. Let $Y$ be the number of $A_{1}$
alleles at locus $A$. The value of $Y$ can take one of $2N+1$ values
ranging from 0 to $2N$. Suppose the distribution of allele frequency in
generation $t$ is given by the vector
${\mbox{\protect{\boldmath $p$}}}_{t}$ with $2N+1$ probabilities
corresponding to each of the $2N+1$ possible values of $Y$. To model
random mating, assume $2N$ gametes are sampled with replacement from the
gametes of generation $t$. Then, ignoring mutation, migration and
selection, the distribution of allele frequency in generation $t+1$ can
be calculated as

$${p}_{t+1}={Bp}_{t},$$

where is a $(2N+1)\times(2N+1)$ matrix with element $i,j$ containing the
probability that a random variable from a Binomial($2N,\frac{j}{2N}$)
distribution would be equal to $i$ for $i,j=0,1,2,\ldots,2N$. If this is
not obvious, section 3.10.2 of the notes given
[here](http://taurus.ansci.iastate.edu/rohan/notes/PopQuantGen.pdf) may
be useful.

To model mutation, assume that an $A_{1}$ allele mutates to an $A_{2}$
with probability $u$ and an $A_{2}$ mutates to an $A_{1}$ with
probability $v$. Now, to accommodate mutation in computing the
distribution of allele frequency, is modified such that column $j$
contains probabilities from the binomial distribution
$$\text{Binomial}[2N,\frac{j}{2N}(1-u)+(1-\frac{j}{2N})v].$$ Selection
can be similarly accommodated by modifying the binomial probabilities
for each $j$. See example
[here](http://taurus.ansci.iastate.edu/rohan/software/Qdist/index.html)

### Computing the joint distribution of allele frequencies at a two linked loci

Consider a locus $A$ with alleles $A_{1}$ and $A_{2}$ and a linked locus
$B$ with alleles $B_{1}$ and $B_{2}$. A population of $2N$ gametes is
now characterized by a vector ${\mbox{\protect{\boldmath $Y$}}}$ with
four elements containing the numbers of gametes with haplotypes:
$A_{1}B_{1}$, $A_{1}B_{2}$, $A_{2}B_{1}$, and $A_{2}B_{2}$. Note that
these four numbers must sum to $2N$. Let be a $k\times4$ matrix with
each row representing a possible value of
${Y}$, where the number $k$ of rows in is
equal to

$$k=\frac{(2N+3)!}{3!(2N)!}.$$

As before let
${p}_{t}$ denote the distribution of
haplotype frequencies in generation $t$. Then, ignoring recombination,
mutation, migration and selection, the distribution of haplotype
frequencies in the next generation are given by

$${p}_{t+1}={Mp}_{t},$$

where is a $k\times k$ matrix with element $i,j$ containing the
probability that a random variable from a
Multinomial$(2N,\frac{{x}_{j}'}{2N})$
distribution would be equal to ${x}_{i}'$,
for $i,j=1,2,\ldots,k$.

To model recombination, consider a population with frequency for
haplotype $A_{i}B_{j}$ given by $f_{ij}$. In gametes produced by this
population, the probability of a non-recombinant $A_{1}B_{1}$ is
$(1-r)\frac{f_{11}}{2N}$. A recombinant $A_{1}B_{1}$ gamete can be
produced in one of four ways. They and their associated probabilities
are:

1.  alleles $A_{1}$ and $B_{1}$ originate from two different
    $A_{1}B_{1}$ haplotypes with associated probability
    $r\frac{f_{11}}{2N}\times\frac{(f_{11}-1)}{2N-1}$;

2.  allele $A_{1}$ originates from an $A_{1}B_{1}$ haplotype and $B_{1}$
    originates from an $A_{2}B_{1}$ with associated probability
    $r\frac{f_{11}}{2N}\times\frac{f_{21}}{2N-1}$;

3.  allele $A_{1}$ originates from an $A_{1}B_{2}$ haplotype and $B_{1}$
    originates from an $A_{1}B_{1}$ with associated probability
    $r\frac{f_{12}}{2N}\times\frac{f_{11}}{2N-1}$; and

4.  allele $A_{1}$ originates from an $A_{1}B_{2}$ haplotype and $B_{1}$
    originates from an $A_{2}B_{1}$ with associated probability
    $r\frac{f_{12}}{2N}\times\frac{f_{21}}{2N-1}$.

Combining these probabilities gives
<br><br>

$$\begin{split}
\Pr(A_{1}B_{1})=
& (1-r)\frac{f_{11}}{2N}+\\
& r\frac{f_{11}}{2N}[\frac{(f_{11}-1)}{2N-1}+\frac{f_{21}}{2N-1}]+\\
& r\frac{f_{12}}{2N}[\frac{f_{11}}{2N-1}+\frac{f_{21}}{2N-1}].
\end{split}
$$

Similarly, probabilities of the remaining three types of gametes are:

$$\begin{split}
\Pr(A_{1}B_{2})=
& (1-r)\frac{f_{12}}{2N}+\\
& r\frac{f_{11}}{2N}[\frac{f_{12}}{2N-1}+\frac{f_{22}}{2N-1}]+\\
& r\frac{f_{12}}{2N}[\frac{(f_{12}-1)}{2N-1}+\frac{f_{22}}{2N-1}],
\end{split}
$$

$$\begin{split}
\Pr(A_{2}B_{1})=
& (1-r)\frac{f_{21}}{2N}+\\
& r\frac{f_{21}}{2N}[\frac{f_{11}}{2N-1}+\frac{(f_{21}-1)}{2N-1}]+\\
& r\frac{f_{22}}{2N}[\frac{f_{11}}{2N-1}+\frac{f_{21}}{2N-1}],
\end{split}
$$

and

$$\begin{split}\Pr(A_{2}B_{2})= & (1-r)\frac{f_{22}}{2N}+\\
 & r\frac{f_{21}}{2N}[\frac{f_{12}}{2N-1}+\frac{f_{22}}{2N-1}]+\\
 & r\frac{f_{22}}{2N}[\frac{f_{12}}{2N-1}+\frac{(f_{22}-1)}{2N-1}].
\end{split}
$$

Let ${\theta}_{j}$ be a vector with the
four probabilities from equations through computed for haplotype
frequencies from ${x}_{j}'$. Now, haplotype
probabilities following mutation can be modeled as

$${\beta}_{j}=T{\theta}_{j}$$

where
$${T}=\begin{bmatrix}(1-u)^{2} & (1-u)v & v(1-u) & v^{2}\\
(1-u)u & (1-u)(1-v) & vu & v(1-v)\\
u(1-u) & uv & (1-v)(1-u) & (1-v)v\\
u^{2} & u(1-v) & (1-v)u & (1-v)^{2}
\end{bmatrix}$$

<br><br>
To accommodate recombination and mutation in computing
the distribution of haplotype frequencies, is modified such that column
$j$ contains probabilities from the
Multinomial$(2N,{\beta}_{j}')$
distribution.

Starting with an allele frequency of 0.5 at each locus and gametic
equilibrium between the two loci, the expected value of $\rho^{2}$ was
computed for 2000 generations given a mutation rate of $u=v=1e^{-9}$, a
recombination rate of $r=0.002$ between the loci and an effective
population size of $N_{e}=5,10,25,$ or 50. The results are shown in the
figures [fig:r2N5] through [fig:r2N50]. In addition to $\rho^{2}$, the
figures also plot the frequencies of three groups of populations. In
populations that belong to group 1, $\rho^{2}<1$. In populations that
belong to group 2, two of the haplotypes are lost such that $\rho^{2}=1$
(for example, when haplotypes $A_{1}B_{2}$ and $A_{2}B_{1}$ are lost and
only haplotypes $A_{1}B_{1}$ and $A_{2}B_{2}$ are segregating,
$\rho=1$). In populations of group 3, one of the loci is fixed, and
therefore $\rho$ is not defined.

In all these cases, group 1 starts out having frequency close to 1.0.
Due to drift, however, the frequency in group 1 drops rapidly and
frequencies in groups 2 and 3 rise. The expected value of $\rho^{2}$
depends to a large extent on the relative magnitudes of the frequencies
of groups 1 and 2. Also, drift seems to reduce the frequency of group 1
faster than that of group 2. Therefore, $\rho^{2}$ rises rapidly and
then stays high for a period. Once frequencies in groups 1 and 2 drop
sufficiently low, changes in group frequencies due to mutation become
significant. Mutation in group 3, which has the highest frequency, adds
to group 1 faster than to group 2. Further, recombination in group 2
also contributes to group 1. The balance between these forces and drift
back into group 3 from groups 1 and 2 determines the equilibrium value
for $\rho^{2}$, which is reached by generation 2000 in all four cases.

![image](../images/R2N5R0p002Mu1em9.pdf)

[fig:r2N5]

![image](../images/R2N10R0p002Mu1em9.pdf)

[fig:r2N10]

![image](../images/R2N25R0p002Mu1em9.pdf)

[fig:r2N25]

![image](../images/R2N50R0p001Mu1em9.pdf)

[fig:r2N50]

The relationship between the equilibrium value of $\rho^{2}$, the
recombination rate, mutation rate, and effective population size can be
seen from figures [fig:N5] through [fig:N25], where for comparison
deterministic formulas by Sved @sved1971 and Hill @Hill75 are also
plotted. When mutation rate is zero, there is good agreement with Sved’s
formula! When mutation is present, agreement is better with Hill’s
formula.

![image](../images/N5Plot.pdf)

[fig:N5]

![image](../images/N10Plot.pdf)

[fig:N10]

![image](../images/N25Plot.pdf)

[fig:N25]



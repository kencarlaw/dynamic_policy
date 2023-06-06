Online Appendix for *"Enforcement Policy in a Dynamic Theory of Deterrence"*
============================================================================

This repository consists of material for replicating the core results from *"Enforcement Policy in a Dynamic Theory of Deterrence"* as well as a range of robustness tests. 

We have implemented all of our numerical simulations in two wholly separate codebases, Matlab and Python. All reported results have been independently replicated in each codebase. One member of our team wrote and maintained the Matlab codebase, another wrote and maintained the Python codebase.

# Replication of central results

## 1. The simulator:

The number of states in $\boldsymbol{Z}$ for $z=1$ and $N\leq 100$ is manageable, meaning it is feasible to produce the associated transition matrix $\boldsymbol{T}$ and use it to directly compute $\boldsymbol{D}$. But for $z>1$ and larger $N$ we must use numerical simulation methods. We now describe the simulation and our benchmarking efforts.

In each period of the simulation: Initiating a simulation requires values for: the quantity of enforcement resources ($R$); people's look back horizon ($z$); the number of potential violators ($N$); the mean ($\mu$) and standard deviation ($\sigma$) of the private gains from committing a violation; the probability a violator is apprehended if investigated ($\gamma$); the utility penalty paid by apprehended violators ($F$); and, the population's Bayesian priors regarding apprehension ($\alpha $, $\beta $). Unless noted, we set these parameters according to our baseline model.

* Every potential violator is assigned a $g$ drawn randomly from $\phi(g)$.
*Using the current z-history, a common subjective probability of apprehension $q$ is calculated.
* People choose to violate if and only if their $g^{t}_{i}\geq q^{t} \cdot F$.
* Violators are apprehended with probability $P(R^{t}, v^{t})$, where $v^{t}$ is the number of violators in the current period.
* The number of violations, apprehensions, and the sum of the gains from violating, the sum of the realized draws on $\phi(g)$ for those who violated, are recorded.
* The frequency distribution of violations and apprehensions is calculated for the current period and recorded.
* The z-history is updated and the simulation moves to the next period.

We calculate the frequency distribution of violations and apprehensions. The frequency distribution of violations is $\boldsymbol{FV}(nv_{0}^{t},nv_{1}^{t},\dotsc ,nv_{N}^{t})$, where $nv_{j}^{t}$ is the number of instances in the first $t$ periods of the simulation in which there were exactly $j$ violations. The frequency distribution of apprehensions is $\boldsymbol{FA}(na_{0}^{t},na_{1}^{t},\dotsc ,na_{N}^{t})$, where $na_{j}^{t}$ is the number of instances in the first $t$ periods of the simulation in which there were exactly $j$ apprehensions. 

The simulation proceeds in blocks of $50,000$ ticks. At the end of each block $n$ we use $\boldsymbol{FV}$ to calculate the relative frequency distribution of violations $\boldsymbol{RFV}^{n}=(rv_{0}^{n},\dotsc ,rv_{N}^{n})$, where $rv_{j}^{n}$ is the proportion of periods in which there were exactly $j$ violations in the first $50000\cdot n$ periods of the simulation. As $n$ approaches infinity  $\boldsymbol{RFV}^{n}$ approaches $\widehat{\boldsymbol{D}}$. The question is: How large must $n$ be to ensure that $\boldsymbol{RFV}^{n}$ is a \textit{good} approximation of $\widehat{\boldsymbol{D}}$? To answer this we need a convergence criterion. 

At the end of each simulation block, for all $n\geq 2$, we calculate the following test statistic:
\[TEST^{n}=\sum_{j=0}^{j=N}|rv_{j}^{n}-rv_{j}^{n-1}|\]
$TEST^{n}$ is a measure, decreasing in $n$, of the \textit{distance} between $\boldsymbol{RFV}^{n}$ and $\boldsymbol{RFV}^{n-1}$. $TEST^{n}$ does not decline monotonically, which complicates any test of convergence. After extensive experimentation we settled on the following criterion: we require that $TEST^{n} \leq 0.01$ for $5$ consecutive blocks of $n$. Denoting the minimum $n$ for which the convergence criterion is satisfied by $n^{\ast }$ yields our approximation of the stationary distribution, $\boldsymbol{RFV}^{n^{\ast }} \sim \widehat{\boldsymbol{D}}$.\footnote{\ An important confounding issue for finding $n^*$ is that the \textit{shape} of $\widehat{\boldsymbol{D}}$ changes dramatically as $R$ changes. For low values of $R$, $\widehat{\boldsymbol{D}}$ is unimodal and right skewed, for high values unimodal and left skewed, for intermediate values it becomes bimodal. This renders $n^*$ non-monotonic in R (see the last row of Table A.1).}

Our convergence simulator has some limitations over certain ranges of some key parameter values, in particular those that comprise the composite parameter $\frac{\sigma}{(\gamma \cdot F)}$ defined in Section \ref{section:CLIFFSET}. If $\frac{\sigma}{(\gamma \cdot F)}$ is small the $n$ required to achieve a good approximation of the stationary distribution is actually much larger than the $n^{\ast}$ generated from the convergence criteria. This is because of the persistence of the q-attractors. For these parameterizations the simulator is highly sensitive to initial conditions, staying for a very long time (sometimes the entire simulation) in only one of the q-attractors and rarely (if ever for some convergence simulations around the cliff) transitions to the other q-attractor, even though the stationary distribution is bi-model in these parameter values. This yields a poor estimate of the stationary distribution even though the convergence simulator may have taken a very long time to converge. More accurate estimations of the stationary distribution for such parameter values would take orders of magnitude longer, longer than is feasible to make simulation useful. It is important to note that this is not a limitation of the model as such. But rather a limitation of our technology of simulation to approximate the stationary distribution for the regular finite Markov chain in these parameter values. 

### 1.a. Benchmarking the simulator

We employ Monte Carlo simulations and the relationship between $\widehat{\boldsymbol{D}}$ and $\boldsymbol{RFV}^{n^{\ast}}$ to: (i) cross validate the output from our Python and Matlab codebases; and, (ii) benchmark the accuracy of the simulator against a known transition matrix.\footnote{\ We are, of course, implicitly testing the validity of our convergence criterion.} The computational intensity of producing both $\widehat{\boldsymbol{D}}$ and $\boldsymbol{RFV}^{n^{\ast}}$ increases geometrically with the size of the state space. It is only practical, therefore, to implement both approaches for a relatively small number of agents ($N=50$) and a single look back period ($z=1$). For these parameters, and for each element of $R=\{5,21,45\}$, we generate $1000$ instances of $\boldsymbol{RFV}^{n^{\ast }}$. We then construct the following distance metric: 
\[ATEST=\widehat{\boldsymbol{D}}-\boldsymbol{RFV}^{n^{\ast }} = \sum_{j=0}^{j=N}|\widehat{d}_{j}-rv_{j}^{n^{\ast }}|\]

In Table 4.1 we report the mean, standard deviation, and maximal value of $ATEST$ along with the mean $n^*$. It is apparent that the two platforms deliver results that are essentially the same, and that $\boldsymbol{RFV}^{n^{\ast }}$ very closely approximates $\widehat{\boldsymbol{D}}$.


| Syntax      |
| Syntax      | Description | Test Text     |
| :---        |    :----:   |          ---: |
| Header      | Title       | Here's this   |
| Paragraph   | Text        | And more      |

\begin{table}[h]
\centering
\caption{Benchmarking the simulator}
\begin{tabular}{lccc;{1pt/1pt}ccc} 
\toprule
 & \multicolumn{3}{c} \textbf{Matlab Codebase} & \multicolumn{3}{c} \textbf{Python Codebase} \\ 
\hline
 & \multicolumn{6}{c}{Value of R} \\
 & \textbf{5} & \textbf{21} & \textbf{45} & \textbf{5} & \textbf{21} & \textbf{45} \\
 \hline
Mean of ATEST & 0.0017 & 0.0098 & 0.0078 & 0.0018 & 0.0115 & 0.0087 \\
Std Dev of ATEST & 0.0009 & 0.0017 & 0.0011 & 0.0009 & 0.0019 & 0.0011 \\
Max value ATEST & 0.0061 & 0.0193 & 0.0116 & 0.0057 & 0.0198 & 0.0128 \\
Mean $n^*$ & 6.0010 & 7.8430 & 7.0870 & 6.0010 & 7.8120 & 6.0010 \\
\bottomrule
\end{tabular}
\end{table}

# Robustness tests

1. **Search algorithm for active policies*

We use a two stage procedure to identify the optimal crackdown and refined crackdown policies.  In the first stage, we use a directed search routine with a loose convergence criterion to identify the neighborhood in which the optimal policy lies, and in the second stage we use Monte Carlo methods to generate tight estimates of the cost of policies in the neighborhood. 

The algorithm at the core of the directed search routine is this. For a given crackdown policy at stage $s,$ $CD^{s}=(\mathbb{BE}^{s},\mathbb{R_{GB}}^{s},\mathbb{R_{BB}}^{s})$, a modified policy $CD^{s+1}=(\mathbb{BE}^{s+1},\mathbb{R_{GB}}^{s+1},\mathbb{R_{BB}}^{s+1})$, is generated using a three step procedure: in the first step, $\mathbb{BE}^{s+1}$ is chosen from the set $\{\mathbb{BE}^{s}-5,\mathbb{BE}^{s}-4,\dotsc ,\mathbb{BE}^{s}+5\}$ to minimize the cost of ASB for policy $(\mathbb{BE}^{s+1},\mathbb{R_{GB}}^{s},\mathbb{R_{BB}}^{s})$; in the second step, $\mathbb{R_{GB}}^{s+1}$ is chosen from the set $\{\mathbb{R_{GB}}^{s}-5,\mathbb{R_{GB}}^{s}-4,\dotsc ,\mathbb{R_{GB}}^{s}+5\}$ to minimize the cost of ASB for policy $(\mathbb{BE}^{s+1},\mathbb{R_{GB}}^{s+1},\mathbb{R_{BB}}^{s})$; in the third step, $\mathbb{R_{BB}}^{s+1}$ is chosen from the set $\{\mathbb{R_{BB}}^{s}-5,\mathbb{R_{BB}}^{s}-4,\dotsc ,\mathbb{R_{BB}}^{s}+5\}$ to minimize the cost of ASB for policy $(\mathbb{BE}^{s+1},\mathbb{R_{GB}}^{s+1},\mathbb{R_{BB}}^{s+1})$. For each of the 33 policies considered, we ran a convergence simulation and used it to estimate cost. The search is terminated when the criteria $(E(C|CD^{s})-E(C|CD^{s+1}))/E(C|CD^{s})<.02$ is met. Each of the 1000 simulations was seeded with an initial policy randomly chosen from the set $\{(\mathbb{BE},\mathbb{R_{GB}},\mathbb{R_{BB}})| 0\leq \mathbb{BE}\leq 100, 0\leq \mathbb{R_{GB}}\leq 100, \mathbb{R_{GB}}\leq \mathbb{R_{BB}}\leq 100\}$ and for each we used the directed search algorithm to a generate a terminal policy. From initial policy to terminal policy, the average number of policies evaluated was 165. To identify the neighborhood of the optimal crackdown policy, we ranked the 1000 terminal policies from lowest to highest cost. The first 78 policies in this ranking, and 98 of the first 100, were in the set $\{(\mathbb{BE},\mathbb{R_{GB}},\mathbb{R_{BB}})| 30\leq \mathbb{BE}\leq 39, 29\leq \mathbb{R_{GB}}, 3\leq \mathbb{R_{BB}}\leq 65\}$. This is neighborhood of the optimal policy. 

The search for the optimal refined crackdown policies was analogous to that for the crackdown policy, with the differences being the need for two additional steps in the procedure to search over the two additional parameters, $\mathsf{BE2}$ and $\mathsf{R_{BB1}}$, of the refined crackdown policy. The algorithm for the refined crackdown policy search is, given a refined crackdown policy at stage $s,$ $CD^{s}=(\mathsf{BE}^{s},\mathsf{BE2}^{s},\mathsf{R_{GB}}^{s},\mathsf{R_{BB1}}^{s},\mathsf{R_{BB2}}^{s})$, a modified policy $CD^{s+1}=(\mathsf{BE}^{s+1},\mathsf{BE2}^{s+1},\mathsf{R_{GB}}^{s+1},\mathsf{R_{BB1}}^{s+1},\mathsf{R_{BB1}}^{s+1})$, is generated using a five step procedure analogous to that of the three step procedure of the crackdown policy search. For each of the 55 policies considered, we again ran a convergence simulation and used it to estimate the cost. The termination criteria for the refined crackdown policy search was the same as that used for the crackdown policy search. The 1000 simulations were randomly seeded with initial refined crackdown policies chosen from the set $\{(\mathsf{BE},\mathsf{BE2},\mathsf{R_{GB}},\mathsf{R_{BB1}},\mathsf{R_{BB2}})|0\leq \mathsf{BE}\leq 100, \mathsf{BE}\leq \mathsf{BE2}\leq 100, 0\leq \mathsf{R_{GB}}\leq 100, \mathsf{R_{GB}}\leq \mathsf{R_{BB1}}\leq \mathsf{R_{BB2}}, \mathsf{R_{BB1}}\leq \mathsf{R_{BB2}}\leq 100\}$. From initial policy to terminal policy, the average number of policies evaluated was 440 for each of the 1000 simulations. The neighborhood of the optimal refined crackdown policies was determined in the same way as for the crackdown policy.  

The estimate of cost that comes out of one convergence simulation is, obviously, a random variable. In the neighborhood of the optimal policy, the cost function is so flat that in order to reliably identify the optimum policy, many independent cost estimates are needed. Accordingly, for every policy in the neighborhood of the optimum we ran 150 convergence simulations and calculated the cost of ASB in the steady state for each of them. Our cost estimate is the mean value of these 150 cost estimates. Typically the standard deviation of cost for the 150 simulations is close to 0.20, and hence the standard error of the estimate is close to $0.016=.2/\sqrt{150}$.

2. **Finding optimal refined crackdown**

In the first stage we identify the neighborhood of the optimal refinement by using a straight forward extension of the directed search routine used in Section \ref{section:optimalcrackdown}. For each of the 1000 randomly chosen policies, the directed search routine identifies a candidate policy. In the second stage, we rank the 1000 candidate policies in ascending order of their cost estimates. We then generate tight costs estimates for the first 25 policies in this ranking, running 150 independent convergence simulations for each policy and using the mean cost as our cost estimate for the policy.

You may have noticed that the first two elements in our notation for a refined crackdown policy, $\mathsf{BE}$ and $\mathsf{R_{GB}}$, are the same as those used for the first two elements of a crackdown policy. This reflects a conscious choice to frame a refined crackdown policy as a refinement of a crackdown policy with the same $\mathbb{BE}$ and $\mathbb{R_{GB}}$ in which the bad bin, $BB$, is partitioned into two bins, $BB_1$ and $BB_2$, with different quantities of the enforcement resource. We chose to frame them in this way to highlight the fact the 25 least costly refined crackdown policies we identify in stage 1 of our search procedure are in fact refinements of one of the least costly crackdown policies we identified in the previous section, and presented in Table \ref{table: activecostestimates}. Key results for the 25 refined crackdown policies with the lowest tight cost estimates are reported in Table \ref{table: refinedcostestimates}. Notice that the skeletons of Tables \ref{table: activecostestimates} and \ref{table: refinedcostestimates} are identical. 

\begin{table}[h]
\centering
\caption{The Value of Refinement.}
\label{table: refinedcostestimates}
\begin{threeparttable}
\scriptsize
\begin{tabular}{c c c c c c c c c c c c c}
\toprule
\multicolumn{1}{l}{}              &\multicolumn{1}{l}{}              & \multicolumn{10}{c}{$\mathbf{\mathsf{BE}}$} &\multicolumn{1}{l}{}\\ 
\multicolumn{1}{l}{}              &\multicolumn{1}{l}{}              & \textbf{30}     & \textbf{31}     & \textbf{32}      & \textbf{33}      & \textbf{34}      & \textbf{35}      & \textbf{36}      & \textbf{37}      & \textbf{38}      & \textbf{39} &\multicolumn{1}{l}{} \\ \cline{1-13} 
\multicolumn{1}{l}{}              &\multicolumn{1}{c}{} & & & & & & & & & &\multicolumn{1}{l}{} \\
\multicolumn{1}{l}{}              &\multicolumn{1}{c}{\textbf{29}}  & \begin{tabular}[c]{@{}c@{}}1.07\\ \tiny{9,18}\end{tabular}  & \begin{tabular}[c]{@{}c@{}}0.73\\ \tiny{11}\end{tabular} & \begin{tabular}[c]{@{}c@{}}0.62\\ \tiny{6}\end{tabular} & \begin{tabular}[c]{@{}c@{}}XXX\\  \end{tabular} & \begin{tabular}[c]{@{}c@{}}XXX\\  \end{tabular} & \begin{tabular}[c]{@{}c@{}}XXX\\  \end{tabular} & \begin{tabular}[c]{@{}c@{}}XXX\\  \end{tabular} & \begin{tabular}[c]{@{}c@{}}XXX\\  \end{tabular} & \begin{tabular}[c]{@{}c@{}}XXX\\  \end{tabular} & \begin{tabular}[c]{@{}c@{}}XXX\\  \end{tabular} &\multicolumn{1}{l}{} \\
\multicolumn{1}{l}{}              &\multicolumn{1}{c}{\textbf{}}    & & & & & & & & & &\multicolumn{1}{l}{} \\
\begin{tabular}[c]{@{}c@{}} \\ $\mathbf{\mathsf{R_{GB}}}$ \end{tabular}  &\multicolumn{1}{c}{\textbf{30}}  & \begin{tabular}[c]{@{}c@{}}XXX\\  \end{tabular}  & \begin{tabular}[c]{@{}c@{}}XXX\\  \end{tabular}  & \begin{tabular}[c]{@{}c@{}}0.68\\ \tiny{2,7,15}\end{tabular}  & \begin{tabular}[c]{@{}c@{}}0.49\\ \tiny{3,4,10}\end{tabular} & \begin{tabular}[c]{@{}c@{}}0.42\\ \tiny{1,5,8,13,14,21,22}\end{tabular} & \begin{tabular}[c]{@{}c@{}}0.17\\ \tiny{20}\end{tabular} & \begin{tabular}[c]{@{}c@{}}XXX\\  \end{tabular} & \begin{tabular}[c]{@{}c@{}}XXX\\  \end{tabular} & \begin{tabular}[c]{@{}c@{}}XXX\\  \end{tabular} & \begin{tabular}[c]{@{}c@{}}XXX\\  \end{tabular} &\multicolumn{1}{l}{} \\
\multicolumn{1}{l}{}              &\multicolumn{1}{c}{\textbf{}}    & & & & & & & & & &\multicolumn{1}{l}{} \\
\multicolumn{1}{l}{}              &\multicolumn{1}{c}{\textbf{31}}  & \begin{tabular}[c]{@{}c@{}}XXX\\  \end{tabular}  & \begin{tabular}[c]{@{}c@{}}XXX\\  \end{tabular} & \begin{tabular}[c]{@{}c@{}}XXX\\  \end{tabular}  & \begin{tabular}[c]{@{}c@{}}0.55\\ \tiny{17}\end{tabular} & \begin{tabular}[c]{@{}c@{}}0.42\\ \tiny{12,24}\end{tabular} & \begin{tabular}[c]{@{}c@{}}0.20\\ \tiny{19,25}\end{tabular} & \begin{tabular}[c]{@{}c@{}}0.24\\ \tiny{16,23}\end{tabular} & \begin{tabular}[c]{@{}c@{}}XXX\\  \end{tabular} & \begin{tabular}[c]{@{}c@{}}XXX\\  \end{tabular} & \begin{tabular}[c]{@{}c@{}}XXX\\  \end{tabular} &\multicolumn{1}{l}{}\\
\multicolumn{1}{l}{}              &\multicolumn{1}{c}{\textbf{}}    & & & & & & & & & &\multicolumn{1}{l}{} \\
\multicolumn{1}{l}{}              &\multicolumn{1}{c}{\textbf{32}}  & \begin{tabular}[c]{@{}c@{}}XXX\\  \end{tabular}  & \begin{tabular}[c]{@{}c@{}}XXX\\  \end{tabular} & \begin{tabular}[c]{@{}c@{}}XXX\\  \end{tabular}  & \begin{tabular}[c]{@{}c@{}}XXX\\  \end{tabular} & \begin{tabular}[c]{@{}c@{}}XXX\\  \end{tabular} & \begin{tabular}[c]{@{}c@{}}XXX\\  \end{tabular} & \begin{tabular}[c]{@{}c@{}}XXX\\  \end{tabular} & \begin{tabular}[c]{@{}c@{}}XXX\\  \end{tabular} & \begin{tabular}[c]{@{}c@{}}XXX\\  \end{tabular} & \begin{tabular}[c]{@{}c@{}}XXX\\  \end{tabular} &\multicolumn{1}{l}{} \\ \bottomrule \end{tabular}
\end{threeparttable}
\end{table}

All 25 of the best refined crackdown policies appear in just 11 of the 40 cells in Table \ref{table: refinedcostestimates}, and 13 of them are in cells (32, 30), (33, 30) and (34, 30). In the top line of each cell we have listed the percentage difference in cost between the best refined crackdown policy and crackdown policy for that cell (from Table \ref{table: activecostestimates}). The second line of each cell identifies the refined crackdown policy by its rank---with rank 1 being the lowest cost estimate. Policies 1, 5, 8, 13, 14, 21 and 22 are listed in cell (34, 30). All are refinements of the crackdown policy (34, 30, 56) from Table \ref{table: activecostestimates}, and all produce a lower cost. This is not surprising -- with more instruments to control ASB we expect to see a lower cost of ASB. What may be surprising is that the reduction in cost is relatively small --- that is what the entry in the first line of each cell addresses. The cost estimate for the refined crackdown policy with rank 1 is only 0.41\% lower than the costs estimate for optimal crackdown policy (34, 30, 56). 

In Table \ref{table: stderrorscostestimates} we report the policies that are ranked 1, 2 and 3, our costs estimates for them, and the standard errors of the cost estimates. When we look at the policies themselves we see that the refinements make intuitive sense. The lowest cost refined crackdown policy is $(\mathsf{BE},\mathsf{R_{GB}},\mathsf{SBE},\mathsf{R_{BB1}},\mathsf{R_{BB2}})=(34,30,44,47,64)$ and it is a refinement of crackdown policy $(\mathbb{BE},\mathbb{R_{GB}},\mathbb{R_{BB}})=(34,30,56)$. When crackdown policy $(34,30,56)$ is used, 56 units of enforcement resource are deployed in the bad bin, $\{v^{t-1}|34 < v^{t-1}\leq 100\}$. When refined crackdown policy $(34,30:44,47,64)$ is used the bad bin is partitioned into a two two bins, the first is $\{v^{t-1}|34 < v^{t-1}\leq 44\}$ and the second is $\{v^{t-1}|44 < v^{t-1}\leq 100\}$, and quantity of enforcement resources deployed in the first is 47, somewhat less 56, and quantity deployed in the second is 64, somewhat greater than 56. In the bad bin, the proximate objective of policy is to drive violations down and back into $GB$, and it makes intuitive sense that the quantity of resources needed to do that efficiently is smaller when violations are close to the lower bound of $BB$ than when violations are close to the upper bound.

\begin{figure}[h]
    \centering
    \includegraphics[width=0.80\textwidth]{Penultimate Draft/Images/Fig 17.png}
    \caption{Cross sections of costs in the neighborhood of $RCD\#1$}
    \label{fig:RCD}
\end{figure}

In Section \ref{section:optimalcrackdown}, to identify the optimal crackdown policy, we generated a
tight cost estimate for every policy in the neighborhood that we identified using the stage one directed search routine. That approach is not feasible for refined crackdown policies. We can use the directed search routine to identify the neighborhood of the optimum, but there are so many policies in that neighborhood that it is not feasible to generate a tight cost estimate for each of them. In addition, in three of the five relevant dimensions, $\mathsf{SBE}$, $\mathsf{R_{BB_{1}}}$ and $\mathsf{R_{BB_2}}$, the cost function is so flat that it would take thousands of convergence simulations to get cost estimates that would allow us to identify the optimal policy. 

Nevertheless, the results reported in Figure \ref{fig:RCD} suggest that the cost of ASB for the optimal policy is very close to the cost estimates in Table \ref{table: stderrorscostestimates}. The cost estimates reported in the figure are based on 50 independent convergence simulations for each policy examined (the standard errors of these estimates are roughly 0.03). In each panel, cost estimates are reported for 19 policies centered on policy \#1: four of the five parameters of the policy are fixed at their values in policy \#1, and the fifth varies up and down from its value in policy \#1 by 9 units. The asterisk indicates policy \#1, and the dot the policy with the lowest cost estimate. There is considerable curvature in the cost function when $\mathsf{R_{GB}}$ and $\mathsf{BE}$ vary, which means that the optimal values of these parameters are much easier to pin down, and that is exactly what we saw in Table \ref{table: stderrorscostestimates}. But over the range of values examined for parameters $\mathsf{SBE}$, $\mathsf{R_{BB_{1}}}$ and $\mathsf{R_{BB_{2}}}$, the cost function is very flat and it is much more difficult to pin down the optimal values of theses parameters. This flatness reflects the fact that there are many ways to extinguish $BB$, and, it does matter much in terms of costs precisely how that is achieved because so little time is spent there. The inserts in each panel present the same information but the scale on the vertical axis is different. The more fine grained scale allows us to see how the cost estimates for the refined crackdown policies compare to our best estimate of the minimum cost achievable with the simpler crackdown policy---the dotted reference line in each insert represents that cost.


\section{Generality of deterrence policy analysis} \label{section:genODP}
In this concluding section we start by considering the generality and determinants of three central results for the baseline parameterization: (i) The reserve capacity rate with the optimal crackdown policy is higher than 60\%, meaning that in a typical period only 40\% of available enforcement resources are used to investigate violations; (ii) replacing the optimal passive policy with the optimal crackdown policy reduces the expected cost of ASB by 12\%; and, (iii) replacing the optimal crackdown policy with the optimal refined crackdown policy reduces the expected cost of ASB by less than 1\%. Naturally, these specific results, the reserve capacity rate and the measured reduction in costs, will differ from one parameterization to another. We show that they are driven primarily by $RAT$ and $DROP$.  

$RAT$ is effectively the ratio of the external cost imposed on society by a violation to the cost of investigating a violation, so it is intuitive that it is an important determinant of these three numbers. As we have seen, managing ASB means managing positive feedback. Since $DROP$ is a measure of the potential for disruptive feedback it is intuitive that it too is important to understanding these results.

\begin{figure}[h]
    \centering
    \includegraphics[width=0.95\textwidth]{Penultimate Draft/Images/Fig 16d.png}
    \caption{Reserve capacity and cost savings as a function of $DROP$.}
    \label{fig:dropscatter}
\end{figure}

In Figures \ref{fig:dropscatter} and \ref{fig:ratscatter} we report the reserve capacity rate and the two cost reduction percentages for  simulations in which $N$ and $z$ are fixed at their baseline values.\footnote{\ The estimates of optimal crackdown and refined crackdown policies reported in this section are generated using simpler versions of the directed search algorithms introduced in Section \ref{section:activepolicies}. The number of randomly seeded searches was reduced to 50 in each case implying that the errors in the estimation of the optimal policy will be slightly higher. But as we observed in Section \ref{section:activepolicies} the space in which the optimum occurs is very flat and so there is very little difference between estimates in the neighborhood of the optimum. Reducing the number of searches increased computational efficiency at nearly zero cost to estimation precision.} In these figures we are treating $DROP$ as something like a composite parameter. Given $N=100$ and $z=2$, as we saw in Section \ref{section:robustness}, $DROP$ is determined by $\frac{\mu}{\gamma \cdot F}$ and $\frac{\sigma}{\gamma \cdot F}$. In selecting the parameterizations used in constructing these figures, we used two criteria. First, that the distribution of $DROP$ be approximately uniform over the interval (0.275, 0.725). Second, subject to this constraint, that we include parameterizations distributed across the entire $(\frac{\mu}{\gamma \cdot F},\frac{\sigma}{\gamma \cdot F})$ parameter space. For each parameterization we calculated reserve capacity and the two cost reduction numbers. Then, in Figures \ref{fig:dropscatter} and \ref{fig:ratscatter}, we plot the values for our three central results for each parameterization with their values of $DROP$ and $RAT$. 

The first thing to notice is the magnitude of the correlation coefficients reported in the four panels of Figure \ref{fig:dropscatter}. The fact that they are relatively high is telling us that $DROP$ does function as a composite of the parameters that determine it, and the fact they are not $1$ is telling us that the effects of variation in these parameters are not completely captured by $DROP$. 

In Panels 1 and 3 of Figure \ref{fig:dropscatter}, we report the reserve capacity rate for the optimal crackdown policy for a number of DROP values in the interval (0.25, 0.75). In Panel 1, $RAT$ is fixed at 2 and in 3 it is fixed at 6. In both panels, reserve capacity is increasing in $DROP$, from a low of 20\% up to high of 70\%. In Panels 2 and 4 the circles represent the reduction in cost when the optimal passive policy is replaced with the optimal crackdown policy. In Panel 2 $RAT$ is 2, and in 4 it is 6. The cost reduction numbers are  positively associated with $DROP$ and  negatively associated with $RAT$. The crosses in those panels represent the reduction in cost when the optimal crackdown policy is replaced with the optimal refined crackdown policy. These cost reduction numbers are negligible and uncorrelated with $DROP$.

\begin{figure}[h]
    \centering
    \includegraphics[width=0.95\textwidth]{Penultimate Draft/Images/Fig 18d.png}
    \caption{Reserve capacity and cost savings as a function of $RAT$.}
    \label{fig:ratscatter}
\end{figure}

In Panel 1 of Figure \ref{fig:ratscatter}, we report the reserve capacity rate with the optimal crackdown policy as a function of $RAT$ for three $DROP$ values. For all three values, there is massive discontinuity in the relationship between reserve capacity and $RAT$. (Given the coarse grid of $RAT$ values, we have not identified any of the three discontinuities precisely.) When $RAT$ is very small, we are in the spitting on the sidewalk case and the best enforcement policy is to do nothing. As $RAT$ passes through the point of discontinuity the reserve capacity rate jumps up, and beyond the discontinuity it increases slowly with $RAT$. For the two larger $DROP$ values to the right of the discontinuity the rate exceeds 60\% and for the lowest value it is much lower. In Panel 2, we report the reduction in cost when the optimal passive policy is replaced with the optimal crackdown policy as a function of $RAT$ for the same three values of $DROP$. Again there is massive discontinuity in the relationship. When $RAT$ is small do nothing is the optimal enforcement policy regardless of whether we are using a passive or crackdown policy, so the cost reduction is 0. Beyond the discontinuity the cost reduction is inversely related to $RAT$.

The results reported in Figures \ref{fig:dropscatter} and \ref{fig:ratscatter} are based on parameterizations in which $N=100$ and $z=2$. Because simulations with larger values of these parameters are very time intensive, we have run just a few of them. Based on those simulations it seems clear that the patterns seen in these figures and the magnitudes of the numbers are robust to increases in the values of the parameters.\footnote{ \ Variations of $N=160$ and $z=3$ were run on two values each of $DROP\approx\{0.7, 0.4\}$ and $RAT=\{2, 10\}$. The pattern of cost savings for crackdown and refined crackdown was preserved. Costs savings for crackdown over passive increased slightly for both the larger $N$ and $z$ values, with the increases ranging from 0.1\% to 1.8\% over the baseline values with the maximum cost savings being 15.2\%. The cost savings for refined crackdown was again negligible with the total savings ranging from 0\% to 0.4\%. The reserve capacity values were almost identical for each for each of the four cases. For example, for $DROP\approx0.7$ and $RAT=2$, reserve capacity = (0.4260, 0.4240, 0.4302) for the $(N, z)$ tuples of (100, 2) (160, 2) and (100, 3) respectively.}

Our results support the following generalizations.

\begin{proposition} \label{DROP}
When the potential for positive feedback as measured by DROP is large: (i) to manage antisocial behavior in a cost effective way substantial reserve enforcement capacity is required; (ii) the cost of antisocial behavior with the optimal crackdown policy is substantially lower than it is with the optimal passive policy; and, (iii) the cost of antisocial behavior with the optimal refined crackdown policy is negligibly less than it is with the optimal crackdown policy.
\end{proposition}

\section{General features of the stochastic dynamic process}  \label{section:robustness}
As we have seen, our dynamic model produces a suite of features including the cliff of positive feedback. Proposition \ref{Genfeatures} compactly describes these features. 

\begin{proposition}\label{Genfeatures}
\begin{enumerate}
\onehalfspacing
\item[]
\textit{\item[(i)] As $R$ increases from 0 one can distinguish three distinct regimes: initially there is an unruly regime where $E(v|R)$ is relatively high and slowly decreasing in $R$, followed by a transitional regime in which $E(v|R)$ drops precipitously as $R$ increases, and finally a compliant regime where $E(v|R)$ is relatively low and virtually constant.}
\textit{\item[(ii)] There is a single q-attractor in the unruly and compliant regimes. In the unruly regime the number of violations at the focal point of the q-attractor $BB$ is relatively high and decreases slowly with $R$; in the compliant regime the focal point of the q-attractor $GB$ is relatively low and unresponsive to $R$.}
\textit{\item[(iii)] There are two q-attractors in the transitional regime: $GB$ where the focal point violations are low and $BB$ where the focal point violations are high. As $R$ transits between the unruly and compliant regimes, the persistence of $BB$ diminishes from approximately 1 to roughly 0 while the persistence of $GB$ increases from approximately 0 to roughly 1.}
\textit{\item[(iv)] Items (ii) and (iii) produce a prominent cliff-like structure in the relationship between the quantity of enforcement resources $R$ and the expected number of violations in the stationary distribution $E(v|R)$.}
\textit{\item[(v)] In all three regimes there is significant positive autocorrelation in the time series of violations for lags of one, two and three periods and for much longer lags in the transitional regime.}
\end{enumerate}
\end{proposition}

Our focus for the remainder of the paper is on enforcement policies that minimize the total costs of ASB. We discover that managing the dynamics of positive feedback is the primary policy problem. Showing that the cliff is a general feature of our model establishes that the dynamic features which produce it are themselves general, and, therefore, that our policy insights generalize to a set of parameterizations beyond the baseline.

In Section \ref{section:CLIFFSET} we develop a measure of the prominence of the cliff called $DROP$. Prominent cliffs produce large values of $DROP$, inconspicuous cliffs produce low values of $DROP$. We then define $CLIFFSET$ as the set of parameterizations such that $DROP$ exceeds a threshold value. In subsequent sections we report evidence that DROP itself is an economically interesting statistic.

In Section \ref{section:robustnessassumptions} we present suggestive evidence that the cliff of positive feedback and the dynamics that produce it are also robust with respect to changes in the assumptions which instantiate the general Markov process in our specific computational model.

\subsection {CLIFFSET}\label{section:CLIFFSET}
The geometric signature of the cliff is a small range of contiguous values of $R$ such that as $R$ is increased through that range there is a dramatic drop in $E(v|R)$. We set the range of $R$ to $RANGE=INT(0.05\cdot N)$. For the baseline parameterization $RANGE=5$. 

For a given a parameterization we first estimate the relationship $E(v|R)$ for all $R\in \{0,...,N\}$. For every integer $X\in \{0,...,N-RANGE)$ we then calculate the reduction in $E(v|R)$ when $R$ is increased from $X$ to $X+RANGE$. $DROP$ is the maximum value of this difference, normalized by expressing it as a fraction of its maximum possible value, $N$. 
\[DROP=\frac{\underset{\{X\}}{\max }[E(v|R=X)-E(v|R=X+RANGE)]}{N}\]

Intuitively, $DROP$ measures the height of the cliff relative to the size of the population of potential violators. It is bounded below by 0 and above by 1. For example, $DROP=0.5$ indicates that in the vicinity of the cliff an increase of $INT(0.05\cdot N)$ in $R$ reduces the expected number of violations in the stationary distribution by $0.5\cdot N$, or that the slope of the cliff is approximately $-10$. In defining $CLIFFSET$, we set the threshold value of $DROP$ to $0.35$, producing a slope of $-7$. The value of $DROP$ for the baseline parameterization is roughly $0.7$, producing a cliff with a slope of $-14$.

Figure \ref{fig:thresholdcliff} plots three different cliffs: one for the baseline parameterization, one in which $DROP = 0.35$, and in which $DROP = 0.1$. The baseline produces the most prominent cliff and the largest value of DROP.   

\begin{figure}[h]
    \centering
    \includegraphics[width=0.65\textwidth]{Penultimate Draft/Images/Fig 8.png}
    \caption{Three cliffs of positive feedback and their measured $DROP$.}
    \label{fig:thresholdcliff}
\end{figure}

$DROP$ is a function of six parameters, $\mu ,\sigma ,\gamma ,F,N$ and $z$, so $CLIFFSET$  is an object in a six dimensional space. 
\[CLIFFSET=\{\mu ,\sigma ,\gamma ,F,N,z|DROP\geq 0.35\}\]

While a complete description is not feasible, we can describe some of the most important aspects of $CLIFFSET$. Since the underlying behavior of the model depends only on the value of $\gamma \cdot F$, and not independently on $\gamma $ or $F$, we are able to reduce $CLIFFSET$'s dimensionality by treating $\gamma \cdot F$ as a composite parameter. Think of $\gamma \cdot F$ as the size of the deterrence stick. We start by describing a subset of $CLIFFSET$, 
\[CLIFFSET^\prime= \{\mu ,\sigma \ | DROP\geq 0.35, \ \gamma \cdot F=0.80, \ N=100, \ z=2\}\] 

In Figure \ref{fig:dropsurface}, we have constructed a surface plot of $DROP$ for $\mu \in\{0,0.05,0.1,...,1.00\}$ and $\sigma \in\{0.15,0.20,...,0.55\}$ (so the set of $(\mu, \sigma)$ pairs has $189$ elements).\footnote{\ We have not included estimates of DROP for values of $\sigma $ less than $0.10$ because our convergence simulation routine is not reliable for small values of $\sigma$. See Appendix \ref{appendix:simulator}.} We have included a shaded reference plane representing our threshold value of $DROP$. The outer boundary of $CLIFFSET^\prime$ is defined by the locus of points where the surface plot and the reference plane intersect. We have projected the iso-DROP locus for $DROP=0.35$ onto the baseplane, $CLIFFSET^\prime$ is the set of $(\mu ,\sigma)$ pairs on and inside the iso-DROP locus.

\begin{figure}[h]
    \centering
    \includegraphics[width=0.95\textwidth]{Penultimate Draft/Images/Fig 9.png}
    \caption{$DROP$ for a range of $\mu$ and $\sigma$.}
    \label{fig:dropsurface}
\end{figure}

There is another interpretation of $CLIFFSET^\prime$. Setting $s>0$ and applying this scalar to three key parameters of $CLIFFSET'$, $(s \cdot \mu ,s \cdot \sigma, s \cdot \gamma \cdot F)$, the ratio of the new set of values for $DROP^\prime$ and the original $DROP$ values are, to a very close approximation, equal to $1$. In other words, $DROP$ is homogeneous of degree 0 ($HD0$) in $\mu ,\sigma$ and $\gamma \cdot F$. This follows from the fact that the probability $\pi$ that any person chooses violate in any period is $HD0$ in these same parameters. In every period, every person has a well defined subjective probability of apprehension, $q$. When their opportunity $g$ is revealed they choose violate if and only if $g>q\cdot \gamma \cdot F$. Prior to the realization of $g$, the probability $\pi$ that they choose violate is $1-\Phi (q \cdot \gamma \cdot F)$. Since $\phi (g)$ is normal, $\Phi (q\cdot \gamma \cdot F)$ and hence $\pi$, are themselves $HD0$ in $\mu ,\sigma$ and $\gamma \cdot F$. 

\begin{proposition}\label{DROPHD0}
Given any parameterization of the model $(\widehat{\mu },\widehat{\sigma },\widehat{\gamma }\cdot \widehat{F},\widehat{z},\widehat{N},\widehat{R})$ and the constructed parameterization $(s\cdot \widehat{\mu },s\cdot \widehat{\sigma },s\cdot \widehat{\gamma }\cdot \widehat{F},\widehat{z},\widehat{N},\widehat{R})$, where $s>0$, the expected choices of potential violators are identical.
\end{proposition}

In light of Proposition \ref{DROPHD0}, a simple rescaling of the $\mu $ and $\sigma $ axes and a redefinition of the variables the axes represent allows us to reinterpret Figure \ref{fig:dropsurface}. The surface plot in the figure represents the DROP statistics for pairs $(\frac{\mu }{\gamma \cdot F},\frac{\sigma }{\gamma \cdot F})$, and the set identifies the subset of $CLIFFSET$ for any parameterization in which $N=100$ and $z=2$.

$DROP$ is increasing in both $N$ and $z$, and hence for larger values of either or both of these parameters the associated subset of $CLIFFSET$ is somewhat larger than the one identified in Figure \ref{fig:dropsurface}. To give a sense of how sensitive the subset is to these variables we have calculated DROP statistics for $z=3$ and for $N=160$ for three points on the outer boundary of $CLIFFSET^\prime$. 

\begin{table}[h]
\centering
\caption{DROP statistic for N=160 and z=3.}
\label{table: DROPNz}
\begin{tabular}{ccc}
\toprule
\multicolumn{1}{c}{\textbf{Threshold tuples  $\mathbf{(\frac{\mu}{\gamma \cdot F}, \frac{\sigma}{\gamma \cdot F})}$}} & \multicolumn{1}{c}{$\mathbf{DROP (N=160)}$} & \multicolumn{1}{c}{$\mathbf{DROP (z = 3)}$}  \\ 
\hline
\multicolumn{1}{c |}{(0.9125, 0.4088)} & 0.3781 & 0.3650 \\ 
\multicolumn{1}{c |}{(0.9981, 0.2993)} & 0.3870 & 0.3763 \\ 
\multicolumn{1}{c |}{(0.7825, 0.4856)} & 0.3968 & 0.3711 \\
\bottomrule
\end{tabular}
\end{table}

\subsection{Generalizing the modeling assumptions}\label{section:robustnessassumptions}
The geometry of the cliff is robust to a range of alternative assumptions in how $p$, $q$ and $g$ are determined. In Panels 1-3 of Figure \ref{fig:comparatives} we present comparative visualizations of the cliff for three illustrative variations in the assumptions that determine $p$, $q$ and $g$. In each case the solid line plots the baseline $E(v|R)$ cliff in Figure \ref{fig:thecliff}, and the dashed line plots the relevant comparison. In each case there is a distinct unruly regime, a compliant regime, a cliff-like transition between them, and a measured DROP exceeding our threshold.\footnote{\ In each case the value of DROP is sensitive to the assumed values of the parameters involved.}

In Panel 1 we alter  the technology of apprehension, which determines $p$. Two assumptions determine the objective probability of apprehension in the computable model: (i) exactly one unit of the enforcement resource, no more and no less, is required to investigate a violation; and, (ii) if a violation is investigated the violator is apprehended with probability $\gamma$. An alternative specification assumes that the probability that a violator is apprehended is continuous in the quantity of resources, $r$, devoted an investigation, $\gamma \cdot (1-1/\epsilon ^{r})$, where $\epsilon >1$ and $0<\gamma <1$. This probability is increasing and concave in $r$ and bounded above by $\gamma$. Given resources $R$, to maximize the expected number of apprehensions we must allocate $\frac{R}{v}$ to each violation. With this allocation, the probability that any violator is apprehended is $\gamma \cdot (1-\frac{1}{\epsilon^{\nicefrac{R}{v}}})$. In Panel 1 we assume $\epsilon=8, \gamma=0.8$

\begin{figure}[h]
    \centering
    \includegraphics[width=0.95\textwidth]{Penultimate Draft/Images/Fig 7.png}
    \caption{Robustness of the cliff, for variations in $P(v,R)$, $q$, and $\phi(g)$.}
    \label{fig:comparatives}
\end{figure}

In Panel 2 we alter the specifications determining peoples' subjective probabilities of apprehension. The computable model assumes that all people share the same mapping from the z-history to their subjective probability of apprehension, $q_{i}^{t}=Q(zh^{t}) \ \forall \ i$. Alternatively, we assume $q_{i}^{t}=Q(zh^{t})+\delta _{i}^{t}$, where the noise term, $\delta _{i}^{t}$, is a person specific random draw from the uniform distribution with support $[-0.2,+0.2]$. 

In Panel 3 we alter the process determining $g_{i}^{t}$. In the computable model, the $g_{i}^{t}$s are random draws from a normal distribution, $\phi (g_{i}^{t})$. Alternatively, the $g_{i}^{t}$s are drawn from a uniform distribution with support $[0.2,1.1]$.

\subsection{The space of passive policies} \label{section:passivepolicies}
Another issue concerns the robustness of our results to alternative levels of sanction. In the passive policies considered to this point the sanction F was fixed, and R varied. We now examine passive policies in which both R and F vary. We restrict our attention to passive policies for which $R \in \{0, ..., N \}$ and $F \in \{0, 0.05, ..., 1.45, 1.5\}$. This results in a set of passive polices with 101 · 31 = 3131 elements. Figure \ref{fig:policysurface} is a surface plot of $E(v|R, F )$ for this set of polices. Two cross-sections are highlighted. The first is produced by holding the baseline sanction constant and varying R. This reproduces the “cliff of positive feedback” from Figure \ref{fig:thecliff}. The second is produced by holding enforcement spending constant (see Section \ref{section:optpass}) and varying $F$. This produces another cliff of positive feedback.

\begin{figure}[h]
    \centering
    \includegraphics[width=0.95\textwidth]{figures/fig_43.png}
    \caption{The space of passive policies.}
    \label{fig:policysurface}
\end{figure}

Both cross-sections, as well as the surface itself, reveal the same regimes we saw in Figure \ref{fig:thecliff}: (i) the unruly regime consists of a plateau around the top edge of the figure; (ii) the compliant regime is the rectangular plateau at the bottom; and, (iii) the transitional regime is a cliff separating the two. The cliff of positive feedback is a pervasive feature of our dynamic process.

\section{Generality of deterrence policy analysis} \label{section:genODP}
In this concluding section we start by considering the generality and determinants of three central results for the baseline parameterization: (i) The reserve capacity rate with the optimal crackdown policy is higher than 60\%, meaning that in a typical period only 40\% of available enforcement resources are used to investigate violations; (ii) replacing the optimal passive policy with the optimal crackdown policy reduces the expected cost of ASB by 12\%; and, (iii) replacing the optimal crackdown policy with the optimal refined crackdown policy reduces the expected cost of ASB by less than 1\%. Naturally, these specific results, the reserve capacity rate and the measured reduction in costs, will differ from one parameterization to another. We show that they are driven primarily by $RAT$ and $DROP$.  

$RAT$ is effectively the ratio of the external cost imposed on society by a violation to the cost of investigating a violation, so it is intuitive that it is an important determinant of these three numbers. As we have seen, managing ASB means managing positive feedback. Since $DROP$ is a measure of the potential for disruptive feedback it is intuitive that it too is important to understanding these results.

\begin{figure}[h]
    \centering
    \includegraphics[width=0.95\textwidth]{Penultimate Draft/Images/Fig 16d.png}
    \caption{Reserve capacity and cost savings as a function of $DROP$.}
    \label{fig:dropscatter}
\end{figure}

In Figures \ref{fig:dropscatter} and \ref{fig:ratscatter} we report the reserve capacity rate and the two cost reduction percentages for  simulations in which $N$ and $z$ are fixed at their baseline values.\footnote{\ The estimates of optimal crackdown and refined crackdown policies reported in this section are generated using simpler versions of the directed search algorithms introduced in Section \ref{section:activepolicies}. The number of randomly seeded searches was reduced to 50 in each case implying that the errors in the estimation of the optimal policy will be slightly higher. But as we observed in Section \ref{section:activepolicies} the space in which the optimum occurs is very flat and so there is very little difference between estimates in the neighborhood of the optimum. Reducing the number of searches increased computational efficiency at nearly zero cost to estimation precision.} In these figures we are treating $DROP$ as something like a composite parameter. Given $N=100$ and $z=2$, as we saw in Section \ref{section:robustness}, $DROP$ is determined by $\frac{\mu}{\gamma \cdot F}$ and $\frac{\sigma}{\gamma \cdot F}$. In selecting the parameterizations used in constructing these figures, we used two criteria. First, that the distribution of $DROP$ be approximately uniform over the interval (0.275, 0.725). Second, subject to this constraint, that we include parameterizations distributed across the entire $(\frac{\mu}{\gamma \cdot F},\frac{\sigma}{\gamma \cdot F})$ parameter space. For each parameterization we calculated reserve capacity and the two cost reduction numbers. Then, in Figures \ref{fig:dropscatter} and \ref{fig:ratscatter}, we plot the values for our three central results for each parameterization with their values of $DROP$ and $RAT$. 

The first thing to notice is the magnitude of the correlation coefficients reported in the four panels of Figure \ref{fig:dropscatter}. The fact that they are relatively high is telling us that $DROP$ does function as a composite of the parameters that determine it, and the fact they are not $1$ is telling us that the effects of variation in these parameters are not completely captured by $DROP$. 

In Panels 1 and 3 of Figure \ref{fig:dropscatter}, we report the reserve capacity rate for the optimal crackdown policy for a number of DROP values in the interval (0.25, 0.75). In Panel 1, $RAT$ is fixed at 2 and in 3 it is fixed at 6. In both panels, reserve capacity is increasing in $DROP$, from a low of 20\% up to high of 70\%. In Panels 2 and 4 the circles represent the reduction in cost when the optimal passive policy is replaced with the optimal crackdown policy. In Panel 2 $RAT$ is 2, and in 4 it is 6. The cost reduction numbers are  positively associated with $DROP$ and  negatively associated with $RAT$. The crosses in those panels represent the reduction in cost when the optimal crackdown policy is replaced with the optimal refined crackdown policy. These cost reduction numbers are negligible and uncorrelated with $DROP$.

\begin{figure}[h]
    \centering
    \includegraphics[width=0.95\textwidth]{Penultimate Draft/Images/Fig 18d.png}
    \caption{Reserve capacity and cost savings as a function of $RAT$.}
    \label{fig:ratscatter}
\end{figure}

In Panel 1 of Figure \ref{fig:ratscatter}, we report the reserve capacity rate with the optimal crackdown policy as a function of $RAT$ for three $DROP$ values. For all three values, there is massive discontinuity in the relationship between reserve capacity and $RAT$. (Given the coarse grid of $RAT$ values, we have not identified any of the three discontinuities precisely.) When $RAT$ is very small, we are in the spitting on the sidewalk case and the best enforcement policy is to do nothing. As $RAT$ passes through the point of discontinuity the reserve capacity rate jumps up, and beyond the discontinuity it increases slowly with $RAT$. For the two larger $DROP$ values to the right of the discontinuity the rate exceeds 60\% and for the lowest value it is much lower. In Panel 2, we report the reduction in cost when the optimal passive policy is replaced with the optimal crackdown policy as a function of $RAT$ for the same three values of $DROP$. Again there is massive discontinuity in the relationship. When $RAT$ is small do nothing is the optimal enforcement policy regardless of whether we are using a passive or crackdown policy, so the cost reduction is 0. Beyond the discontinuity the cost reduction is inversely related to $RAT$.

The results reported in Figures \ref{fig:dropscatter} and \ref{fig:ratscatter} are based on parameterizations in which $N=100$ and $z=2$. Because simulations with larger values of these parameters are very time intensive, we have run just a few of them. Based on those simulations it seems clear that the patterns seen in these figures and the magnitudes of the numbers are robust to increases in the values of the parameters.\footnote{ \ Variations of $N=160$ and $z=3$ were run on two values each of $DROP\approx\{0.7, 0.4\}$ and $RAT=\{2, 10\}$. The pattern of cost savings for crackdown and refined crackdown was preserved. Costs savings for crackdown over passive increased slightly for both the larger $N$ and $z$ values, with the increases ranging from 0.1\% to 1.8\% over the baseline values with the maximum cost savings being 15.2\%. The cost savings for refined crackdown was again negligible with the total savings ranging from 0\% to 0.4\%. The reserve capacity values were almost identical for each for each of the four cases. For example, for $DROP\approx0.7$ and $RAT=2$, reserve capacity = (0.4260, 0.4240, 0.4302) for the $(N, z)$ tuples of (100, 2) (160, 2) and (100, 3) respectively.}


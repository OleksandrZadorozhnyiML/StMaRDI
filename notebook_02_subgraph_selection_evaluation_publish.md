Graph subset selection
================
David Reiffenscheidt and Oleksandr Zadorozhnyi

**Setup of the problem**

Given a "large" graph *G* = (*V*,*E*), where *V* represents the set of
nodes and *E* represents the set of edges, the goal is to select a
subgraph *S* = (*V*<sub>*s*</sub>,*E*<sub>*s*</sub>) from G such that S
is a meaningful and informative representation of the original graph
*G*. The subgraph selection problem involves finding an optimal or
near-optimal subgraph that satisfies certain criteria or objectives.

In the context of graphical modelling the problem of subgraph selection
corresponds to the problem of preserving the conditional independence
relationship in the subgraph *S* which has to be transferred from the
underlying graph *G*.

Furthermore, it is important to be able to evaluate the results of the
structure learning algorithms (learned on the data) taking as a ground
truth the selected subgraph.

``` r
### loading packages
library("bnlearn")
library("qgraph")
library("igraph")
```

    ## 
    ## Attaching package: 'igraph'

    ## The following objects are masked from 'package:bnlearn':
    ## 
    ##     as.igraph, compare, degree, subgraph

    ## The following objects are masked from 'package:stats':
    ## 
    ##     decompose, spectrum

    ## The following object is masked from 'package:base':
    ## 
    ##     union

``` r
data("alarm")
alarm_df <- as.data.frame(na.omit(alarm))

p = length(names(alarm))
n = dim(alarm)[1]
for (i in c(1:p)) {
  alarm_df[,i]<-as.numeric(alarm_df[,i])
}
```

Applying nonparanormal transformation to standardize the data.

``` r
library("huge")
alarm_df <- huge.npn(alarm_df)
```

Conducting the nonparanormal (npn) transformation via shrunkun ECDF....done.

``` r
#####
```

Defining “true” graph as proposed for the ALARM dataset in bnlearn

``` r
# "True" Graph ALARM
dag_alarm <- empty.graph(names(alarm))
modelstring(dag_alarm) <- paste0("[HIST|LVF][CVP|LVV][PCWP|LVV][HYP][LVV|HYP:LVF]",                    "[LVF][STKV|HYP:LVF][ERLO][HRBP|ERLO:HR][HREK|ERCA:HR][ERCA][HRSA|ERCA:HR]",
"[ANES][APL][TPR|APL][ECO2|ACO2:VLNG][KINK][MINV|INT:VLNG][FIO2]",            "[PVS|FIO2:VALV][SAO2|PVS:SHNT][PAP|PMB][PMB][SHNT|INT:PMB][INT]",              "[PRSS|INT:KINK:VTUB][DISC][MVS][VMCH|MVS][VTUB|DISC:VMCH]",                    "[VLNG|INT:KINK:VTUB][VALV|INT:VLNG][ACO2|VALV][CCHL|ACO2:ANES:SAO2:TPR]",
"[HR|CCHL][CO|HR:STKV][BP|CO:TPR]", sep = "")

qgraph(dag_alarm, legend.cex = 0.3,
       asize=2,edge.color="black", vsize= 4)
```

![](notebook_02_subgraph_selection_evaluation_publish_files/figure-markdown_github/unnamed-chunk-4-1.png)

Selection of the set of nodes for subsetting

``` r
### Subgraph 
subgraph_nodes <- c("INT","VALV","LVF","LVV","PCWP","HR","CCHL","CVP","HYP","HRSA","ERCA")
```

First procedure (bnlearns subgraph function). It is implemented through
a simple subsetting of the edges which are adjacent to the vertices
contained in the subgraph_nodes set.

``` r
procedure1 <- bnlearn::subgraph(dag_alarm, subgraph_nodes)
qgraph(procedure1, legend.cex = 0.3,
       asize=2,edge.color="black", vsize= 5)
```

![](notebook_02_subgraph_selection_evaluation_publish_files/figure-markdown_github/unnamed-chunk-6-1.png)

Second procedure selects the subgraph based on the following euristics.
Given the ground truth DAG *G* = (*V*,*E*) and subset of vertices
*V*<sub>*s*</sub> ⊂ *V* the goal is to find the corresponding set of
vertices *E*<sub>*s*</sub> such that for
(*V*<sub>*s*</sub>,*E*<sub>*s*</sub>) the structure of the distribution
for *P*<sub>*{V_*<sub>*s</sub>}</sub>* does not contradict the structure of the distribution
*P*<sub>*V*</sub> so that we the task of structure estimation (and then
benchmarking) on (*V*<sub>*s*</sub>,*E*<sub>*s*</sub>) can be done with
a new ground truth.

Let *G* = (*V*,*E*) be the original directed acyclic graph, and
*P*<sub>*V*</sub> is the joint distribution of random variables from
*V*. We take two and check whether these vertices are *d*−connected
given all others in *V*<sub>*s*</sub>. If they are not *d*-connected,
there is no association (no arrow in any direction). Otherwise, we have
correlation (with unknown direction). If one of the directions leads to
the cycle in the original graph, we resolve it and keep the other
direction, otherwise we keep both directions and then get the CPDAG.

``` r
######### Extract subgraph function

combn(subgraph_nodes, 2)[2,1]
```

    ## [1] "VALV"

``` r
dim(combn(subgraph_nodes,2))[2]
```

    ## [1] 55

``` r
extract_subgraph <- function(dag, nodes){
  sg <- bnlearn::subgraph(dag,nodes) # procedure 1 (to be discussed)
  combinations <- combn(nodes,2) # all combinations of 2 distinct nodes in "nodes"
  n <- dim(combinations)[2]
  for (i in 1:n){
    observed <- nodes[nodes!=combinations[1,i] & nodes!=combinations[2,i]] # V'\{v,w}
    if (!is.element(combinations[1,i], nbr(sg, combinations[2,i])) & # check if there exists an edge already
        !dsep(dag,combinations[1,i],combinations[2,i])){ ### check if d-connected
      sg <- set.edge(sg, from = combinations[1,i], to = combinations[2,i]) ### undirected edge in case d-connected
    }
  }
  return(cpdag(sg)) ### to be discussed: return(cpdag(sg))
}

procedure2 <- extract_subgraph(dag_alarm, subgraph_nodes)
qgraph(procedure2, legend.cex = 0.3,
       asize=2,edge.color="black", vsize= 5)
```

![](notebook_02_subgraph_selection_evaluation_publish_files/figure-markdown_github/unnamed-chunk-7-1.png)
Subsetting the dataset according to “subgraph_nodes” selection

``` r
alarm_dfSubset <-as.data.frame(alarm_df[,subgraph_nodes])
```

Applying constraint-based algorithms

``` r
Res_stable=pc.stable(alarm_dfSubset)

Res_iamb=iamb(alarm_dfSubset)

Res_gs=gs(alarm_dfSubset)

Res_fiamb=fast.iamb(alarm_dfSubset)

Res_mmpc=mmpc(alarm_dfSubset)
```

Applying score-based algorithms

``` r
Res_hc = hc(alarm_dfSubset)
Res_tabu = tabu(alarm_dfSubset)
```

Visualize resulting subgraph

``` r
par(mfrow=c(3,4))

ig_proc1 <- as.igraph(procedure1)
ig_proc2 <- as.igraph(procedure2)
ig_pc <- as.igraph(Res_stable)
ig_iamb <- as.igraph(Res_iamb)
ig_gs <- as.igraph(Res_gs)
ig_fiamb <- as.igraph(Res_fiamb)
ig_mmpc <- as.igraph(Res_mmpc)

ig_hc <- as.igraph(Res_hc)
ig_tabu <- as.igraph(Res_tabu)

plot(ig_proc1, main = "proc1", frame = T, layout=layout_with_fr, vertex.size=6,
     vertex.label.dist=0.5, vertex.color="red", edge.arrow.size=0.1, arrow.size=0.1, vertex.label.cex = 0.5)
plot(ig_proc2, main = "proc2", frame = T, layout=layout_with_fr, vertex.size=6,
     vertex.label.dist=0.5, vertex.color="red", edge.arrow.size=0.1, arrow.size=0.1, vertex.label.cex = 0.5)
plot(ig_pc, main = "pc", frame = T, layout=layout_with_fr, vertex.size=6,
     vertex.label.dist=0.5, vertex.color="red", edge.arrow.size=0.1, arrow.size=0.1, vertex.label.cex = 0.5)
plot(ig_iamb, main = "iamb", frame = T, layout=layout_with_fr, vertex.size=6,
     vertex.label.dist=0.5, vertex.color="red", edge.arrow.size=0.1, arrow.size=0.1, vertex.label.cex = 0.5)
plot(ig_gs, main = "gs", frame = T, layout=layout_with_fr, vertex.size=6,
     vertex.label.dist=0.5, vertex.color="red", edge.arrow.size=0.1, arrow.size=0.1, vertex.label.cex = 0.5)
plot(ig_fiamb, main = "fiamb", frame = T, layout=layout_with_fr, vertex.size=6,
     vertex.label.dist=0.5, vertex.color="red", edge.arrow.size=0.1, arrow.size=0.1, vertex.label.cex = 0.5)
plot(ig_mmpc, main = "mmpc", frame = T, layout=layout_with_fr, vertex.size=6,
     vertex.label.dist=0.5, vertex.color="red", edge.arrow.size=0.1, arrow.size=0.1, vertex.label.cex = 0.5)
plot(ig_hc, main = "hc", frame = T, layout=layout_with_fr, vertex.size=6,
     vertex.label.dist=0.5, vertex.color="red", edge.arrow.size=0.1, arrow.size=0.1, vertex.label.cex = 0.5)
plot(ig_tabu, main = "tabu", frame = T, layout=layout_with_fr, vertex.size=6,
     vertex.label.dist=0.5, vertex.color="red", edge.arrow.size=0.1, arrow.size=0.1, vertex.label.cex = 0.5)
```

<img src="notebook_02_subgraph_selection_evaluation_publish_files/figure-markdown_github/unnamed-chunk-11-1.png" width="20%" />

Selfdefined function for different measures

``` r
measure = function(estim, true){
  result <- matrix(1,4)
  com <- bnlearn::compare(estim, true)
  shd <- shd(estim,true)
  result[1,] <- com$tp
  result[2,] <- com$fp
  result[3,] <- com$fn
  result[4,] <- shd
  rownames(result) <- c("true positives","false positives","false negatives","structural hamming distance")
  colnames(result) <- deparse(substitute(estim))
  return(result)
}
```

Metric evaluation for procedure 1

``` r
measure(Res_stable, procedure1)
```

    ##                             Res_stable
    ## true positives                       6
    ## false positives                      2
    ## false negatives                      5
    ## structural hamming distance          5

``` r
measure(Res_iamb, procedure1)
```

    ##                             Res_iamb
    ## true positives                     5
    ## false positives                    3
    ## false negatives                    6
    ## structural hamming distance        4

``` r
measure(Res_gs, procedure1)
```

    ##                             Res_gs
    ## true positives                   5
    ## false positives                  3
    ## false negatives                  6
    ## structural hamming distance      4

``` r
measure(Res_fiamb, procedure1)
```

    ##                             Res_fiamb
    ## true positives                      5
    ## false positives                     3
    ## false negatives                     6
    ## structural hamming distance         4

``` r
measure(Res_mmpc, procedure1)
```

    ##                             Res_mmpc
    ## true positives                     0
    ## false positives                    8
    ## false negatives                   11
    ## structural hamming distance        9

``` r
measure(Res_hc, procedure1)
```

    ##                             Res_hc
    ## true positives                   5
    ## false positives                  3
    ## false negatives                  7
    ## structural hamming distance      9

``` r
measure(Res_tabu, procedure1)
```

    ##                             Res_tabu
    ## true positives                     5
    ## false positives                    3
    ## false negatives                    7
    ## structural hamming distance        9

Metric evaluation for procedure 2

``` r
measure(Res_stable, procedure2)
```

    ##                             Res_stable
    ## true positives                       4
    ## false positives                     16
    ## false negatives                      7
    ## structural hamming distance         16

``` r
measure(Res_iamb, procedure2)
```

    ##                             Res_iamb
    ## true positives                     6
    ## false positives                   14
    ## false negatives                    5
    ## structural hamming distance       14

``` r
measure(Res_gs, procedure2)
```

    ##                             Res_gs
    ## true positives                   6
    ## false positives                 14
    ## false negatives                  5
    ## structural hamming distance     14

``` r
measure(Res_fiamb, procedure2)
```

    ##                             Res_fiamb
    ## true positives                      6
    ## false positives                    14
    ## false negatives                     5
    ## structural hamming distance        14

``` r
measure(Res_mmpc, procedure2)
```

    ##                             Res_mmpc
    ## true positives                     4
    ## false positives                   16
    ## false negatives                    7
    ## structural hamming distance       16

``` r
measure(Res_hc, procedure2)
```

    ##                             Res_hc
    ## true positives                   3
    ## false positives                 17
    ## false negatives                  9
    ## structural hamming distance     17

``` r
measure(Res_tabu, procedure2)
```

    ##                             Res_tabu
    ## true positives                     5
    ## false positives                   15
    ## false negatives                    7
    ## structural hamming distance       17

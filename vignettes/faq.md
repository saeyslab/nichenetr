FAQ NicheNet
================
Robin Browaeys
2020-05-20

<!-- github markdown built using
rmarkdown::render("vignettes/faq.Rmd", output_format = "github_document")
-->

This document tries to give an extensive answer to some questions we got
in the past months.

## When going to the vignettes, I see you require differential expression between two conditions. What if I only have steady-state data and want to use NicheNet? Can’t NicheNet be used to find ligand-receptor pairs in steady-state?

NicheNet is a tool that let you study how ligands affect gene expression
in putatively neighboring/interacting cells, and prioritize the most
important ligands based on their effect (in other words: prioritizing
expressed ligand-receptor interactions based on their observed target
genes). But to do this you need to have data about this effect in gene
expression you want to study. So, there need to be ‘some kind of’
differential expression in a receiver cell population, caused by ligands
from one of more interacting sender cell populations. Concretely, this
can be differential expression in one cell type between two conditions
caused by cell-cell interactions, or this can also be differential
expression between two cell types if for example these are a progenitor
and differentiated cell type where the differentiation is influenced by
the microenvironment. In that case you don’t necessarily need two
“conditions”.

If you would just be interested in the ligand-receptor interactions in
homeostatic conditions (and you only have homeostatic condition data),
you can always just infer ligand-receptor interactions based on their
expression level. This is the end goal in some other tools like
CellphoneDB (see also further question), or an intermediary step in the
NicheNet pipeline: (cf step 3 in the Seurat steps vignette
<https://github.com/saeyslab/nichenetr/blob/master/vignettes/seurat_steps.md>).

However, prioritizing expressed ligand-receptor interactions based on
observed target genes by applying NicheNet on steady-state data will be
harder because you don’t have a clear set of genes affected by the
intercellular communication process only. You could still do it, by
using the total set of genes in the genome as the background, and all
genes expressed in your cell type as gene set of interest (so you need
to change the geneset\_oi to all expressed genes, and the background to
all genes in the ligand-target matrix). But I don’t really recommend
that because not all expressed genes will be influenced by cell-cell
interactions, on the contrary most genes will probably be expressed
because they are part of general and cell-type specific gene expression
programs. And if you can’t discriminate between cell-intrinsic and
cell-extrinsic effects in your data, and you will still use NicheNet,
then you risk that ligands will be falsely linked to some
‘cell-intrinsic’ genes. And this will confound the analysis and lead
to irrelevant and noisy predictions.

So the take-away message is: NicheNet predicts upstream ligand-receptor
pairs that might regulate genes in your gene set of interest. The more
this gene set of interest consists of genes that you expect to be
regulated by the extracellular microenvironment, the better. If this
gene set would also consist of genes not necessarily influenced by the
environment (which is the case if using all expressed or cell-type
specific genes as gene set of interest), NicheNet can still provide good
predictions sometimes, but only if a substantial part of these genes
would affected by cell-cell interactions. Therefore, be cautious when
doing that and remember: garbage in, garbage out.

## As gene set of interest, I am using the Seurat cluster markers found via FindMarkers function. Is this a good idea?

As discussed in the question above, this is not a good idea in general.
If you use the marker genes of a (receiver) cluster as gene set of
interest, you assume that the cluster-specific genes are induced by
cell-cell interactions with other cell types. However, I think in most
cases, many cluster-specific genes reflect a cell-type specific gene
expression program that is part of the ‘cell-instrinsic’ gene expression
program, and not entirely influenced by environmental factors. Or in
other words, most genes of the cluster markers are probably not really
induced/regulated by the ligands from the interacting cells, and
therefore using NicheNet on this gene set is not ideal. Using NicheNet
in this case will increase the chance of falsely linking ligands to some
‘cell-intrinsic’ genes.

There might be some exceptions for which using cluster-specific genes as
gene set of interest might work though. Maybe you have 2 clusters of the
same cell type, but one cluster is located differently in the tissue,
and thus has different extracellular influences: then DE genes between
these clusters do not reflect differences in cell type, but in cell-cell
interactions affecting that cell type. In that case, you can use
NicheNet. Or another example possibility is: a progenitor cell cluster
differentiates into a another cell cluster under influence of cell-cell
interactions. Then you can use NicheNet on the cluster-specific genes of
that differentiated cell population to check influence of cell-cell
interactions on this differentiation process.

## If I compare the prioritized ligand-receptor interactions of NicheNet with those of CellphoneDB I find little overlap. How can that be explained?

A main reason is that both tools have a different goal in mind to
prioritize ligand-receptor interactions. The goal of CellphoneDB is to
find cell-type specific ligand-receptor pairs. So CellphoneDB is ideal
if you want to investigate the differential interaction potential of
different cell types. CellphoneDB first looks for all expressed
ligand-receptor pairs between interacting cells. Then it finds cell-type
specific ligand-receptor pairs by looking at the expression of both the
ligand and the receptor. The stronger and more specific
ligands/receptors are expressed in the sender/receiver cell populations,
the better it will be prioritized via a CellphoneDB analysis.

The goal of NicheNet is complementary to that of CellphoneDB. NicheNet
wants to find ligand-receptor pairs that are most likely to regulate
gene expression in the receiver cell. So it looks for ligand-receptor
pairs for which evidence for signaling interactions exist in the data.
In contrast to CellphoneDB this might give more functional information
(you have an idea about the downstream targets of a ligand-receptor
interaction), but also some clues about which interactions might really
be active. This because expression of the ligand and the receptor at RNA
level (as found by CellphoneDB), does not necessarily mean that they
interact in reality. If downstream signaling effects of such an
interaction are observed, e.g. via NicheNet, there is a bit more
evidence that this interaction might be happening in reality. So just
like CellphoneDB, NicheNet starts to look for all expressed
ligand-receptor pairs between the interacting cells of interest. But
instead of prioritizing these by looking at expression strength,
NicheNet will prioritize them based on observed target genes in the
receiver cell. A possible disadvantage of this approach is that some
ligand-receptor pairs with relatively low expression can be returned.
Therefore we also recommend to check expression of the ligands and their
receptors after the NicheNet prioritization. To give an example: If a
ligand ranked 7th out of 200 would be much stronger expressed than a
ligand that would be 1st or 2nd ranked according to the ligand activity,
this more strongly expressed candidate ligand might be more
interesting\!

So when comparing NicheNet to CellphoneDB output following differences
can be expected to be present. Ligad-receptor pairs picked up by
NicheNet but not by CellphoneDB might be or generally expressed or
rather lowly expressed, but have some ‘signaling evidence’. Pairs picked
by CellphoneDB but not by NicheNet will be strongly and cell-type
specific expressed, but there is no evidence based on prior knowledge on
signaling pathways that these pairs will actually be functional. Of
course, a lack of evidence does not mean there is no signaling is going
on.

## When I check the expression of some top ligands according to NicheNet, I see that some of these ligands and/or their receptor are really lowly expressed? How is that possible? And does NicheNet not take into account expression data of the cells?

The prioritization of ligands by NicheNet (ligand activity analysis)
will only occur based on enrichment of their target genes in the set of
genes that are differentially expressed in the receiver cell. So there
is no prioritization based on the strength of expression of the ligand
in the sender cell or strength of expression of the receptor(s) in the
receiver cell. Expression in sender cells is only used to determine
which ligands are expressed in a sender cell, and expression in receiver
cells is used to determine which receptors are expressed in the receiver
cell. The default definition of ‘being expressed’ is that a gene should
be expressed in 10% of cells in the cluster of interest. This is not so
high (you can put a more stringent cutoff if you want), resulting in the
possible outcome that a ligand, top-ranked according to the enrichment
of its target genes, is actually not very highly expressed. So what you
observe, can be expected based on how NicheNet prioritizes ligands.

In the current version of NicheNet, expression strength is thus not
directly included because we find it hard to formalize the tradeoff
between ligand activity and expression. But because expression level is
important, we recommend to check expression of the ligands and their
receptors after the NicheNet prioritization. To give an example: If a
ligand ranked 7th out of 200 would be much stronger expressed than a
ligand that would be 1st or 2nd ranked according to the ligand activity,
this more strongly expressed candidate ligand might be more
interesting\!

## I observe that some top ligands according to NicheNet to be inducing in my condition of interest compared to control are higher expressed in control than condition of interest. How is this possible?

The ligand prioritization is currently done without taking into account
to expression value of the ligand, the only condition is that the ligand
is expressed. This means that the algorithm does not consider whether a
ligand is more strongly expressed in the case vs control, as long as it
is expressed at a sufficient level. Of course, this might be interesting
information to use for further prioritization after the NicheNet ligand
activity analysis\! Therefore we recommend checking the expression of
your ligands after NicheNet prioritization. The most interesting hits
might be the ones where ligand and/or receptor is upregulated in the
condition of interest compared to the control. The opposite pattern
could also be interesting if you could consider following hypothesis.
Some ligands might downregulate many genes (instead of, or in addition
to, upregulate many genes). This means that if ligand X is active in the
control condition, it can downregulate genes there. When ligand X itself
is repressed in the case condition, and has a lower expression, its
downregulatory effect might disappear, leading to upregulation of genes
that are a downregulated target gene of the ligand X, and high ligand
activity of ligand X according to NicheNet.

If you would only want to consider ligands upregulated in the
case-vs-control condition, you can always change the pipeline by first
performing a differential expression analysis on your sender cells, and
considering only upregulated ligands, and not all expressed ligands as
‘potential ligands’ for the NicheNet analysis. But I would recommend
only filtering on upregulation afterwards, because ligands can be
differentially active due to other processes than upregulation at the
RNA level.

## Can I use NicheNet if I want to find which potential ligands were important in cells of interest, even though I don’t have expression data of the possible interacting/sender cells?

It is perfectly possible to apply NicheNet in case you don’t have a
‘sender cell type’. In that case, your analysis is more of an
‘upstream regulator analysis’ than an ‘intercellular communication
analysis’. You can do this by looking at all ligands in the NicheNet
database (or all ligands for which a potential receptor is expressed)
instead of only looking ligands that are expressed by a specific sender
cell type. Only a small change in the code is required for this by using
`ligands` instead of `expressed_ligands` in the basic and Seurat steps
vignette for defining `potential_ligands`. If using the Seurat wrapper
function you can set the parameter function of sender to `sender =
“undefined”`.

You just need to be sure that you can use the expression data of the
responder/receiver cell type to define a set of genes that are likely
affected by the upstream regulatory intercellular signaling process.

In case you don’t have data about the sender cell type(s) of interest,
but you know there is some publicly available expression data available
of that cell type or the tissue of interest, it might be helpful to use
this dataset as proxy for the expression of the ligands.

## Is the NicheNet prior model based on human or mouse data or is there a separate model for human and mouse?

The model is primarily build based on human data sources. However, we
also included some mouse data sources into the general model (at least
for genes that have a mouse one-to-one ortholog), but these were
strongly in the minority. So the NicheNet model is tailored more for
human than mouse, and the mouse model that we use in some of the
vignettes is just obtained by converting human gene symbols to their
mouse one-to-one orthologs.

Because this model is indeed tailored more for human, you can expect
that it would perform a bit better for human than mouse, but I don’t
think the difference is very large. In our paper, we evaluated the
ligand-target predictions based on \>100 ligand treatment datasets of
both human and mouse origin, and there is no indication on these
datasets that human would work better than mouse.

However, there is a disadvantage: for some mouse-specific genes, we
don’t have information in the NicheNet ligand-target model. This
implies that NicheNet might miss some patterns if well-known mouse-only
target genes would be DE in the dataset you are working on. A solution
to this would be to construct a mouse model with additional mouse data
sources yourself. This is some work, but is not so difficult. In the
vignette
<https://github.com/saeyslab/nichenetr/blob/master/vignettes/model_construction.md>,
we show how to build a model yourself.

## Can I use NicheNet on data other than mouse or human data?

Yes, there are two options here.

First, you can explore the possibility of constructing an
organism-specific model yourself with only data sources of your organism
of interest. This is some work, but is not so difficult. In the vignette
<https://github.com/saeyslab/nichenetr/blob/master/vignettes/model_construction.md>,
we show how to build a model yourself. The main thing to consider
though, is that you should have enough primary data sources of the
organism to do this.

The second option would be to use the current human NicheNet model and
convert the human gene symbols to their one-to-one orthologs of the
organism of interest. The main thing to consider here is that the
organism of interest should not be too dissimilar from human.

So, you should make the tradeoff between the homology between species,
and the number and quality of data sources that are available to build a
new species-specific model. For example, for primates, you could use the
human model and convert gene symbols, but for Drosophila making a new
model from own data sources is probably more appropriate.

## I decided that I want to construct a model with my own data sources. Do I really need to optimize the data source weights?

As shown in our paper, data source weight optimization does improve the
performance the model a bit, but this effect is not large. A model
without data source weight optimization also seemed to work well.
Because the optimization is a difficult and time-consuming procedure, we
don’t see a problem in not optimizing your data source weights, but only
if you are confident that your data sources are of high quality and not
too noisy\!

## Can NicheNet say which cell populations interact with each other and which don’t?

No, NicheNet can’t give you a clear-cut answer to that question. But, it
can suggest this. If for many ligands from a certain cell types, many
target genes are DE in the receiver cell type, this might suggest that
some intercellular signaling processes are going on between both cell
types. In this vignette:
<https://github.com/saeyslab/nichenetr/blob/master/vignettes/target_prediction_evaluation_geneset.md>,
we show how you can calculate what fraction of genes DE in the receiver
cell type, might be a target of the top prioritized ligands of NicheNet.
This way, we found that in one of the case studies described in the
paper 50% of the DE genes were a strongly predicted target gene of the
prioritized ligands. Such a high fraction suggests that intercellular
signaling might indeed be going on, but this should of course be still
validated.

## Can I use NicheNet as a gene set based approach, e.g. to extract a gene set consisting of target genes of ligand X?

Yes, you can use the functions `extract_top_n_targets()` and
`extract_top_fraction_targets()` for this. In the near future, we will
also add a vignette that can be used to score individual cells based on
their expression of the top n target genes of a ligand of interest.

## How can I generate a double layer ligand-receptor-target circos plot as shown in the Immunity paper (<https://www.cell.com/immunity/fulltext/S1074-7613(19)30368-1>)?

As far as we know, there is no straightforward way to directly make this
kind of “ligand-receptor-target” circos plot via the circlize R package.
Therefore we made this “ligand-receptor-target” circos plot by making
first two separate circos plots: the ligand-target and ligand-receptor
circos plot. We made sure that the ligand-receptor circos plot was
bigger than the ligand-target plot. Then we overlayed them in Inkscape
(with the center of the two circles at the same location) and removed
everything from the ligand-receptor circos plot except the outer
receptor layer. Making the ligand-receptor circos plot is very similar
to making the ligand-target circos plot, as shown in the circos plot
vignette.

If you would want to split up target genes and receptors in different
groups according to signaling pathway (as done in that paper), then you
first need to define these groups in a specific data frame in advance
(cf what is shown for ligands in the `ligand_type_indication_df`in the
vignette). When you then want to overlay receptors in this case, you
need to make sure that the ligand-receptor weights of receptors in one
group are proportional to the ligand-target weights of the targets in
that group (to generate the nice overlay effect). So in that case, the
ligand-receptor weights are proportional to the ‘underlying’
ligand-target regulatory potential scores and not reflective of prior
information supporting the specific ligand-receptor interaction (as
shown in the current vignette for ligand-receptor circos plots).

We see that this is indeed a cumbersome and far from ideal situation.
Ideally we would be working on a more straightforward solution in the
future, but we can’t promise anything. This is not so high priority for
us because we think that the combined ligand-activity, ligand-receptor
and ligand-target heatmaps are (at least) equally informative.

## Can I use my own ligand-receptor network instead of the one included in the NicheNet model?

Yes, you definitely can\! In this vignette
<https://github.com/saeyslab/nichenetr/blob/master/vignettes/model_construction.md>,
you can see how to incorporate your new data source with new
ligand-receptor pairs and build a customized NicheNet prior model. I
would suggest you just add your new ligand-receptor data source to the
existing ones (or alternatively: remove the existing ones), and give
your new data source as weight 1, and use the optimized weights for the
other data sources.

## Will you update the data sources behind NicheNet and the NicheNet model?

We are indeed planning to update the prior model once by incorporating
some new and updated databases, but we can’t pin a date on this now.

## When interpreting the output ligand-target matrix, I see that not all my genes of my gene set of interest are in the final ligand-target matrix? Why is this?

The reason for this is that the ligand-target matrix only shows gene
that are a top predicted target of at least one of the top ligands. You
can tweak some of the parameters of the function to be more permissive
if you want to include more genes (by allowing lower regulatory
potential scores / setting the n of top n targets higher).

## Can I use NicheNet on bulk RNAseq data as well, or only on single-cell data?

Yes, you can if you are working on cell population – sorted bulk data.
Because of the better sequencing depth, this might even work better (of
course depending on the setting and whether you have clean data of the
cell populations of interest). We already applied NicheNet on bulk
RNA-seq data, as shown in this paper:
<https://www.cell.com/immunity/fulltext/S1074-7613(19)30368-1> .

We did this analysis by following the steps explained in the basic
vignette:
<https://github.com/saeyslab/nichenetr/blob/master/vignettes/ligand_activity_geneset.md>;
Some adaptations to the basic vignette are needed, though. For example,
you need another way for defining the set of expressed genes in sender
and receiver, and another way of doing the DE analysis for defining the
gene set of interest. To define which genes are expressed you could use
`filterByExpr` from the package `edgeR` ; or check the distribution of
expression values and choosing a cutoff based on this distribution. If
you doubt about choosing the cutoff, I would recommend to consider more
genes expressed than removing some genes from the analysis. The number
of expressed genes differs per cell type and organism, but should be in
the range of 10000-15000 for bulk rna-seq data (so not 1000 or 2000
genes).

In the Materials and Methods of this paper, you can see as example what
we did to determine expressed genes and perform the DE analysis.

## Does NicheNet return all significant ligand-receptor interactions?

No, you should realize that NicheNet does not give a clear-cut answer
about which ligand-receptor pairs are significant and which not. First,
in the NicheNet pipeline, we determine all possible expressed
ligand-receptor pairs. Then, a ligand activity analysis is performed to
rank ligands based on their activity (this accords to the enrichment of
strongly predicted target genes of a ligand in the gene set of interest
in the receiver cell). Based on that, the top n (default 20) ligands are
selected for further analyses and visualizations.

## If I look at the ligand-receptor network, I see some low interaction scores between highly expressed ligand-receptor interactions, or some high scores between lowly expressed pairs? What is the reason for this?

The ligand-receptor interaction scores shown in the ligand-receptor
matrix/heatmap are a proxy for the confidence that this ligand interact
with that receptor based on prior knowledge. The ligand-receptor network
at the basis of NicheNet is the result of integrating multiple single
ligand-receptor data sources. The score you see on this plot, the prior
ligand-receptor interaction potential, is a proxy of the number and
‘quality’ of the ligand-receptor and PPI data sources that report a
certain ligand-receptor interaction. If a ligand-receptor interaction is
described in many, confident databases, this score will be higher than
if the interaction is predicted via just one protein-protein interaction
database.

This score is thus solely based on prior knowledge and not based on the
expression in your dataset\! The only way expression is used, is to
define whether a ligand/receptor is expressed yes or no. So
ligands/receptor that are not expressed at all won’t appear in this
heatmap, but it is indeed possible that some lowly expressed
ligand-receptor pairs will be shown here.

In the near future, we will provide code to show how to make this
heatmap with scores that reflect the expression strength of both ligand
and receptor (product of expression / average of expression).

## In the Seurat vignette, you show some code to define bona fide ligand-receptor interactions. But I see that some ligand-receptor pairs that are not considered bona fide have higher prior interaction scores than some bona fide pairs. How is this possible if the prior interaction score is based on the confidence of the prior information? Should bona fide pairs not have the highest scores?

We can indeed see why this confusing\! We categorize bona-fide ligands
based on whether they were coming from a validated curated
ligand-receptor database (like KEGG, guide2Pharmacology, …) or not. This
because we have also ligand-receptor interactions not present in these
databases, but predicted based on protein-protein interaction networks:
e.g. if based on annotation we know that gene X encodes a ligand, and
gene Y a receptor, we predict a ligand-receptor interaction between X
and Y if there is a PPI database (but not ligand-receptor database)
supporting X-Y interaction. Because they are predicted as
ligand-receptor interaction (and not part of a curated ligand-receptor
database), they are not bona fide. But they can have higher scores than
some bona fide ones, if there are many confident PPI and pathway
databases reporting this interaction. So basically: a bona fide
ligand-receptor interaction is known as ligand-receptor interaction, a
non-bona fide interaction is a PPI (with much evidence if high scores)
but not yet documented as LR interaction in a curated database. Focusing
on bona fide LR interactions will give you well-validated pairs, whereas
including non-bona fide ones can be more novel, but are less confident.

## Although NicheNet already prioritizes ligand-receptor pairs strongly, I still find there are too many possible pairs to experimentally validate. What types of information would you recommend to consider for even further prioritization?

I think you could take many factors into consideration. First, you can
consider the NicheNet ligand activity (as you already did). But instead
of having a strong preference for the 1st ranked ligand versus the 10th
ligand, I suggest also taking into account the expression level of both
ligand and receptor. For this you can look at the expression value, the
fraction of cells expressing the ligand/receptor and whether the
ligand/receptor is cell-type specific or not. In the ideal case, you
would have case-vs-control data and you could also look at
ligands/receptor pairs for which the ligand and/or receptor is
upregulated in the case vs control condition. Next to this, you could
check whether the ligand-receptor interaction is ‘bona fide’ or not (cf
previous question) and whether the ligand-receptor prior interaction
potential weight is high or not. But I suggest giving the highest weight
to the ligand activity and expression/upregulation information\!

## How should I interpret the pearson correlation values of the ligand activity? What does it mean when I have negative values?

In the ligand\_pearson\_matrix, you find the pearson correlation values
between NicheNet prior ligand-target predictions and the gene set of
interest. A high score means that top predicted target genes of a ligand
are enriched in the gene set of interest compared to the background. We
would say that a good Pearson score for top ligands is around 0.10 and
higher. If you have low non-zero scores (such as 0.04), it might still
be that the ranking of ligands is valuable, although this ranking will
only be based on a few top predicted target genes of a ligand that are
in the gene set of interest. Scores around zero (such as 0.008) mean
that top predicted targets of ligands are not enriched in the gene set
of interest compared to a background, thus that ranking of ligands is
not valuable in that case. These ligands will typically have an AUROC
around 0.50. This means that NicheNet does not find evidence that the
sender cell is regulating gene expression changes in the receiver cell
(although the lack of evidence does not demonstrate the lack of).

If these values are negative, this means that there is an enrichment of
the top predicted target genes in the background compared to the gene
set of interest. When these values are negative, but very close to zero,
then there is nothing to worry about. But when these values are really
negative, this can indicate that something went wrong in the analysis,
probably in the definition of background and gene set of interest.
Because the gene set of interest is typically much smaller than the
background, a stronger anti-correlation score for a ligand would only
occur when the gene set of interest would almost entirely consist of the
genes with the lowest regulatory potential scores for a ligand. And this
would only rarely occur. One possible explanation for finding strongly
negative scores, could be that your background gene set is too small,
and thus not a good representation of the genomic background. Your gene
set of interest should never be larger than the background. We recommend
to have at least around 5000-6000 genes in your background. If you would
still have issues, or not enough expressed genes to have a large
background, you could use all genes present in the ligand-target matrix
as genomic background.

## If I look at the combined ligand-activity-target heatmap, I see that some ligands seem to have higher activity, but less target genes than other ligands. How is this possible if the ligand activity is based on the observed target genes in the gene set of interest?

The ranking of the ligands is based on how well top-ranked target genes
of a ligand (based on regulatory potential) are enriched in the gene set
of interest (= Pearson correlation between whether a gene belongs to the
gene set of interest and the regulatory potential scores of target genes
of a ligand). For highly ranked ligands, this means that many top-ranked
target genes are in the gene set of interest compared to the background.
Note that the Pearson correlation to rank ligands is calculated for each
ligand separately. As a result, this correlation only depends on the
ranking of target genes for one ligand ( ‘relative’ ligand-target
regulatory potential scores for that ligand, without looking at other
ligands).

The regulatory potential scores of ligand-target links that are shown in
the heatmap visualize the ‘absolute’ ligand-target regulatory potential.
This is a confidence score related to how many data sources confirm the
regulatory interaction between a ligand and a target. A higher score
means that there is more evidence for the specific ligand-target
interaction. It is important to note that in this heatmap, you only show
the genes that are part of your gene set as possible target genes.

The reason why you can see a contradiction between the Pearson
correlation and regulatory potential scores is thus the following:
Ligand X is ranked higher than Ligand Y, because its top-ranked target
genes are more enriched than the target genes of Ligand Y, although the
absolute regulatory potential scores of the Ligand Y-target links are
higher (this can be the case for very well-studied ligands like TNF,
IFNG, IL1B, …) . Ligand Y targets are however probably less enriched
because there will might be more high scoring target genes of Ligand Y
that are not in the gene set of interest compared to Ligand X targets.
That can explain why there is stronger enrichment for Ligand X targets,
although there is more information for Ligand Y.

## Can NicheNet also calculate ligand-receptor interactions that specifically downregulates target genes? Or can NicheNet distinguish between target genes that are upregulated or downregulated due to a ligand-receptor interaction?

NicheNet can indeed work with both upregulated and downregulated target
genes. In default cases, I would recommend considering both upregulated
and downregulated genes as potential target genes in the gene set of
interest, because some ligand-receptor interactions will indeed lead to
both up- and downregulation of genes as you mention. But it is perfectly
possible to look at down- or upregulated genes only. Depending on the
research question, this could sometimes be interesting. It could also be
interesting to find which ligands are more responsible for upregulation
and which more for downregulation.

This might be interesting because the NicheNet prioritization algorithm
does not directly distinguish between up- and downregulated genes. This
has an important consequence to keep in mind: if there are more
upregulated genes than downregulated genes, ligands that only upregulate
genes will be more likely active than downregulatory genes.

NicheNet cannot directly distinguish between up- and downregulated genes
because the regulatory potential accords to the
evidence/potential/probability that a ligand might regulate the
expression of a target gene. And regulate can mean both upregulation and
downregulation. It’s not that NicheNet knows a priori which target genes
will be upregulated and which downregulated. This is because of the
following reason: during construction of the model (so the calculate the
prior regulatory potential scores), we considered all retrieved
interactions as active because most databases at the basis of NicheNet
don’t provide information about whether an interaction is inducing or
repressing. So the only information we have in the final model is: there
is no evidence for the regulation of target X by ligand Y (reg potential
of 0) or there is evidence for the regulation of target X by ligand Y
(reg potential \> 0; with more evidence leading to higher scores; but a
repressive regulatory interaction leading to downregulation will also
have a score higher than 0).

What you can do as a user is visualizing the expression of the target
genes (in a combined heatmap for example, as shown in the basic
vignette). This will allow you to check which specific genes are up- and
downregulated.

## If I have two main cell types, divided in two subpopulations based on the condition (case-vs-control): what do I need to use as receiver and sender cell type of interest in e.g. the Seurat vignette?

I recommend using all cells of the receiver cell type of interest as
receiver (so including both conditions). For the senders, it might
sometimes be better to use only the cells in the case condition.
However, I think it is perfectly fine to include all sender cells as
well (thus from both conditions) and check after the NicheNet analyses
which ligands might be upregulated in the case condition compared to
control.

## My question is not in these list? What should I do now?

First, you can check the open and closed issues
(<https://github.com/saeyslab/nichenetr/issues>) of this package on
github to see whether your question might be addressed in one of these.
If not, don’t hesitate to open a new issue. If you would prefer to keep
the discussion private, you can also send me an email
(<robin.browaeys@ugent.be>), but I prefer that you open an issue so
other users can learn from it as well\!

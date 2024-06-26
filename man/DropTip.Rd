% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/DropTip.R
\name{DropTip}
\alias{DropTip}
\alias{DropTip.phylo}
\alias{DropTip.Splits}
\alias{DropTipPhylo}
\alias{DropTip.multiPhylo}
\alias{DropTip.list}
\alias{DropTip.NULL}
\alias{KeepTipPreorder}
\alias{KeepTipPostorder}
\alias{KeepTip}
\title{Drop leaves from tree}
\usage{
DropTip(tree, tip, preorder = TRUE, check = TRUE)

\method{DropTip}{phylo}(tree, tip, preorder = TRUE, check = TRUE)

\method{DropTip}{Splits}(tree, tip, preorder, check = TRUE)

DropTipPhylo(tree, tip, preorder = TRUE, check = TRUE)

\method{DropTip}{multiPhylo}(tree, tip, preorder = TRUE, check = TRUE)

\method{DropTip}{list}(tree, tip, preorder = TRUE, check = TRUE)

\method{DropTip}{`NULL`}(tree, tip, preorder = TRUE, check = TRUE)

KeepTipPreorder(tree, tip)

KeepTipPostorder(tree, tip)

KeepTip(tree, tip, preorder = TRUE, check = TRUE)
}
\arguments{
\item{tree}{A tree of class \code{\link[ape:read.tree]{phylo}}.}

\item{tip}{Character vector specifying labels of leaves in tree to be
dropped, or integer vector specifying the indices of leaves to be dropped.
Specifying the index of an internal node will drop all descendants of that
node.}

\item{preorder}{Logical specifying whether to \link{Preorder} \code{tree} before
dropping tips.  Specifying \code{FALSE} saves a little time, but will result in
undefined behaviour if \code{tree} is not in preorder.}

\item{check}{Logical specifying whether to check validity of \code{tip}. If
\code{FALSE} and \code{tip} contains entries that do not correspond to leaves of the
tree, undefined behaviour may occur.}
}
\value{
\code{DropTip()} returns a tree of class \code{phylo}, with the requested
leaves removed. The edges of the tree will be numbered in preorder,
but their sequence may not conform to the conventions of \code{\link[=Preorder]{Preorder()}}.

\code{KeepTip()} returns \code{tree} with all leaves not in \code{tip} removed,
in preorder.
}
\description{
\code{DropTip()} removes specified leaves from a phylogenetic tree, collapsing
incident branches.
}
\details{
This function differs from \code{\link[ape:drop.tip]{ape::drop.tip()}}, which roots unrooted trees,
and which can crash when trees' internal numbering follows unexpected schema.
}
\section{Functions}{
\itemize{
\item \code{DropTipPhylo()}: Direct call to \code{DropTip.phylo()}, to avoid overhead of
querying object's class.

\item \code{KeepTipPreorder()}: Faster version with no checks.
Does not retain labels or edge weights.
Edges must be listed in preorder.
May crash if improper input is specified.

\item \code{KeepTipPostorder()}: Faster version with no checks.
Does not retain labels or edge weights.
Edges must be listed in postorder.
May crash if improper input is specified.

}}
\examples{
tree <- BalancedTree(9)
plot(tree)
plot(DropTip(tree, c("t5", "t6")))

unrooted <- UnrootTree(tree)
plot(unrooted)
plot(DropTip(unrooted, 4:5))

summary(DropTip(as.Splits(tree), 4:5))
}
\seealso{
Other tree manipulation: 
\code{\link{AddTip}()},
\code{\link{CollapseNode}()},
\code{\link{ConsensusWithout}()},
\code{\link{EnforceOutgroup}()},
\code{\link{ImposeConstraint}()},
\code{\link{KeptPaths}()},
\code{\link{KeptVerts}()},
\code{\link{LeafLabelInterchange}()},
\code{\link{MakeTreeBinary}()},
\code{\link{Renumber}()},
\code{\link{RenumberTips}()},
\code{\link{RenumberTree}()},
\code{\link{RootTree}()},
\code{\link{SortTree}()},
\code{\link{Subtree}()},
\code{\link{TipTimedTree}()},
\code{\link{TrivialTree}}

Other split manipulation functions: 
\code{\link{Subsplit}()},
\code{\link{TrivialSplits}()}
}
\author{
\href{https://orcid.org/0000-0001-5660-1727}{Martin R. Smith}
(\href{mailto:martin.smith@durham.ac.uk}{martin.smith@durham.ac.uk})
}
\concept{split manipulation functions}
\concept{tree manipulation}

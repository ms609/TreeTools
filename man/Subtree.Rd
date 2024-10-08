% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/phylo.R
\name{Subtree}
\alias{Subtree}
\title{Extract a subtree}
\usage{
Subtree(tree, node)
}
\arguments{
\item{tree}{A tree of class \code{\link[ape:read.tree]{phylo}}, with internal
numbering in cladewise order (use \code{\link{Preorder}(tree)} or (slower)
\code{\link{Cladewise}(tree)}) .}

\item{node}{The number of the node at the base of the clade to be extracted.}
}
\value{
\code{Subtree()} returns a tree of class \code{phylo} that represents a
clade extracted from the original tree.
}
\description{
\code{Subtree()} safely extracts a clade from a phylogenetic tree.
}
\details{
Modified from the \pkg{ape} function \code{\link[ape]{extract.clade}}, which
sometimes behaves unpredictably.
Unlike extract.clade, this function supports the extraction of "clades"
that constitute a single tip.
}
\examples{
tree <- Preorder(BalancedTree(8))
plot(tree)
ape::nodelabels()
ape::nodelabels(13, 13, bg="yellow")

plot(Subtree(tree, 13))

}
\seealso{
Other tree manipulation: 
\code{\link{AddTip}()},
\code{\link{CollapseNode}()},
\code{\link{ConsensusWithout}()},
\code{\link{DropTip}()},
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
\code{\link{TipTimedTree}()},
\code{\link{TrivialTree}}
}
\author{
\href{https://orcid.org/0000-0001-5660-1727}{Martin R. Smith}
(\href{mailto:martin.smith@durham.ac.uk}{martin.smith@durham.ac.uk})
}
\concept{tree manipulation}

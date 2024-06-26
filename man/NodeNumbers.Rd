% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/tree_properties.R
\name{NodeNumbers}
\alias{NodeNumbers}
\title{Numeric index of each node in a tree
\code{NodeNumbers()} returns a sequence corresponding to the nodes in a tree}
\usage{
NodeNumbers(tree, tips = FALSE)
}
\arguments{
\item{tree}{A tree of class \code{\link[ape:read.tree]{phylo}}.}

\item{tips}{Logical specifying whether to also include the indices of leaves.}
}
\value{
\code{NodeNumbers()} returns an integer vector corresponding to the
indices of nodes within a tree.
}
\description{
Numeric index of each node in a tree
\code{NodeNumbers()} returns a sequence corresponding to the nodes in a tree
}
\seealso{
Other tree properties: 
\code{\link{ConsensusWithout}()},
\code{\link{MatchEdges}()},
\code{\link{NSplits}()},
\code{\link{NTip}()},
\code{\link{PathLengths}()},
\code{\link{SplitsInBinaryTree}()},
\code{\link{TipLabels}()},
\code{\link{TreeIsRooted}()}

Other tree navigation: 
\code{\link{AncestorEdge}()},
\code{\link{CladeSizes}()},
\code{\link{DescendantEdges}()},
\code{\link{EdgeAncestry}()},
\code{\link{EdgeDistances}()},
\code{\link{ListAncestors}()},
\code{\link{MRCA}()},
\code{\link{MatchEdges}()},
\code{\link{NDescendants}()},
\code{\link{NodeDepth}()},
\code{\link{NodeOrder}()},
\code{\link{RootNode}()}
}
\author{
\href{https://orcid.org/0000-0001-5660-1727}{Martin R. Smith}
(\href{mailto:martin.smith@durham.ac.uk}{martin.smith@durham.ac.uk})
}
\concept{tree navigation}
\concept{tree properties}

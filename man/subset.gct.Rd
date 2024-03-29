% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/ptm_utils.R
\name{subset.gct}
\alias{subset.gct}
\title{Subset a gct object using the provided row and column ids}
\usage{
subset.gct(g, rid = NULL, cid = NULL)
}
\arguments{
\item{g}{a gct object}

\item{rid}{a vector of character ids or integer indices for ROWS}

\item{cid}{a vector of character ids or integer indices for COLUMNS}
}
\description{
Subset a gct object using the provided row and column ids
}
\seealso{
Other GCT utilities: 
\code{\link{annotate.gct}()},
\code{\link{melt.gct}()},
\code{\link{merge.gct}()},
\code{\link{rank.gct}()}
}
\concept{GCT utilities}

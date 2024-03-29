% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/ptm_utils.R
\name{annotate.gct}
\alias{annotate.gct}
\title{Add annotations to a GCT object}
\usage{
annotate.gct(g, annot, dimension = "row", keyfield = "id")
}
\arguments{
\item{g}{a GCT object}

\item{annot}{a \code{\link{data.frame}} or path to text table of annotations}

\item{dimension}{either 'row' or 'column' indicating which dimension
of \code{g} to annotate}

\item{keyfield}{the character name of the column in \code{annot} that
matches the row or column identifiers in \code{g}}
}
\value{
a GCT object with annotations applied to the specified
  dimension
}
\description{
Given a GCT object and either a \code{\link{data.frame}} or
a path to an annotation table, apply the annotations to the
gct using the given \code{keyfield}.
}
\seealso{
Other GCT utilities: 
\code{\link{melt.gct}()},
\code{\link{merge.gct}()},
\code{\link{rank.gct}()},
\code{\link{subset.gct}()}
}
\concept{GCT utilities}

% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/ptm_io.R
\name{.append.dim}
\alias{.append.dim}
\title{Append matrix dimensions to filename}
\usage{
.append.dim(ofile, mat, extension = "gct")
}
\arguments{
\item{ofile}{the file name}

\item{mat}{the matrix}

\item{extension}{the file extension}
}
\value{
a character string of the filename with
  matrix dimensions appended
}
\description{
Append matrix dimensions to filename
}
\details{
This is a helper function that most users
  will not use directly
}
\seealso{
Other GCTX parsing functions: 
\code{\link{.fix.datatypes}()},
\code{\link{.parse.gctx}()},
\code{\link{.process_ids}()},
\code{\link{.read.gctx.ids}()},
\code{\link{.read.gctx.meta}()},
\code{\link{.write.gct}()},
\code{\link{.write.gctx}()},
\code{\link{.write.gctx.meta}()}
}
\concept{GCTX parsing functions}
\keyword{internal}

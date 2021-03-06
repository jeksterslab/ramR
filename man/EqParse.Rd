% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/EqParse.R
\name{EqParse}
\alias{EqParse}
\title{Equations Parser}
\usage{
EqParse(eq)
}
\arguments{
\item{eq}{Character string. Equations. See Details.}
}
\value{
Returns a data.frame with the following columns

\describe{
\item{lhs}{
is the variable on the \strong{left-hand side},
}
\item{rhs}{
is the variable on the \strong{right-hand side},
}
\item{op}{
is the \strong{operation} between \code{lhs} and \code{rhs},
}
\item{par.label}{
is the column of \strong{parameter label},
}
\item{par.start}{
is the column of \strong{starting values} for estimation
if \code{eq} has a fifth colulmn, and
}
\item{par.names}{
is the column of \strong{parameter label}
with \code{NAs} on fixed parameters.
}
}
}
\description{
Parse equations and return a parameter table (\code{par.table}).
The argument \code{eq} is a character string
that specifies the associations between the variables.
See \code{Syntax}, \code{Operations}, \code{par.label},
\code{par.start}, \verb{Line breaks}, and \code{Comments} below.
}
\section{Syntax}{

Each line should follow the syntax below

\verb{lhs <space> op <space> rhs <space> par.label <\\n> or <;>}

\describe{
\item{lhs}{
is the variable on the \strong{left-hand side},
}
\item{rhs}{
is the variable on the \strong{right-hand side},
}
\item{op}{
is the \strong{operation} between \code{lhs} and \code{rhs},
}
\item{par.label}{
is the column of \strong{parameter label},
}
\item{par.start}{
is the column of \strong{starting values} for estimation
(only for functions that require it), and
}
\item{\\n or ;}{
are \strong{line breaks}.
\strong{Each line should end with a line break.}
}
}
}

\section{Operations}{

The associations are defined by the following operations

\describe{
\item{by}{
\verb{left-hand side} measured \strong{by} \verb{right-hand side},
}
\item{on}{
\verb{left-hand side} regressed \strong{on} \verb{right-hand side},
}
\item{with}{
\verb{left-hand side} covarying \strong{with} \verb{right-hand side},
}
\item{on 1}{
\verb{left-hand side} regressed \strong{on 1} for mean structure.
}
}
}

\section{par.label}{

Each parameter should be labeled.
The \code{par.label} should be a number for fixed parameters
and a character string for free parameters.
Equality contraints can be imposed by using the same \code{par.label}.
}

\section{par.start}{

Numerical values as starting values for estimation.
Starting values for fixed parameters should be \code{NA}.
Starting values should be identical
for parameters constrained to be equal,
that is, parameters with the same \code{par.label}.
}

\section{Line breaks}{

The characters \verb{\\n} and \verb{;} can be used as line breaks.
\strong{Each line should end with a line break.}
}

\section{Comments}{

Comments can be written after a hash (\verb{#}) sign.
}

\examples{
eq <- "
  # lhs op   rhs par.label
    e   by   y   1
    y   on   x   beta
    e   with e   sigmae2
    x   with x   sigmax2
    y   on   1   alpha
    x   on   1   mux
"
EqParse(eq)

eq <- "
  # lhs op   rhs value
    e   by   y   1
    y   on   x   1
    e   with e   1
    x   with x   0.25
    y   on   1   0
    x   on   1   0.50
"
EqParse(eq)

eq <- "
  # lhs op   rhs par.label par.start
    e   by   y   1         NA
    y   on   x   beta      0.00
    e   with e   sigmae2   1.25
    x   with x   sigmax2   0.25
    y   on   1   alpha     0.00
    x   on   1   mux       0.50
"
EqParse(eq)

# \n as line breaks------------------------------------------------

eq <- "
  e by y 1 NA \n y on x beta 0.00
  e with e sigmae2 1.25 \n x with x sigmax2 0.25
  y on 1 alpha 0.00 \n x on 1 mux 0.50
"
EqParse(eq)

# ; as line breaks------------------------------------------------

eq <- "
  e by y 1 NA; y on x beta 0.00
  e with e sigmae2 1.25; x with x sigmax2 0.25
  y on 1 alpha 0.00; x on 1 mux 0.50
"
EqParse(eq)
}
\seealso{
Other eq functions: 
\code{\link{Eq2Expectations}()},
\code{\link{Eq2RAM}()}
}
\author{
Ivan Jacob Agaloos Pesigan
}
\concept{eq functions}
\keyword{eq}

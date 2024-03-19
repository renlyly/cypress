\name{cypress-methods}
\alias{cypress-methods}
\alias{getcypress}
\alias{setcypress}
\alias{show,cypress_out-method}
\alias{show,est_out-method}

\alias{setcypress,cypress_out-method}
\alias{setcypress,design_in-method}
\alias{setcypress,est_out-method}
\alias{getcypress,cypress_out-method}
\alias{getcypress,design_in-method}
\alias{getcypress,est_out-method}

\title{Get a slot from cypress output}
\usage{
getcypress(object, name)
setcypress(object, name, value)

}

\arguments{
\item{object}{object from cypress.}
\item{name}{name of the slot in cypress object.}
\item{value}{value of the slot in cypress object.}
}

\value{
Methods for cypress.
}

\description{
Accessor function for getting or replace slots.
Show methods for cypress object.
}


\examples{
data(quickPowerIBD)

getcypress(ibd_propPower, "ct_TDR_bio_smry")

}

% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/utils.R
\name{prepare_box_data}
\alias{prepare_box_data}
\title{Prepare Box Data}
\usage{
prepare_box_data(data, cat_a, cat_b, group, cat_a_order, cat_b_order)
}
\arguments{
\item{data}{A data frame containing the variables.}

\item{cat_a}{The name of the column representing category A.}

\item{cat_b}{The name of the column representing category B.}

\item{group}{The name of the column representing the grouping variable.}

\item{cat_a_order}{A vector specifying the order of category A.}

\item{cat_b_order}{A vector specifying the order of category B.}
}
\value{
A data frame with box boundaries for plotting.
}
\description{
Prepares data for plotting boxes by calculating box boundaries based on category positions.
}
\examples{
library(dplyr)
data <- data.frame(
  cat_a = rep(letters[1:3], each = 2),
  cat_b = rep(LETTERS[1:2], times = 3),
  group = rep(c("G1", "G2"), times = 3)
)
cat_a_order <- c("a", "b", "c")
cat_b_order <- c("A", "B")
prepare_box_data(data, "cat_a", "cat_b", "group", cat_a_order, cat_b_order)
}

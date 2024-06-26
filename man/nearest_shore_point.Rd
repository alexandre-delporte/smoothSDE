% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/utility.R
\name{nearest_shore_point}
\alias{nearest_shore_point}
\title{Compute nearest point on the land (shoreline)}
\usage{
nearest_shore_point(point, coastline)
}
\arguments{
\item{point}{the location from which wa want to find the nearest point on land}

\item{coastline}{sf object (list of polygons) defining the land}
}
\value{
matrix with on row and two columns that are the coordinates of the nearest point
}
\description{
Compute nearest point on the land (shoreline)
}

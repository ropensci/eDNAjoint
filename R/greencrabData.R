#' greencrabData
#' @srrstats {G5.1} data that is exported and is available for the the user can
#'   confirm tests and run examples
#'
#' European green crab data from eDNA and traditional (baited trap) surveys
#' collected from 20 locations in Washington state. Data can be used to
#' demonstrate functions in the R package eDNAjoint. Data is retrieved from
#' Keller et al. 2022.
#'
#' @title greencrabData
#' @docType data
#' @keywords data
#' @format A list with four matrices representing eDNA sampling data (qPCR.N and
#'   qPCR.K) and trap sampling data (count and count.type).
#' \describe{
#'   \item{qPCR.N}{Total number of eDNA qPCR replicates at each site (row) and
#'     eDNA sample replicate (column). Data includes 20 total sites and 5 eDNA
#'     sample replicates.}
#'   \item{qPCR.K}{Total number of positive eDNA qPCR detections at each site
#'     (row) and eDNA sample replicate (column). Data includes 20 total sites
#'     and 5 eDNA sample replicates.}
#'   \item{count}{Count of green crab individuals in trap samples at each site
#'     (row) and trap sample replicate (column). Data includes 20 total sites
#'     and a maximum of 420 trap replicates. NA indicates that fewer trap
#'     samples were collected than the maximum at a site.}
#'   \item{count.type}{Integer indicating the traditional gear type used at each
#'     site (row) and trap sample replicate (column). '1' refers to Fukui traps,
#'     and '2' refers to Minnow traps. Data includes 20 total sites and a
#'     maximum of 420 trap replicates. NA indicates that fewer trap samples
#'     were collected than the maximum at a site.}
#'}
#' @source \doi{10.6084/m9.figshare.15117102.v2}
#' @references Keller, A.G., Grason, E.W., McDonald, P.S., Ramon-Laca, A.,
#'   Kelly, R.P. (2022). Tracking an invasion front with environmental DNA.
#'   Ecological Applications. 32(4): e2561. https://doi.org/10.1002/eap.2561
"greencrabData"

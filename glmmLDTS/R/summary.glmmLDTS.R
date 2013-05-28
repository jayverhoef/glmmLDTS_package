#-------------------------------------------------------------------------------
#
#           summaryglmmLDTS
#
#-------------------------------------------------------------------------------

#' summary
#'
#' Tests for Tag Loss/Malfunction at end of telemetry data sequence.  Tests
#' are performed subject by subject, so a subject column must be included in
#' the data.frame
#'
#' @param datain a data frame of haul-out data 
#' @param SubjectColName column containing subject IDs 
#' @param HauloutColName column containing haul-out data
#'
#' @return a data.frame with the subject ID and the number of days where the
#'   subject has 100% haulout.  If the number of days is large, then it is
#'   likely that the tag fell of on the surface and kept transmitting dry times.
#'
#' @author Jay Ver Hoef
#' @export

summary.glmmLDTS <- function(x) {
  list(
    WARNINGS = x$WARNINGS,
    fixed.formula = x$fixed.formula,
    random.formula = x$random.formula,
    sample.size = x$sample.size,
    timecol = x$timecol,
    trialscol = x$trialscol,
    group.vec = x$group.vec,
    ridge.reg = x$ridge.reg,
    lambda = x$lambda,
    start.time = x$start.time,
    end.time = x$end.time,
    R.cov.parameters = x$R.cov.parameters,
    G.cov.parameters = x$G.cov.parameters,
    typeIII.hypoth = x$typeIII.hypoth,
    fixed.effects = x$fixed.effects,
    random.effects = x$random.effects,
    outer.iterations = x$outer.iterations,
    inner.interations = x$inner.iterations2)
}


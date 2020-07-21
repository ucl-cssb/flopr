### FUNCTIONS FOR APPENDING COLUMNS TO FLOW FRAME
### from updated version of flowCore not available to R < 4.0

#' Append data columns to a flowFrame
#'
#' Append data columns to a flowFrame
#'
#' It is used to add extra data columns to the existing flowFrame.  It handles
#' keywords and parameters properly to ensure the new flowFrame can be written
#' as a valid FCS through the function \code{write.FCS} .
#'
#' @param fr A \code{\link[flowCore:flowFrame-class]{flowFrame}}.
#' @param cols A numeric matrix containing the new data columns to be added.
#' Must has column names to be used as new channel names.
#'
#' @name ff_append_cols
#' @aliases ff_append_cols
#' @usage
#' ff_append_cols(fr, cols)
#' @return
#'
#' A \code{\linkS4class{flowFrame}}
#' @author Mike Jiang
#' @keywords IO
ff_append_cols <- function(fr, cols){
  new_pd <- cols_to_pd(fr, cols)
  pd <- flowCore::pData(flowCore::parameters(fr))
  pd <- rbind(pd, new_pd)
  #add to exprs
  fr@exprs <- cbind(flowCore::exprs(fr), cols)
  flowCore::pData(flowCore::parameters(fr)) <- pd

  update_kw_from_pd(fr, new_pd)
}

update_kw_from_pd <- function(fr, new_pd)
{
  new_pid <- rownames(new_pd)
  #take care of flowCore_$PnRmax
  trans <- flowCore::keyword(fr)[["transformation"]]
  if(!is.null(trans) && trans == "custom"){
    flowCore::keyword(fr)[paste0("flowCore_", new_pid, "Rmax")] <- new_pd[new_pid, "maxRange"]
    flowCore::keyword(fr)[paste0("flowCore_", new_pid, "Rmin")] <- new_pd[new_pid, "minRange"]
  }
  fr
}
#' generate new pData of flowFrame based on the new cols added
#' @param fr A \code{\link[flowCore:flowFrame-class]{flowFrame}}.
#' @param cols A numeric matrix containing the new data columns to be added.
#' Must has column names to be used as new channel names.
#' @noRd
cols_to_pd <- function(fr, cols){
  checkmate::checkClass(cols, "matrix")
  ncol <- ncol(cols)
  cn <- colnames(cols)
  if(length(cn) != ncol)
    stop("All columns in 'cols' must have colnames!")
  #add to pdata
  pd <- flowCore::pData(flowCore::parameters(fr))
  ncol_old <- ncol(fr)
  new_pid <- max(as.integer(gsub("\\$P", "", rownames(pd)))) + 1
  new_pid <- seq(new_pid, length.out = ncol)
  new_pid <- paste0("$P", new_pid)

  new_pd <- do.call(rbind, lapply(cn, function(i){
    vec <- cols[,i]
    rg <- range(vec)
    data.frame(name = i, desc = NA, range = diff(rg) + 1, minRange = rg[1], maxRange = rg[2])
  }))
  rownames(new_pd) <- new_pid
  new_pd
}

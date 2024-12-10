
#' Create block diagonal matrix (safe version)
#' 
#' This version of bdiag checks whether the matrices passed as
#' arguments are NULL. This avoids errors that would arise if
#' using bdiag directly.
#' 
#' @param ... Matrix or list of matrices (only the first argument is used)
#' 
#' @return Block diagonal matrix
#' 
#' @importFrom Matrix bdiag
#' @importFrom stringr str_extract
bdiag_check <- function(...) {
    # Only use first argument (matrix or list of matrices)
    args <- list(...)[[1]]
    
    # Check which matrices are non-NULL
    # (Needs two conditions to keep non-empty vectors which don't have
    # dimensions, and also keep empty matrices with set dimensions)
    check <- sapply(args, function(arg) {
        !is.null(dim(arg)) | length(arg) > 0
    })
    
    if(length(check) == 0)
        return(NULL)
    else
        return(bdiag(args[check]))
}

#' Grid of covariates
#' 
#' @param var Name of variable
#' @param data Data frame containing the covariates
#' @param covs Optional data frame with a single row and one column
#' for each covariate, giving the values that should be used. If this is
#' not specified, the mean value is used for numeric variables, and the
#' first level for factor variables.
#' @param formulas List of formulas used in the model
#' 
#' @return Data frame of covariates, with 'var' defined over a grid,
#' and other covariates fixed to their mean (numeric) or first level
#' (factor).
cov_grid <- function(var, data, covs = NULL, formulas) {
    # Get covariate names
    var_names <- unique(rapply(formulas, all.vars))
    
    # pi might appear in the formulas (e.g. used in periodic terms),
    # in which case it needs to be added to the data frame
    if(any(var_names == "pi")) {
        data$pi <- pi
    }
    
    # Get data frame of covariates
    all_vars <- data[, var_names, drop = FALSE]
    
    # Grid of covariate
    if(is.factor(all_vars[, var])) {
        n_grid <- length(unique(all_vars[, var]))
        grid <- unique(all_vars[, var])
    } else {
        n_grid <- 1e3
        grid <- seq(min(all_vars[, var]), max(all_vars[, var]), length = n_grid)
    }
    
    # New data frame for covariate grid
    new_data <- matrix(NA, nrow = n_grid, ncol = ncol(all_vars))
    colnames(new_data) <- colnames(all_vars)
    new_data <- as.data.frame(new_data)
    new_data[, var] <- grid
    
    # Select value for other covariates
    covs_list <- lapply(seq_along(all_vars), function(i) {
        if(!is.null(covs[[var_names[i]]])) {
            # Set to user-provided value if possible
            return(covs[[var_names[i]]])
        } else {
            # No user-provided value
            col <- all_vars[, i]
            if(is.numeric(col)) {
                # If numeric, use mean value
                return(mean(col, na.rm = TRUE)) 
            } else {
                # If factor, use first factor level
                return(unique(col)[1])
            }
        }
    })
    covs <- as.data.frame(covs_list)
    colnames(covs) <- colnames(all_vars)
    
    # Fill columns for other covariates
    for(var_name in colnames(new_data)) {
        if(var_name != var)
            new_data[, var_name] <- covs[, var_name]
    }
    
    return(new_data)
}

#' logLik function for SDE objects
#' 
#' This function makes it possible to call generic R methods such
#' as AIC and BIC on SDE objects. It is based on the number of
#' degrees of freedom of the *conditional* AIC (rather than
#' marginal AIC), i.e., including degrees of freedom from the
#' smooth/random effect components of the model.
#' 
#' @param object SDE model object
#' @param ... For compatibility with S3 method
#' 
#' @return Maximum log-likelihood value for the model, with attributes
#' \code{df} (degrees of freedom) and \code{nobs} (number of observations)
#' 
#' @export
logLik.SDE <- function(object, ...) {
    par_all <- c(object$tmb_rep()$par.fixed, 
                 object$tmb_rep()$par.random)
    val <- - object$tmb_obj_joint()$fn(par_all)
    attributes(val)$df <- object$edf_conditional()
    attributes(val)$nobs <- nrow(object$data())
    class(val) <- "logLik"
    return(val)
}

#' Indices of coefficients for given model term
#' 
#' @param names_fe Names of fixed effect coefficients
#' @param names_re Names of random effect coefficients
#' @param term Name of term as character string, e.g. "time", 
#' or "s(time)"
#' 
#' @return List with elements \code{fe} and \code{re}, which are vectors
#' of the indices in coeff_fe and coeff_re (respectively) corresponding
#' to the term.
#' 
#' @export
term_indices <- function(names_fe, names_re, term) {
    # Find indices of coeff_fe and coeff_re that we want to keep
    # (this is simply done by finding those that have 'term' in their names)
    wh_keep_fe <- grep(term, names_fe, fixed = TRUE)
    wh_keep_re <- grep(term, names_re, fixed = TRUE)
    
    return(list(fe = wh_keep_fe, re = wh_keep_re))
}

#' Get covariance matrix from precision matrix
#' 
#' The covariance matrix is the inverse of the precision matrix. By default,
#' the function \code{solve} is used for inversion. If it fails (e.g.,
#' singular system), then \code{MASS::ginv} is used instead, and returns the
#' Moore-Penrose generalised inverse of the precision matrix.
#' 
#' @param prec_mat Precision matrix (either of 'matrix' type
#' or sparse matrix on which as.matrix can be used)
#' 
#' @return Precision matrix
#' 
#' @importFrom MASS ginv
#' @export
prec_to_cov <- function(prec_mat) {
    cov_mat <- try(as.matrix(solve(prec_mat)), silent = TRUE)
    if(inherits(cov_mat, "try-error")) {
        message <- attr(cov_mat, 'condition')$message
        cov_mat <- MASS::ginv(as.matrix(prec_mat))
        warning(paste0("Inversion of precision matrix using 'solve' failed: ", 
                       message, ". Using 'MASS::ginv' instead (uncertainty ",
                       "estimates may be unreliable)."))
    }
    colnames(cov_mat) <- colnames(prec_mat)
    rownames(cov_mat) <- colnames(prec_mat)
    return(cov_mat)
}

#' Covariance matrix of CTCRW transition density
#' 
#' The derivation of these equations is for example included in
#' Michelot et al. (2019), Stochastic models of animal movement and
#' habitat selection (Sections 6.2.2.3-6.2.2.4).
#' 
#' @param beta Parameter beta of CTCRW model (mean reversion)
#' @param sigma Parameter sigma of CTCRW model (diffusion)
#' @param dt Time interval
#' 
#' @return Covariance matrix of position/velocity joint 
#' process
#' 
#' @export
CTCRW_cov <- function(beta, sigma, dt) {
    Q <- matrix(0, 2, 2)
    Q[1,1] <- sigma^2/(2*beta) * (1-exp(-2*beta*dt))
    Q[2,2] <- (sigma/beta)^2 * (dt + (1-exp(-2*beta*dt))/(2*beta) -
                                    2*(1-exp(-beta*dt))/beta)
    Q[1,2] <- sigma^2/(2*beta^2) * (1 - 2*exp(-beta*dt) + exp(-2*beta*dt))
    Q[2,1] <- Q[1,2]
    return(Q)   
}


#' Covariance matrix of RACVM transition density
#' 
#' 
#' @param beta Parameter beta of RACVM model (mean reversion)
#' @param sigma Parameter sigma of RACVM model (diffusion)
#' @param omega Parameter omega of RACVM (angular speed)
#' @param dt Time interval
#' 
#' @return Covariance matrix of position/velocity joint 
#' process
#' 
#' @export
RACVM_cov <- function(beta, sigma,omega,dt) {
  
  #relevant matrices
  C=beta^2+omega^2
  A=matrix(c(beta,-omega,omega,beta),nrow=2,byrow=TRUE)
  invA=1/C*matrix(c(beta,omega,-omega,beta),nrow=2,byrow=TRUE)
  R=matrix(c(cos(omega*dt),sin(omega*dt),-sin(omega*dt),cos(omega*dt)),byrow=TRUE,nrow=2)
  expAdt=exp(-beta*dt)*R
  
  # Covariance of next state vector
  var_xi=sigma^2/C*(dt+(omega^2-3*beta^2)/(2*beta*C)-exp(-2*dt*beta)/(2*beta)+
                      2*exp(-dt*beta)*(beta*cos(omega*dt)-omega*sin(omega*dt))/C)
  var_zeta=sigma^2/(2*beta)*(1-exp(-2*dt*beta))
  cov1=sigma^2/(2*C)*(1+exp(-2*dt*beta)-2*exp(-dt*beta)*cos(omega*dt))
  cov2=sigma^2/C*(exp(-dt*beta)*sin(omega*dt)-omega/(2*beta)*(1-exp(-2*dt*beta)))
  Q=matrix(c(var_xi,0,cov1,cov2,0,var_xi,-cov2,cov1,cov1,-cov2,var_zeta,0,cov2,cov1,0,var_zeta),nrow=4,byrow=TRUE)
  
  return(Q)   
}


#' Link matrix of RACVM transition density
#' 
#' 
#' @param beta Parameter beta of RACVM model (mean reversion)
#' @param sigma Parameter sigma of RACVM model (diffusion)
#' @param omega Parameter omega of RACVM (angular speed)
#' @param dt Time interval
#' 
#' @return Link matrix of position/velocity joint 
#' process
#' 
#' @export
RACVM_link <- function(beta,omega,dt) {
  
  #relevant matrices
  C=beta^2+omega^2
  A=matrix(c(beta,-omega,omega,beta),nrow=2,byrow=TRUE)
  invA=1/C*matrix(c(beta,omega,-omega,beta),nrow=2,byrow=TRUE)
  R=matrix(c(cos(omega*dt),sin(omega*dt),-sin(omega*dt),cos(omega*dt)),byrow=TRUE,nrow=2)
  expAdt=exp(-beta*dt)*R
  
  # link matrix
  L <- matrix(0, nrow = 4, ncol = 4)
  L[1:2, 1:2] <- diag(2)
  L[1:2, 3:4] <- invA%*%(diag(2)-expAdt)
  L[3:4, 1:2] <- matrix(c(0,0,0,0),nrow=2)
  L[3:4, 3:4] <- expAdt
  return(L)   
}

#' Drift matrix of RACVM transition density
#' 
#' 
#' @param beta Parameter beta of RACVM model (mean reversion)
#' @param sigma Parameter sigma of RACVM model (diffusion)
#' @param omega Parameter omega of RACVM (angular speed)
#' @param dt Time interval
#' 
#' @return Link matrix for drift term of position/velocity joint 
#' process
#' 
#' @export
RACVM_drift <- function(beta,omega,dt) {
  
  #relevant matrices
  C=beta^2+omega^2
  A=matrix(c(beta,-omega,omega,beta),nrow=2,byrow=TRUE)
  invA=1/C*matrix(c(beta,omega,-omega,beta),nrow=2,byrow=TRUE)
  R=matrix(c(cos(omega*dt),sin(omega*dt),-sin(omega*dt),cos(omega*dt)),byrow=TRUE,nrow=2)
  expAdt=exp(-beta*dt)*R
  
  #matrix for drift term
  B<- matrix(0, nrow = 4, ncol = 2)
  B[1:2, 1:2] <- dt*diag(2)-invA%*%(diag(2)-expAdt)
  B[3:4,1:2] <- diag(2)-expAdt
  
  return(B)   
}


#' Transforms matrix to dgTMatrix
#' 
#' @param x Matrix or vector. If this is a vector, it is formatted into
#' a single-column matrix.
#' 
#' @return Sparse matrix of class dgTMatrix
as_sparse <- function(x) {
    if(length(dim(x)) < 2) {
        x <- matrix(x, ncol = 1)
    }
    # # This is the syntax recommended by Matrix > 1.5.0, but doesn't seem
    # # to be compatible with earlier versions of Matrix.
    # mat <- as(as(as(x, "dMatrix"), "generalMatrix"), "TsparseMatrix")
    mat <- suppressMessages(as(x, "dgTMatrix"))
    return(mat)
}


#' Split string with parenthesis
#' @param string a string to split based on character sep
#' @param sep separation parameter
#' 
#' @return list of string split based on sep, but ignoring sep that are within parenthesis

split_term <- function(string,sep) {
    parts <- character(0)
    current_part <- ""
    parenthesis_count <- 0
    
    for (char in unlist(strsplit(string, ""))) {
        if (char == sep && parenthesis_count == 0) {
            parts <- c(parts, trimws(current_part))
            current_part <- ""
        } else {
            current_part <- paste(current_part, char, sep = "")
            if (char == "(") {
                parenthesis_count <- parenthesis_count + 1
            } else if (char == ")") {
                parenthesis_count <- parenthesis_count - 1
            }
        }
    }
    
    if (nchar(current_part) > 0) {
        parts <- c(parts, trimws(current_part))
    }
    
    return(parts)
}

#' Get covariable names from a formula
#' @param formula a formula (usually for one parameter of the SDE) 
#'
#'@return Vector of names of the covariables in the formula
get_variables=function(formula) {
  
  #list of strings for each term in the formula
  terms=colnames(attr(terms(formula),which="factors"))
  
  #initialize vector of covariates
  covariates=c()
  
  #loop over the formula terms
  for (term in  terms) {
      
    # Metacharacters to get only what is after the first parenthesis
    pattern="\\(.*"
    
    # Find the substring within parentheses (ex : "(X,k=5,bs='cs')" )
    match=str_extract(term,pattern)
    
    #if no parentheses, assume it is a linear term
    if (is.na(match)) {
        covariates=c(covariates,term)
    }
    else {
        # Remove the parenthesis
        match=substring(match,2,nchar(match)-1)
      
        if (length(match) > 0) {
            new_covs <- split_term(match,sep=",")
            # Filter out terms containing "="
            new_covs <- new_covs[!grepl("=", new_covs)]
            covariates=c(covariates,new_covs)
        }
    }
  }

  return (unique(covariates))
}



#' Compute signed angle in [-pi,pi] that rotates first vector into second vector
#' 
#' @param u matrice of bivariate vectors with 2 columns
#' @param v matrice of bivariate vectors with 2 columns
#' @return vector of signed angles in [-pi,pi] that rotates first vector into second vector
#' # https://math.stackexchange.com/questions/529555/signed-angle-between-2-vectors
#' 
#' @export
signed_angle <- function(u, v) {
  # u, v: matrices of bivariate vectors with 2 columns
  u <- matrix(u, ncol = 2)
  v <- matrix(v, ncol = 2)
  if (nrow(u) != nrow(v)) stop("u and v must have the same number of 
                                  rows")
  result <- as.numeric(atan2(v[,2], v[,1]) - atan2(u[,2], u[,1]))
  ind1 <- which(result > pi)
  ind2 <- which(result <= -pi)
  result[ind1] <- result[ind1] - 2*pi
  result[ind2] <- result[ind2] + 2*pi
  return(result) 
} 

#' Check whether a point is on land or not
#' 
#' @param point the location we want to tell if it is on land or not
#' @param land sf object (list of polygons) defining the land
#' @return boolean TRUE if it is on land, FALSE otherwise
#' 
#' @export

is_in_land=function(point,land) {
    return (length(st_intersects(point,land$geometry)[[1]])>0)
}


#' Compute nearest point on the shoreline
#' 
#' @param point the location from which we want to find the nearest point on land
#' @param land sf object (list of polygons) defining the land
#' @return matrix with on row and two columns that are the coordinates of the nearest point in land
#' 
#' @export

nearest_boundary_point <- function(point, land) {
    # Initialize variables to store the minimum distance and nearest boundary point
    min_distance <- Inf
    nearest_point <- NULL
    
    # Iterate over each polygon in land
    for (j in 1:length(land$geometry)) {
        
        # Get the boundary of the current polygon
        boundary <- st_boundary(land$geometry[[j]])
        
        # Calculate distance and nearest point to the boundary
        distance <- st_distance(point, boundary)
        candidate <- st_nearest_points(point, boundary)
        
        # Ensure that the candidate is not empty
        if (!st_is_empty(candidate)) {
            coordinates <- st_coordinates(candidate)[2, c("X", "Y")]
            
            # Update the nearest point if this boundary point is closer
            if (distance < min_distance) {
                min_distance <- distance
                nearest_point <- coordinates
            }
        }
    }
    
    return(as.matrix(nearest_point))
}

#' Check if two columns of a dataframe are orthogonal
#' @param data dataframe
#' @param C1 column name
#' @param C2 column name
#' 
#' @return boolean
#' @export
are_orthogonals=function(data,C1,C2){
    inner_product=data[,C1]*data[,C2]
    return (all(inner_product==0))
}



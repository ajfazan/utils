least_squares_network <- function( obs, ref, alpha = 1.0, sigma_priori = 1.0 )
{
  # print( is.data.frame( obs) )

  stations <- unique( sort( c( as.vector( obs$source ), as.vector( obs$target ) ) ) )

  u <- 3 * length( stations )

  n <- 3 * ( nrow ( obs ) + nrow( ref ) )

  lb <- matrix( rep( 0.0, n ), nrow = n, ncol = 1 )

  a <- matrix( rep( 0.0, n * u ), nrow = n, ncol = u )

  cov_obs <- matrix( rep( 0.0, n * n ), nrow = n, ncol = n )

  j <- 1 : 3
  # s <- 0.0 # debug

  for( i in 1 : nrow( obs ) ) {

    lb[j[1],1] <- obs[i,"DX"]
    lb[j[2],1] <- obs[i,"DY"]
    lb[j[3],1] <- obs[i,"DZ"]

    cov_obs[j[1],j[1]] <- ( obs[i,"sigmaDX"] ^ 2 )
    cov_obs[j[2],j[2]] <- ( obs[i,"sigmaDY"] ^ 2 )
    cov_obs[j[3],j[3]] <- ( obs[i,"sigmaDZ"] ^ 2 )

    k <- 3 * which( stations == obs[i,"source"] )
    # s <- s + length( k ) # debug

    a[j[1],k-2] <- -1.0
    a[j[2],k-1] <- -1.0
    a[j[3],k  ] <- -1.0

    k <- 3 * which( stations == obs[i,"target"] )
    # s <- s + length( k ) # debug

    a[j[1],k-2] <- 1.0
    a[j[2],k-1] <- 1.0
    a[j[3],k  ] <- 1.0

    j <- j + 3
  }

  for( i in 1 : nrow( ref ) ) {

    lb[j[1],1] <- ref[i,"X"]
    lb[j[2],1] <- ref[i,"Y"]
    lb[j[3],1] <- ref[i,"Z"]

    cov_obs[j[1],j[1]] <- ( ( alpha * ref[i,"sigmaX"] ) ^ 2 )
    cov_obs[j[2],j[2]] <- ( ( alpha * ref[i,"sigmaY"] ) ^ 2 )
    cov_obs[j[3],j[3]] <- ( ( alpha * ref[i,"sigmaZ"] ) ^ 2 )

    k <- 3 * which( stations == ref[i,"station"] )
    # s <- s + length( k ) # debug

    a[j[1],k-2] <- 1.0
    a[j[2],k-1] <- 1.0
    a[j[3],k  ] <- 1.0

    j <- j + 3
  }

  # print( s == ( 2 * nrow( obs ) + nrow( ref ) ) ) # debug

  w <- sigma_priori * solve( cov_obs )

  atw <- t( a ) %*% w

  inm <- solve( atw %*% a )

  x <- inm %*% ( atw %*% lb )

  v <- a %*% x - lb

  delta <- n - u ## degrees of freedom

  sigma_post <- as.double( t( v ) %*% w %*% v )

  critical_chisq <- qchisq( 0.95, df = delta )

  if( sigma_post > critical_chisq ) {

    print( "Chi-square test failed" )
  }

  sigma_post <- sigma_post / delta

  cov_x <- sigma_post * inm

  targets <- unique( sort( obs$target ) )

  n <- length( targets )

  zeros <- rep( 0.0, n )

  obs <- data.frame( targets, zeros, zeros, zeros, zeros, zeros, zeros )
  colnames( obs ) <- c( "station", "X", "Y", "Z", "sigmaX", "sigmaY", "sigmaZ" )

  # d1 <- NULL
  # d2 <- NULL

  for( i in 1 : n ) {

    k <- 3 * which( stations == obs[i,"station"] ) - c( 2, 1, 0 )

    # d1 <- c( d1, k )
    # d2 <- c( d2, obs[i,"station"] )

    obs[i,"X"] <- x[k[1],1]
    obs[i,"Y"] <- x[k[2],1]
    obs[i,"Z"] <- x[k[3],1]

    obs[i,"sigmaX"] <- sqrt( cov_x[k[1],k[1]] )
    obs[i,"sigmaY"] <- sqrt( cov_x[k[2],k[2]] )
    obs[i,"sigmaZ"] <- sqrt( cov_x[k[3],k[3]] )
  }

  for( i in 1 : nrow( ref ) ) {

    n <- n + 1

    marker <- ref[i,"station"]
    obs[n,"station"] <- marker

    k <- 3 * which( stations == marker ) - c( 2, 1, 0 )

    # d1 <- c( d1, k )
    # d2 <- c( d2, marker )

    obs[n,"X"] <- x[k[1],1]
    obs[n,"Y"] <- x[k[2],1]
    obs[n,"Z"] <- x[k[3],1]

    obs[n,"sigmaX"] <- sqrt( cov_x[k[1],k[1]] )
    obs[n,"sigmaY"] <- sqrt( cov_x[k[2],k[2]] )
    obs[n,"sigmaZ"] <- sqrt( cov_x[k[3],k[3]] )
  }

  # print( sum( sort( d1 ) == seq.int( 1, nrow( x ) ) ) == nrow( x ) )
  # print( sum( sort( d2 ) == stations ) == length( stations ) )

  result <- list()
  result$stations <- obs
  result$sigma_post <- sigma_post
  result$critical_chisq <- critical_chisq
  result$delta <- delta
  # result$v <- v
  # result$w <- w

  return ( result )
}

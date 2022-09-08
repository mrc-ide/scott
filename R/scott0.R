#' Simulate a coalescent transmission tree 
#'
#' This function will return a phylogeny in the `ape` package format with additional attributes showing who infected whom. 
#' 
#' @param tfgy output of demo/epi simulation
#' @param s vector of sample times 
#' @param X matrix of sample states (probabilities)
#' @return A phylogenetic tree which contains a $waifw attribute. This is a data frame containing transmission pairs and time of transmission. 
#' @export 
scott <- function(
  tfgy
  , s
  , X 
  , min_dt = NULL 
  , max_dt = NULL 
  , event_tol = 1 
)
{
	times <- tfgy[[1]]
	Fs <- tfgy[[2]]
	Gs <- tfgy[[3]]
	Ys <- tfgy[[4]]
	fint <- tail( times, 1) #? 
	if ( is.null( min_dt ))
		min_dt <- abs( diff(times)[1] /100  )
	if ( is.null( max_dt ))
		max_dt <- median( abs( diff(times) ) )
	
	stopifnot( nrow(X) == length( s ))
	stopifnot( min_dt >= 0 )
	stopifnot( event_tol >= 0 )
	
	m <- nrow( Fs[[1]] )
	n <- length( s ) 
	nnode <- n -1 
	N <- n + nnode 
	
	states <- rep( NA, N) 
	states[1:n] <- apply( X, 1, function(x){
		sample.int( m, size =1 , prob = x )
	})
	extant <- rep( FALSE, N )
	sampled <- rep( FALSE, n )
	
	
	edge <- matrix( integer(0), nrow = N-1, ncol = 2 )
	edge.length <- integer( N -1 )
	tip.label <- paste0( 'tip_', 1:n)
	
	nodeHeights <- rep(Inf, N )
	nodesAdded <- 0 
	edgesAdded <- 0 
	
	pids <- 1:N
	nwaifw <- 0 
	waifw <- list() 
	
	ncos <- 0 
	tstart <- max( s) 
	xt = tstart 
	dxt <- -0 
	while(  (ncos < nnode)  &  (xt > fint)  ){
		# update xt  
		xt <- xt  + dxt  # dxt negative 
		xi <- which.min( xt <= times ) #??should get closest to xt 
		if ( length( xi ) == 0 )
			xi <- 1 
		h <- tstart - xt #?? tstart - xt; xt < tstart ; h > 0 
		# which are extant 
		# incorporate new samples 
		newsamps <- which( (!sampled) & (xt<=s) )
		sampled[newsamps] <- TRUE
		extant[newsamps] <- TRUE
		nodeHeights[newsamps] <- h 
		we = which( extant ) 
		A <- length( we ) 
		# protect against n nextant == 1
		# update pairwise co rates 
		Y <- Ys[[xi]]
		rF <- t(Fs[[xi]] / Y ) / Y 
		rF[ is.na( rF ) ] <- 0 
		extantgrid = as.matrix(  expand.grid( we , we ) ) 
		extantgrid <- extantgrid[ extantgrid[,1] != extantgrid[,2] , ]
		if ( A  >  1)
		{
			pwcorates <- rF[  matrix( states[ as.matrix( extantgrid )  ] , ncol=2 )  ] 
			corate <- sum( pwcorates )
			ncos <- rpois(1, corate * abs(dxt) )
		} else{ 
			ncos <- 0 
		}
		# update indiv mig rates 
		rG <- t(Gs[[xi]]) / Y 
		diag( rG ) <- 0 
		rG[ is.na( rG )] <- 0
		migratesmat <- rG[ states[ we ] , ]
		if ( A > 1 ){
			migrates <- rowSums( migratesmat )
		} else if ( A == 1 ){
			migrates <- sum( migratesmat )
		}
		migrate <- sum( migrates ) 
		nmigs <- rpois( 1, migrate * abs(dxt) ) 
		# update transm mig rates 
		rFmig <- t(Fs[[xi]]) / Y  # 
		rFmig[ is.na( rFmig ) ] <- 0 
		migratesmatF <- rFmig[ states[ we ] , ]
		if ( A > 1 ){
			migratesF <- rowSums( migratesmatF )
		} else if ( A == 1 ){
			migratesF <- sum( migratesmatF )
		}
		migrateF <- sum( migratesF ) 
		nmigsF <- rpois( 1, migrateF * abs(dxt) ) 
		# update time step 
		# update dxt so that few events per time step 
		dxt  <- -( min(max_dt, max( min_dt, event_tol / (migrateF + migrate + corate) ) ) )
		stopifnot( dxt < 0 )
		
		# do co's 
		if ( ncos >0 ) for ( i in 1:ncos )
		{
			j <- sample.int( length( pwcorates ), prob = pwcorates, size =1  )
			uv = extantgrid[j, ] 
			u <- uv[1] 
			v <- uv[2] 
			if ( all( extant[uv] ) )
			{ # needed since  a node can only co once per time step 
				nodesAdded <- nodesAdded + 1
				a <- nodesAdded + n 
				nodeHeights[a] <- h 
				edgesAdded <- edgesAdded + 1 
				edge[ edgesAdded, 1 ] <- a 
				edge[ edgesAdded, 2 ] <- u
				edge.length[ edgesAdded ] <- h - nodeHeights[u]
				edgesAdded <- edgesAdded + 1 
				edge[ edgesAdded, 1 ] <- a 
				edge[ edgesAdded, 2 ] <- v
				edge.length[ edgesAdded ] <- h - nodeHeights[v]
				extant[u] <- FALSE 
				extant[v] <- FALSE 
				extant[a] <- TRUE 
				# fill waifw , v is the transmitter 
				nwaifw <- nwaifw + 1 
				waifw[[ nwaifw ]] <- c(  pids[v], pids[u], xt  ) 
				pids[a] <- pids[v]
				states[a] <- states[v] 
			}
		}
		# do mig's 
		#we1 <- which( extant ) 
		if ( nmigs > 0 ) for ( i in 1:nmigs )
		{
			j <- sample.int( length( migrates ), prob = migrates, size = 1 )
			u <- we[ j ] 
			if ( extant[u] )# maybe not extant if co happened 
			{
				k <- sample.int( ncol(rG) , prob = rG[ states[u],] , size = 1) 
				states[u] <- k
			}
		}
		# transmission migs 
		if ( nmigsF > 0 )for ( i in 1:nmigsF )
		{
			j <- sample.int( length( migratesF ), prob = migratesF, size=1 )
			u <- we[ j ] 
			if ( extant[u] ) # maybe not extant if co happened 
			{
				# note this doesn't necessarily change state, but will update waifw
				k <- sample.int( ncol(rFmig) , prob = rFmig[ states[u] ,] , size = 1) 
				states[u] <- k # may be same as previous state
				# fill waifw 
				nwaifw <- nwaifw + 1 
				waifw[[ nwaifw]] <- c(  N+nwaifw, pids[u], xt  ) 
				pids[u] <- N + nwaifw 
			}
		}
	}
	if ( nodesAdded < nnode ) {
		# add polytomy root
		we = which( extant ) 
		nodesAdded <- nodesAdded + 1
		a <- n + nodesAdded
		#for ( i in 1:(nnode - ncos) ){
		for ( u in we ){
			edgesAdded <- edgesAdded + 1 
			edge[ edgesAdded, 1 ] <- a 
			edge[ edgesAdded, 2 ] <- u
			edge.length[ edgesAdded] <- h - nodeHeights[u] 
		}
		edge <- na.omit( edge ) 
		edge.length <- edge.length[ 1:nrow(edge) ]
	}
	
	waifw <- do.call( rbind, waifw )
	waifw <- as.data.frame( waifw )
	colnames( waifw ) <- c( 'donor', 'recipient', 'time' )
	tree = structure( list( edge = edge, edge.length = edge.length, tip.label = tip.label, Nnode = nodesAdded  ) , class = 'phylo')
	tree <- read.tree( text = write.tree( tree )) 
	tree = ladderize( tree )
	tree$waifw2 <- waifw 
	tree$waifw <- waifw[ waifw[,1] <= n & waifw[,2] <= n , ]
	tree
}



rm(list = ls())

##### PACKAGES #####

require(adehabitatLT)
require(adehabitatHR)
require(sp)
require(rgeos)
require(spatstat)
require(maptools)
require(raster)
require(CircStats)

##### SIMULATION SETUP #####

# Arena
ENR <- readShapePoly("Northern_Range_Bound.shp")
arena <- ENR
arena@data$AREA <- arena@data$AREA/1000000
arena@data$PERIMETER <- arena@data$PERIMETER/1000
arena@polygons[[1]]@Polygons[[1]]@labpt <- arena@polygons[[1]]@Polygons[[1]]@labpt/1000
arena@polygons[[1]]@Polygons[[1]]@area <- arena@polygons[[1]]@Polygons[[1]]@area/1000000
arena@polygons[[1]]@Polygons[[1]]@coords <- arena@polygons[[1]]@Polygons[[1]]@coords/1000
arena@polygons[[1]]@labpt <- arena@polygons[[1]]@labpt/1000
arena@polygons[[1]]@area <- arena@polygons[[1]]@area/1000000
arena@bbox <- arena@bbox/1000

# Number of steps for all trajectories
nsteps = 2160

# Number of model iterations
reps = 100

# Number of elk and wolves
nprey = 100
npred = 4

# Density of cameras (grid cell area)
n.cameras <- c(50,100,200,400)

# Wolf starting locations
wolf.starts <- spsample(arena,n=4,type="regular",offset=c(0.5,0.5), iter=10)

##### RANDOM MODEL #####

random <- function(
  abPred=npred,                       #wolf abundance 
  abPrey=nprey,                       #elk abundance
  n=nsteps,                           #number of movement steps per trajectory
  C2F=0.134,                          #commuting-->foraging transition probability for elk
  F2C=0.076,                          #foraging-->commuting transition probability for elk
  Fk=0.789,                           #shape parameter for elk foraging speed
  Flam=0.135,                         #scale parameter for elk foraging speed
  Fmu=0.057,                          #mean for elk foraging turning angle
  Fcon=0.414,                         #concentration for elk foraging turning angle
  Ck=0.809,                           #shape parameter for elk commuting speed
  Clam=0.343,                         #scale parameter for elk commuting speed
  Cmu=0.058,                          #mean for elk commuting turning angle
  Ccon=0.057,                         #concentration for elk commuting turning angle
  Sk=0.857,                           #shape parameter for wolf search speed
  Slam=0.856,                         #scale parameter for wolf search speed
  Smu=0.012,                          #mean for wolf search turning angle
  Scon=0.258                          #concentration for wolf search turning angle
){
  
  preytrajs <- list()
  predtrajs <- list()
  
  # ELK MOVEMENT
  
  for (p in 1:abPrey){
    
    # random starting location and direction
    x <- as.numeric(spsample(arena, n = 1, type = "random", iter=20)@coords)
    theta.x <- runif(1,-pi,pi)
    x.t <- x
    
    # random starting behaviour
    behaviours <- c("Foraging", "Commuting")
    behav <- sample(behaviours,1)
    behav2 <- behav
    
    for (i in 1:n){
      
      # change of behaviour based on transition probabilities
      if (behav2 == "Foraging"){
        behav.change <- rbinom(1,1,F2C)
        if (behav.change == 1){ 
          behav <- "Commuting"
        }
      }
      
      if (behav2 == "Commuting") {
        behav.change <- rbinom(1,1,C2F)
        if (behav.change == 1){ 
          behav <- "Foraging"
        }
      }
      
      # Movement step
      if (behav == "Foraging"){
        
        ta <- rwrpcauchy(n=1,location=Fmu,rho=Fcon)
        theta.x <- theta.x + ta
        d.x <- rweibull(1, shape=Fk, scale=Flam)/2.5
        x. <- x + c(d.x*cos(theta.x),d.x*sin(theta.x))
        
      }
      
      if (behav == "Commuting"){
        
        ta <- rwrpcauchy(1,location=Cmu,rho=Ccon)
        theta.x <- theta.x + ta
        d.x <- rweibull(1, shape=Ck, scale=Clam)/2.5
        x. <- x + c(d.x*cos(theta.x),d.x*sin(theta.x))
        
      }
      
      # Build trajectory
      x.t <- rbind(x.t,x.)
      x <- x.
      
    }
    
    # list of all elk trajectories
    preytrajs[[p]] <- data.frame(x.t,row.names=NULL)
    
  }
  
  p.tmp <- do.call("rbind", preytrajs)
  date <- rep(ld(simm.crw(1:(n+1)))$date, abPrey)
  prey.ids <- rep(paste("Prey", 1:abPrey, sep = ""), each = n+1)
  prey.trajs <- as.ltraj(p.tmp,date=date,id=prey.ids)
  
  # WOLF MOVEMENT
  
  for (P in 1:abPred){
    
    # random starting location and direction within territory
    y <- as.numeric(wolf.starts@coords[P,])
    theta.y <- runif(1,-pi,pi)
    y.t <- y
    
    for (i in 1:n){
      
      ysp <- SpatialPoints(matrix(y, nrow = 1))
      ypos <- over(ysp, arena)
      
      if (is.na(ypos[1,2]) == T){
        
        y. <- y.t[nrow(y.t)-1,]
        
      }
      
      else {
        
        ta <- rwrpcauchy(1,location=Smu,rho=Scon) 
        theta.y <- theta.y + ta
        d.y <- rweibull(1, shape=Sk, scale=Slam)
        y. <- y + c(d.y*cos(theta.y),d.y*sin(theta.y))
        
      }
      
      # build trajectory
      y.t <- rbind(y.t,y.)
      y <- y.
      
    }
    
    # list of all wolf trajectories
    predtrajs[[P]] <- data.frame(y.t,row.names=NULL)
    
  }
  
  P.tmp <- do.call("rbind", predtrajs)
  date <- rep(ld(simm.crw(1:(n+1)))$date, abPred)
  pred.ids <- rep(paste("Predator", 1:abPred, sep = ""), each = n+1)
  pred.trajs <- as.ltraj(P.tmp,date=date,id=pred.ids)
  
  return(c(prey.trajs, pred.trajs))
}

##### RPH MODEL #####

risky.places <- function(
  abPred=npred,                       #predator abundance 
  abPrey=nprey,                       #prey abundance
  n=nsteps,                           #number of movement steps per trajectory
  C2F=0.134,                          #commuting-->foraging transition probability for elk
  F2C=0.076,                          #foraging-->commuting transition probability for elk
  Fk=0.789,                           #shape parameter for elk foraging speed
  Flam=0.135,                         #scale parameter for elk foraging speed
  Fmu=0.057,                          #mean for elk foraging turning angle
  Fcon=0.414,                         #concentration for elk foraging turning angle
  Ck=0.809,                           #shape parameter for elk commuting speed
  Clam=0.343,                         #scale parameter for elk commuting speed
  Cmu=0.058,                          #mean for elk commuting turning angle
  Ccon=0.057,                         #concentration for elk commuting turning angle
  Sk=0.857,                           #shape parameter for wolf search speed
  Slam=0.856,                         #scale parameter for wolf search speed
  Smu=0.012,                          #mean for wolf search turning angle
  Scon=0.258                          #concentration for wolf search turning angle
){
  
  preytrajs <- list()
  predtrajs <- list()
  
  # WOLF MOVEMENT
  
  for (P in 1:abPred){
    
    # random starting location and direction within territory
    y <- as.numeric(wolf.starts@coords[P,])
    theta.y <- runif(1,-pi,pi)
    y.t <- y
    
    for (i in 1:n){
      
      ysp <- SpatialPoints(matrix(y, nrow = 1))
      ypos <- over(ysp, arena)
      
      if (is.na(ypos[1,2]) == T){
        
        y. <- y.t[nrow(y.t)-1,]
        
      }
      
      else{
        
        ta <- rwrpcauchy(1,location=Smu,rho=Scon)
        theta.y <- ta + theta.y
        d.y <- rweibull(1, shape=Sk, scale=Slam)
        y. <- y + c(d.y*cos(theta.y),d.y*sin(theta.y))
        
      }
      
      # build trajectory
      y.t <- rbind(y.t,y.)
      y <- y.
      
    }
    
    # list of all wolf trajectories
    predtrajs[[P]] <- data.frame(y.t,row.names=NULL)
    
  }
  
  P.tmp <- do.call("rbind", predtrajs)
  date <- rep(ld(simm.crw(1:(n+1)))$date, abPred)
  pred.ids <- rep(paste("Predator", 1:abPred, sep = ""), each = n+1)
  pred.trajs <- as.ltraj(P.tmp,date=date,id=pred.ids)
  
  predrelocs.df <- ld(pred.trajs)[,1:2]
  predrelocs.sp <- SpatialPoints(predrelocs.df[,1:2])
  predKDE <- kernelUD(predrelocs.sp, grid=UD.grid2(predrelocs.df))
  predcore <- getverticeshr(predKDE, percent=50)
  
  refugia <- arena-predcore
  
  # ELK MOVEMENT
  
  for (p in 1:abPrey){
    
    # random starting location and direction
    x <- as.numeric(spsample(refugia, n = 1, type = "random", iter=20)@coords)
    theta.x <- runif(1,-pi,pi)
    x.t <- x
    
    # random starting behaviour
    behaviours <- c("Foraging", "Commuting")
    behav <- sample(behaviours,1)
    behav2 <- behav
    
    for (i in 1:n){
      
      xsp <- SpatialPoints(matrix(x, nrow = 1))
      xpos <- over(xsp, predcore)
      
      if (is.na(xpos[1,2]) == F){
        
        ta <- rwrpcauchy(1,location=Cmu,rho=Ccon)
        theta.x <- ta + theta.x
        d.x <- rweibull(1, shape=Ck, scale=Clam)/2.5
        x. <- x + c(d.x*cos(theta.x),d.x*sin(theta.x))
        
      }
      
      else{
        # change of behaviour based on transition probabilities
        if (behav2 == "Foraging"){
          behav.change <- rbinom(1,1,F2C)
          if (behav.change == 1){ 
            behav <- "Commuting"
          }
        }
        
        if (behav2 == "Commuting") {
          behav.change <- rbinom(1,1,C2F)
          if (behav.change == 1){ 
            behav <- "Foraging"
          }
        }
        
        # Movement step
        if (behav == "Foraging"){
          
          ta <- rwrpcauchy(n=1,location=Fmu,rho=Fcon)
          theta.x <- ta + theta.x
          d.x <- rweibull(1, shape=Fk, scale=Flam)/2.5
          x. <- x + c(d.x*cos(theta.x),d.x*sin(theta.x))
          
        }
        
        if (behav == "Commuting"){
          
          ta <- rwrpcauchy(1,location=Cmu,rho=Ccon)
          theta.x <- ta + theta.x
          d.x <- rweibull(1, shape=Ck, scale=Clam)/2.5
          x. <- x + c(d.x*cos(theta.x),d.x*sin(theta.x))
          
        }
        
      }
      
      # Build trajectory
      x.t <- rbind(x.t,x.)
      x <- x.
      
    }
    
    # list of all elk trajectories
    preytrajs[[p]] <- data.frame(x.t,row.names=NULL)
    
  }
  
  p.tmp <- do.call("rbind", preytrajs)
  date <- rep(ld(simm.crw(1:(n+1)))$date, abPrey)
  prey.ids <- rep(paste("Prey", 1:abPrey, sep = ""), each = n+1)
  prey.trajs <- as.ltraj(p.tmp,date=date,id=prey.ids)
  
  return(c(prey.trajs, pred.trajs))
}

##### RTH MODEL #####

risky.times <- function(
  abPred=npred,                      #predator abundance 
  abPrey=nprey,                      #prey abundance
  n=nsteps,                          #number of movement steps per trajectory
  C2F=0.134,                          #commuting-->foraging transition probability for elk
  F2C=0.076,                          #foraging-->commuting transition probability for elk
  Fk=0.789,                           #shape parameter for elk foraging speed
  Flam=0.135,                         #scale parameter for elk foraging speed
  Fmu=0.057,                          #mean for elk foraging turning angle
  Fcon=0.414,                         #concentration for elk foraging turning angle
  Ck=0.809,                           #shape parameter for elk commuting speed
  Clam=0.343,                         #scale parameter for elk commuting speed
  Cmu=0.058,                          #mean for elk commuting turning angle
  Ccon=0.057,                         #concentration for elk commuting turning angle
  Sk=0.857,                           #shape parameter for wolf search speed
  Slam=0.856,                         #scale parameter for wolf search speed
  Smu=0.012,                          #mean for wolf search turning angle
  Scon=0.258,                         #concentration for wolf search turning angle
  dist=1                              #response distance for elk
){
  
  preytrajs <- list()
  predtrajs <- list()
  
  # WOLF MOVEMENT
  
  for (P in 1:abPred){
    
    # random starting location and direction
    y <- as.numeric(wolf.starts@coords[P,])
    y.t <- y
    
    for (i in 1:n){
      
      ysp <- SpatialPoints(matrix(y, nrow = 1))
      ypos <- over(ysp, arena)
      
      if (is.na(ypos[1,2]) == T){
        
        y. <- y.t[nrow(y.t)-1,]
        
      }
      
      else{
        
        theta.y <- rwrpcauchy(1,location=Smu,rho=Scon)*(180/pi)                     
        d.y <- rweibull(1, shape=Sk, scale=Slam)
        y. <- y + c(d.y*cos(theta.y),d.y*sin(theta.y))
        
      }
      
      # Build trajectory
      y.t <- rbind(y.t,y.)
      y <- y.
      
    }
    
    # list of all elk trajectories
    predtrajs[[P]] <- data.frame(y.t,row.names=NULL)
    
  }
  
  P.tmp <- do.call("rbind", predtrajs)
  date <- rep(ld(simm.crw(1:(n+1)))$date, abPred)
  pred.ids <- rep(paste("Predator", 1:abPred, sep = ""), each = n+1)
  pred.trajs <- as.ltraj(P.tmp,date=date,id=pred.ids)
  
  pflater <- ld(pred.trajs)
  step <- rep(1:(n+1), npred)
  pflater <- cbind(pflater, step)
  
  # ELK MOVEMENT
  
  for (p in 1:abPrey){
    
    # random starting location and direction
    x <- as.numeric(spsample(arena, n = 1, type = "random", iter=20)@coords)
    x.t <- x
    
    # random starting behaviour
    behaviours <- c("Foraging", "Commuting")
    behav <- sample(behaviours,1)
    behav2 <- behav
    
    for (i in 1:n){
      
      # extract all pred relocs for the same time step
      sim.pred <- subset(pflater, pflater$step >= i-24 & pflater$step <= i)[,1:2]
      sim.pred.pts <- SpatialPoints(sim.pred)
      x.sp <- SpatialPoints(matrix(x, ncol = 2))
      dists <- spDistsN1(sim.pred.pts, x.sp)
      min.dist <- min(dists)                                                
      min.ind <- which(dists == min(dists))
      
      if (min.dist < dist){
        
        next.p <- c(sim.pred[min.ind,1], sim.pred[min.ind,2])                                                                       
        psi <- atan2(next.p[2]-x[2],next.p[1]-x[1])
        if (psi < 2*pi){
          npsi <- psi+pi
        }
        if (psi > 2*pi){
          npsi <- psi-pi
        }
        
        theta.x <- npsi
        #theta.x <- rwrpcauchy(1,location=npsi,rho=Ccon)*(180/pi)                    
        d.x <- rweibull(1, shape=Ck, scale=Clam)
        x. <- x + c(d.x*cos(theta.x),d.x*sin(theta.x))
        
      }
      
      else {
        
        # change of behaviour based on transition probabilities
        if (behav2 == "Foraging"){
          behav.change <- rbinom(1,1,F2C)
          if (behav.change == 1){ 
            behav <- "Commuting"
          }
        }
        
        if (behav2 == "Commuting") {
          behav.change <- rbinom(1,1,C2F)
          if (behav.change == 1){ 
            behav <- "Foraging"
          }
        }
        
        # Movement step
        if (behav == "Foraging"){
          
          theta.x <- rwrpcauchy(n=1,location=Fmu,rho=Fcon)                  
          d.x <- rweibull(1, shape=Fk, scale=Flam)/2.5
          x. <- x + c(d.x*cos(theta.x),d.x*sin(theta.x))
          
        }
        
        if (behav == "Commuting"){
          
          theta.x <- rwrpcauchy(1,location=Cmu,rho=Ccon)*(180/pi)                    
          d.x <- rweibull(1, shape=Ck, scale=Clam)/2.5
          x. <- x + c(d.x*cos(theta.x),d.x*sin(theta.x))
          
        }
        
      }
      
      # build trajectory
      x.t <- rbind(x.t,x.)
      x <- x.
      
    }
    
    # list of all wolf trajectories
    preytrajs[[p]] <- data.frame(x.t,row.names=NULL)
    
  }
  
  p.tmp <- do.call("rbind", preytrajs)
  date <- rep(ld(simm.crw(1:(n+1)))$date, abPrey)
  prey.ids <- rep(paste("Prey", 1:abPrey, sep = ""), each = n+1)
  prey.trajs <- as.ltraj(p.tmp,date=date,id=prey.ids)
  
  return(c(prey.trajs, pred.trajs))
}

##### RAH MODEL #####

risk.allocation <- function(
  abPred=npred,                      #predator abundance 
  abPrey=nprey,                      #prey abundance
  n=nsteps,                          #number of movement steps per trajectory
  C2F=0.134,                          #commuting-->foraging transition probability for elk
  F2C=0.076,                          #foraging-->commuting transition probability for elk
  Fk=0.789,                           #shape parameter for elk foraging speed
  Flam=0.135,                         #scale parameter for elk foraging speed
  Fmu=0.057,                          #mean for elk foraging turning angle
  Fcon=0.414,                         #concentration for elk foraging turning angle
  Ck=0.809,                           #shape parameter for elk commuting speed
  Clam=0.343,                         #scale parameter for elk commuting speed
  Cmu=0.058,                          #mean for elk commuting turning angle
  Ccon=0.057,                         #concentration for elk commuting turning angle
  Sk=0.857,                           #shape parameter for wolf search speed
  Slam=0.856,                         #scale parameter for wolf search speed
  Smu=0.012,                          #mean for wolf search turning angle
  Scon=0.258,                         #concentration for wolf search turning angle
  dist=1,                             #response distance for elk
  tc1=24,                             #temporal scale of response in low background risk areas
  tc2=1                               #temporal scale of response in high background risk areas
){
  
  preytrajs <- list()
  predtrajs <- list()
  
  # WOLF MOVEMENT
  
  for (P in 1:abPred){
    
    # random starting location and direction within territory
    y <- as.numeric(spsample(arena, n = 1, type = "random", iter=20)@coords)
    y.t <- y
    
    for (i in 1:n){
      
      ysp <- SpatialPoints(matrix(y, nrow = 1))
      ypos <- over(ysp, arena)
      
      if (is.na(ypos[1,2]) == T){
        
        y. <- y.t[nrow(y.t)-1,]
        
      }
      
      else{
        
        theta.y <- rwrpcauchy(1,location=Smu,rho=Scon)*(180/pi)                     
        d.y <- rweibull(1, shape=Sk, scale=Slam)
        y. <- y + c(d.y*cos(theta.y),d.y*sin(theta.y))
        
      }
      
      # build trajectory
      y.t <- rbind(y.t,y.)
      y <- y.
      
    }
    
    # list of all wolf trajectories
    predtrajs[[P]] <- data.frame(y.t,row.names=NULL)
    
  }
  
  P.tmp <- do.call("rbind", predtrajs)
  date <- rep(ld(simm.crw(1:(n+1)))$date, abPred)
  pred.ids <- rep(paste("Predator", 1:abPred, sep = ""), each = n+1)
  pred.trajs <- as.ltraj(P.tmp,date=date,id=pred.ids)
  
  pflater <- ld(pred.trajs)
  step <- rep(1:(n+1), npred)
  pflater <- cbind(pflater, step)
  
  predrelocs.df <- ld(pred.trajs)[,1:2]
  predrelocs.sp <- SpatialPoints(predrelocs.df[,1:2])
  predKDE <- kernelUD(predrelocs.sp, grid=UD.grid2(predrelocs.df))
  predcore <- getverticeshr(predKDE, percent=50)
  
  # ELK MOVEMENT
  
  for (p in 1:abPrey){
    
    # random starting location and direction
    x <- as.numeric(spsample(arena, n = 1, type = "random", iter=20)@coords)
    x.t <- x
    
    # random starting behaviour
    behaviours <- c("Foraging", "Commuting")
    behav <- sample(behaviours,1)
    behav2 <- behav
    
    for (i in 1:n){
      
      xsp <- SpatialPoints(matrix(x, nrow = 1))
      xpos.wrt.predcore <- over(xsp, predcore)
      
      if (is.na(xpos.wrt.predcore[1,2]) == F){
        
        # extract all pred relocs for the same time step
        sim.pred <- subset(pflater, pflater$step >= i-tc2 & pflater$step <= i)[,1:2]
        sim.pred.pts <- SpatialPoints(sim.pred)
        x.sp <- SpatialPoints(matrix(x, ncol = 2))
        dists <- spDistsN1(sim.pred.pts, x.sp)
        min.dist <- min(dists)                                                
        min.ind <- which(dists == min(dists))
        
        if (min.dist <= dist){
          
          next.p <- c(sim.pred[min.ind,1], sim.pred[min.ind,2])                                                                       
          psi <- atan2(next.p[2]-x[2],next.p[1]-x[1])
          if (psi < 2*pi){
            npsi <- psi+pi
          }
          if (psi > 2*pi){
            npsi <- psi-pi
          }
          
          theta.x <- npsi
          #theta.x <- rwrpcauchy(1,location=npsi,rho=Ccon)*(180/pi)                    
          d.x <- rweibull(1, shape=Ck, scale=Clam)/2.5
          x. <- x + c(d.x*cos(theta.x),d.x*sin(theta.x))
          
        }
        
        if (min.dist > dist){
          
          # change of behaviour based on transition probabilities
          if (behav2 == "Foraging"){
            behav.change <- rbinom(1,1,F2C)
            if (behav.change == 1){ 
              behav <- "Commuting"
            }
          }
          
          if (behav2 == "Commuting") {
            behav.change <- rbinom(1,1,C2F)
            if (behav.change == 1){ 
              behav <- "Foraging"
            }
          }
          
          # Movement step
          if (behav == "Foraging"){
            
            theta.x <- rwrpcauchy(n=1,location=Fmu,rho=Fcon)*(180/pi)                  
            d.x <- rweibull(1, shape=Fk, scale=Flam)/2.5
            x. <- x + c(d.x*cos(theta.x),d.x*sin(theta.x))
            
          }
          
          if (behav == "Commuting"){
            
            theta.x <- rwrpcauchy(1,location=Cmu,rho=Ccon)*(180/pi)                    
            d.x <- rweibull(1, shape=Ck, scale=Clam)/2.5
            x. <- x + c(d.x*cos(theta.x),d.x*sin(theta.x))
            
          }
          
        }
        
      }
      
      else {
        
        # extract all pred relocs for 24 hours before
        sim.pred <- subset(pflater, pflater$step >= i-tc1 & pflater$step <= i)[,1:2]
        sim.pred.pts <- SpatialPoints(sim.pred)
        x.sp <- SpatialPoints(matrix(x, ncol = 2))
        dists <- spDistsN1(sim.pred.pts, x.sp)
        min.dist <- min(dists)                                                
        min.ind <- which(dists == min(dists))
        
        if (min.dist <= dist){
          
          next.p <- c(sim.pred[min.ind,1], sim.pred[min.ind,2])                                                                       
          psi <- atan2(next.p[2]-x[2],next.p[1]-x[1])
          if (psi < 2*pi){
            npsi <- psi+pi
          }
          if (psi > 2*pi){
            npsi <- psi-pi
          }
          
          theta.x <- npsi
          #theta.x <- rwrpcauchy(1,location=npsi,rho=Ccon)*(180/pi)                    
          d.x <- rweibull(1, shape=Ck, scale=Clam)/2.5
          x. <- x + c(d.x*cos(theta.x),d.x*sin(theta.x))
          
        }
        
        if (min.dist > dist){
          
          # change of behaviour based on transition probabilities
          if (behav2 == "Foraging"){
            behav.change <- rbinom(1,1,F2C)
            if (behav.change == 1){ 
              behav <- "Commuting"
            }
          }
          
          if (behav2 == "Commuting") {
            behav.change <- rbinom(1,1,C2F)
            if (behav.change == 1){ 
              behav <- "Foraging"
            }
          }
          
          # Movement step
          if (behav == "Foraging"){
            
            theta.x <- rwrpcauchy(n=1,location=Fmu,rho=Fcon)*(180/pi)                  
            d.x <- rweibull(1, shape=Fk, scale=Flam)/2.5
            x. <- x + c(d.x*cos(theta.x),d.x*sin(theta.x))
            
          }
          
          if (behav == "Commuting"){
            
            theta.x <- rwrpcauchy(1,location=Cmu,rho=Ccon)*(180/pi)                    
            d.x <- rweibull(1, shape=Ck, scale=Clam)/2.5
            x. <- x + c(d.x*cos(theta.x),d.x*sin(theta.x))
            
          }
          
        }
        
      }
      
      # Build trajectory
      x.t <- rbind(x.t,x.)
      x <- x.
      
    }
    
    # list of all elk trajectories
    preytrajs[[p]] <- data.frame(x.t,row.names=NULL)
    
  }
  
  p.tmp <- do.call("rbind", preytrajs)
  date <- rep(ld(simm.crw(1:(n+1)))$date, abPrey)
  prey.ids <- rep(paste("Prey", 1:abPrey, sep = ""), each = n+1)
  prey.trajs <- as.ltraj(p.tmp,date=date,id=prey.ids)
  
  return(c(prey.trajs, pred.trajs))
}






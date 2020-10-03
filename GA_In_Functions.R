

#*************************************************
#**************GA by Sara Mohammadi***************
#*************************************************

GeneticAlgo <- function(filePath, pop.size,iteration_num,exePath,resultPath){
  set.seed(1234)
  library("dplyr")
  library("readxl")
  library("gdata")
  library(compare)
  setwd(exePath)
  mutation.prob = 0.03
  tblGL <- read.csv(filePath)
  dfGL <- data.frame(tblGL)
  pre.dfGL2 <- dfGL[, -1]
  um.cols <-c("Liq1DensityGL","GasDensityGL","Liq1ViscosityGL","GasViscosityGL","GasLiqSurTensionGL","PipeDiaGL","PipeIncliGL","PipeRoughGL","Liq1SuperficialVelGL","GasSuperficialVelGL","PressureDropGL")
  inval.rows <- get.invalid.rows(pre.dfGL2[,um.cols])
  if (length(inval.rows) > 0) {
    dfGL2 <- pre.dfGL2[-inval.rows,]
  } else {
    dfGL2 <- pre.dfGL2
  }
  # dim(dfGL2)
  # dim(pre.dfGL2)
  #*****************************************************************
  n <- floor((1 + sqrt(1 + 4 * (0.5 * pop.size))) / 2)
  if (n >= 3) {
    direct.n <- floor(n / 3)
  } else
    direct.n <- n
  #*****************************************************************
  N <- floor(1 + sqrt(1 + 4 * 0.2 * pop.size) / 2)
  
  generation.colnames <- vector()
  for (i in 1:pop.size) {
    generation.colnames[i] <- paste("chromosome", i, sep = "")
  }
  is.init.generation <- 1
  
  file.name <-  tools::file_path_sans_ext(basename(filePath))
  #file.name <-  substr(start =1 , stop = nchar(fileName)-4, x = fileName)
  pre_name <- paste(file.name, paste(pop.size, "pop", sep = ""), sep = "_")
  pre_name2 <- paste(pre_name, mutation.prob, sep = "_")
  pre_name3 <- paste(resultPath,pre_name2,sep = "/")
  f.path <- paste(paste(pre_name3, "m", sep = ""), ".txt", sep = "")
  #***************************************************************************************
  #ls.fitness <- list()
  for (p in 1:iteration_num) {
  # stop.ga <- FALSE
  # p = 1
  # while(!stop.ga){
    if (is.init.generation) {
      if (file.exists(f.path)) file.remove(f.path)
      generation <- replicate(pop.size, create.random.genes())
      initial.generation <- data.frame(generation,row.names = c("gene1","gene2","gene3","gene4","gene5","gene6","gene7","gene8","gene9"))
      colnames(initial.generation) <- generation.colnames
      GA.fitness <- sapply(as.list(initial.generation),  run.unifiedmodel,filterTable=dfGL2)
      gene.fit.df <- as.data.frame(rbind(initial.generation, t(GA.fitness)))
      row.names(gene.fit.df) <- c("gene1","gene2","gene3","gene4","gene5","gene6","gene7","gene8","gene9","fitness")
      is.init.generation <- 0
    } else{
      next.generation <- as.data.frame(matrix(NA, ncol = pop.size, nrow = 9))
      rownames(next.generation) <-c("gene1","gene2","gene3","gene4","gene5","gene6","gene7","gene8","gene9")
      colnames(next.generation) <- generation.colnames
      
      #**************************************************************
      #*******n best chromosomes that have the min fitness value*****
      #**************************************************************
    
      # Roulette Wheel Selection
      # r.w.s <- rws(GA.fitness, n-1)
      # r.w <- c(r.w.s, as.integer(which.min(GA.fitness)))
      # next.generation[, 1:direct.n] <- initial.generation[, r.w[1:direct.n]]
      # temp.generation <- as.data.frame(matrix(NA, ncol = length(r.w), nrow = 9), row.names = NULL)
      # temp.generation <- initial.generation[, r.w]
      # 
      # # n.min <- sort(unique(GA.fitness), partial = 1:n)[1:n]
      # # n.min.int <- match(n.min, GA.fitness)
      # # next.generation[, 1:direct.n] <- initial.generation[, n.min.int[1:direct.n]]
      # # temp.generation <- as.data.frame(matrix(NA, ncol = length(n.min.int), nrow = 9), row.names = NULL)
      # # temp.generation <- initial.generation[, n.min.int]
      # crossover.result <- list()
      # for (i in 1:(length(r.w) - 1)) {
      #   crossover.result[[i]] <- apply(temp.generation,2,crossover.chromosome,temp.generation[, 1])
      #   if (i < length(r.w))
      #   {
      #     temp.generation <- temp.generation[,-1]
      #   }
      # }
      # df <- as.data.frame(matrix(data = unlist(crossover.result) ,ncol = length(unlist(crossover.result)) / 9,nrow = 9))
      #new.df <- df[, !(df %in% initial.generation[, n.min.int])]
      n.min <- sort(unique(GA.fitness), partial = 1:n)[1:n]
      n.min.int <- match(n.min, GA.fitness)
      next.generation[, 1:direct.n] <- initial.generation[, n.min.int[1:direct.n]]
      temp.generation <- as.data.frame(matrix(NA, ncol = length(n.min.int), nrow = 9), row.names = NULL)
      temp.generation <- initial.generation[, n.min.int]
      crossover.result <- list()
      for (i in 1:(length(n.min.int) - 1)) {
        crossover.result[[i]] <- apply(temp.generation,2,crossover.chromosome,temp.generation[, 1])
        if (i < length(n.min.int))
        {
          temp.generation <- temp.generation[,-1]
        }
      }
      df <- as.data.frame(matrix(data = unlist(crossover.result) ,ncol = length(unlist(crossover.result)) / 9,nrow = 9))
      #***************************************************
      #***************************************************
      v.extra.index <- vector()
      v.extra.index2 <- vector()
      if (n > 4) {
        v.extra.index[1] <- 1
        v.extra.index[2] <- v.extra.index[1] + 2 * n
        v.extra.index[3] <- v.extra.index[2] + 2 * (n - 1)
        v.extra.index2 <-
          c(v.extra.index[1] + 1, v.extra.index[2] + 1, v.extra.index[3] + 1)
        for (i in 4:(n - 1)) {
          v.extra.index[i] = v.extra.index[i - 1] + 2 * (n - i + 1)
          v.extra.index2[i] = v.extra.index[i] + 1
        }
      } else if (n == 2) {
        v.extra.index[1] <- 1
        v.extra.index2 <- v.extra.index[1] + 1
      } else if (n == 3) {
        v.extra.index[1] <- 1
        v.extra.index[2] <- v.extra.index[1] + 2 * n
        v.extra.index2 <- c(v.extra.index[1] + 1, v.extra.index[2] + 1)
      }
      v.extra.index
      v.extra.index2
      final.v <- sort(c(v.extra.index, v.extra.index2))
      new.df <- df[, setdiff(c(1:(n * (n + 1) - 2)), final.v)]
      
      #************************************************************
      #**************************Mutation**************************
      #************************************************************
      # mutation.prob <- 1/pop.size
      # mutation.prob <- 0.03
      rand.0.1 <- runif(1, 0, 1)
      if (rand.0.1 < mutation.prob)
      {
        chr.to.mutate <- sample(1:dim(new.df)[2], 1)
        gene.to.mutate <- sample(1:9, 1)
        new.gene.val <-
          mutation(df.mutate = new.df, index.gene = gene.to.mutate)
        new.df[gene.to.mutate, chr.to.mutate] <- new.gene.val
      }
      t.unique.new.df <- t(new.df)[!duplicated(t(new.df)), ]
      unique.new.df <- t(t.unique.new.df)
      
      #********************************
      if(is.null(dim(unique.new.df))){
        next.generation[,(direct.n + 1):2] <- unique.new.df
      } else {
        colNumbers <- (direct.n + 1):(dim(unique.new.df)[2] + direct.n )
        next.generation[,colNumbers] <- unique.new.df
        #   next.generation[,(direct.n + 1):(dim(unique.new.df)[2] + 1)] <- unique.new.df
      }
      dim(unique.new.df)[2]== n*(n-1) # number of crossover result of n  must be equal to length of
      rand.from.intial <- sample(setdiff((1:pop.size), n.min.int), N)
      N.crossover <- initial.generation[, rand.from.intial]
      new.cross.over <-
        apply(N.crossover, 2, crossover.chromosome, initial.generation[, n.min.int[length(n.min.int)]])
      N.df <- as.data.frame(matrix(data = unlist(new.cross.over),ncol =  N * 2 ,nrow = 9))
      #we do not want to trasfer all the crossover results
      odd.cols <- seq(1, ncol(N.df), by = 2)
      t.unique.N.df <- t(N.df[, odd.cols])[!duplicated(t(N.df[, odd.cols])), ]
      unique.N.df <- t(t.unique.N.df)
      next.generation[, (dim(unique.new.df)[2] + 2):(dim(unique.new.df)[2] + 1 + length(odd.cols))] <- N.df[, odd.cols]
      temp.gen <-t(t(next.generation)[!duplicated(t(next.generation)), ])
      temp.gen <-temp.gen[, colnames(temp.gen)[colSums(!is.na(temp.gen)) > 0]]
      next.generation[, 1:dim(temp.gen)[2]] <- temp.gen
      next.generation[, (dim(temp.gen)[2] + 1):pop.size] <- NA
      
      #*****************************************************************************************
      #*******Number of columns in next.generation that still are empty and we need to fill*****
      #*****************************************************************************************
      if (sum(is.na(next.generation) != 0)) {
        empty.col.size <- sum(is.na(next.generation)) / 9
        index.no <- (pop.size - empty.col.size) + 1
        next.generation[, index.no:pop.size] <-
          replicate(empty.col.size, create.random.genes())
      }
      GA.fitness <- sapply(as.list(next.generation),  run.unifiedmodel,filterTable=dfGL2)
      initial.generation <- next.generation
      out.result <- vector()
      out.result <- replicate(direct.n, "direct")
      cross.len <- length((direct.n + 1):(n * (n - 1) + 1))
      out.result[(direct.n + 1):(direct.n + cross.len)] <- replicate(cross.len, "crossover_1st")
      out.result[(direct.n + cross.len + 1):(direct.n + cross.len + N)] <- replicate(N, "crossover_2nd")
      out.result[(direct.n + cross.len + N):pop.size] <- replicate(length((direct.n + cross.len + N):pop.size), "Random_genes")
      generation.fit <- rbind(initial.generation, fitness = GA.fitness)
      gene.fit.df <- as.data.frame(rbind(generation.fit, result_of = out.result))
      t.gene.fit.df <- as.data.frame(t(gene.fit.df))
      
      write.table(t(t.gene.fit.df[order(t.gene.fit.df$fitness,decreasing = FALSE), ]),
                  f.path,sep = "\t",append = TRUE,col.names = FALSE,row.names = FALSE)
      #******************************************
      is.init.generation <- 0
    }
  } # ******** End of MAIN for loop ******************
  #write.table((ls.fitness),f.path, col.names = FALSE,row.names = FALSE)
} # ********** End of Function GeneticAlgo ***********

#**************************************************
#**************************************************
#**************************************************



#*************************************************
#**************GA by Sara Mohammadi***************
#*************************************************

GeneticAlgo <- function(filePath, pop.size,iteration_num,exePath,resultPath){
  set.seed(1234)
  library("dplyr")
  library("readxl")
  library("gdata")
  library(compare)
  setwd(exePath)
  mutation.prob = 0.03
  tblGL <- read.csv(filePath)
  dfGL <- data.frame(tblGL)
  pre.dfGL2 <- dfGL[, -1]
  um.cols <-c("Liq1DensityGL","GasDensityGL","Liq1ViscosityGL","GasViscosityGL","GasLiqSurTensionGL","PipeDiaGL","PipeIncliGL","PipeRoughGL","Liq1SuperficialVelGL","GasSuperficialVelGL","PressureDropGL")
  inval.rows <- get.invalid.rows(pre.dfGL2[,um.cols])
  if (length(inval.rows) > 0) {
    dfGL2 <- pre.dfGL2[-inval.rows,]
  } else {
    dfGL2 <- pre.dfGL2
  }
  # dim(dfGL2)
  # dim(pre.dfGL2)
  #*****************************************************************
  n <- floor((1 + sqrt(1 + 4 * (0.5 * pop.size))) / 2)
  if (n >= 3) {
    direct.n <- floor(n / 3)
  } else
    direct.n <- n
  #*****************************************************************
  N <- floor(1 + sqrt(1 + 4 * 0.2 * pop.size) / 2)
  
  generation.colnames <- vector()
  for (i in 1:pop.size) {
    generation.colnames[i] <- paste("chromosome", i, sep = "")
  }
  is.init.generation <- 1
  
  file.name <-  tools::file_path_sans_ext(basename(filePath))
  #file.name <-  substr(start =1 , stop = nchar(fileName)-4, x = fileName)
  pre_name <- paste(file.name, paste(pop.size, "pop", sep = ""), sep = "_")
  pre_name2 <- paste(pre_name, mutation.prob, sep = "_")
  pre_name3 <- paste(resultPath,pre_name2,sep = "/")
  f.path <- paste(paste(pre_name3, "m", sep = ""), ".txt", sep = "")
  #***************************************************************************************
  #ls.fitness <- list()
  for (p in 1:iteration_num) {
    # stop.ga <- FALSE
    # p = 1
    # while(!stop.ga){
    if (is.init.generation) {
      if (file.exists(f.path)) file.remove(f.path)
      generation <- replicate(pop.size, create.random.genes())
      initial.generation <- data.frame(generation,row.names = c("gene1","gene2","gene3","gene4","gene5","gene6","gene7","gene8","gene9"))
      colnames(initial.generation) <- generation.colnames
      GA.fitness <- sapply(as.list(initial.generation),  run.unifiedmodel,filterTable=dfGL2)
      gene.fit.df <- as.data.frame(rbind(initial.generation, t(GA.fitness)))
      row.names(gene.fit.df) <- c("gene1","gene2","gene3","gene4","gene5","gene6","gene7","gene8","gene9","fitness")
      is.init.generation <- 0
    } else{
      next.generation <- as.data.frame(matrix(NA, ncol = pop.size, nrow = 9))
      rownames(next.generation) <-c("gene1","gene2","gene3","gene4","gene5","gene6","gene7","gene8","gene9")
      colnames(next.generation) <- generation.colnames
      
      #**************************************************************
      #*******n best chromosomes that have the min fitness value*****
      #**************************************************************
      
      # Roulette Wheel Selection
      # r.w.s <- rws(GA.fitness, n-1)
      # r.w <- c(r.w.s, as.integer(which.min(GA.fitness)))
      # next.generation[, 1:direct.n] <- initial.generation[, r.w[1:direct.n]]
      # temp.generation <- as.data.frame(matrix(NA, ncol = length(r.w), nrow = 9), row.names = NULL)
      # temp.generation <- initial.generation[, r.w]
      # 
      # # n.min <- sort(unique(GA.fitness), partial = 1:n)[1:n]
      # # n.min.int <- match(n.min, GA.fitness)
      # # next.generation[, 1:direct.n] <- initial.generation[, n.min.int[1:direct.n]]
      # # temp.generation <- as.data.frame(matrix(NA, ncol = length(n.min.int), nrow = 9), row.names = NULL)
      # # temp.generation <- initial.generation[, n.min.int]
      # crossover.result <- list()
      # for (i in 1:(length(r.w) - 1)) {
      #   crossover.result[[i]] <- apply(temp.generation,2,crossover.chromosome,temp.generation[, 1])
      #   if (i < length(r.w))
      #   {
      #     temp.generation <- temp.generation[,-1]
      #   }
      # }
      # df <- as.data.frame(matrix(data = unlist(crossover.result) ,ncol = length(unlist(crossover.result)) / 9,nrow = 9))
      #new.df <- df[, !(df %in% initial.generation[, n.min.int])]
      n.min <- sort(unique(GA.fitness), partial = 1:n)[1:n]
      n.min.int <- match(n.min, GA.fitness)
      next.generation[, 1:direct.n] <- initial.generation[, n.min.int[1:direct.n]]
      temp.generation <- as.data.frame(matrix(NA, ncol = length(n.min.int), nrow = 9), row.names = NULL)
      temp.generation <- initial.generation[, n.min.int]
      crossover.result <- list()
      for (i in 1:(length(n.min.int) - 1)) {
        crossover.result[[i]] <- apply(temp.generation,2,crossover.chromosome,temp.generation[, 1])
        if (i < length(n.min.int))
        {
          temp.generation <- temp.generation[,-1]
        }
      }
      df <- as.data.frame(matrix(data = unlist(crossover.result) ,ncol = length(unlist(crossover.result)) / 9,nrow = 9))
      #***************************************************
      #***************************************************
      v.extra.index <- vector()
      v.extra.index2 <- vector()
      if (n > 4) {
        v.extra.index[1] <- 1
        v.extra.index[2] <- v.extra.index[1] + 2 * n
        v.extra.index[3] <- v.extra.index[2] + 2 * (n - 1)
        v.extra.index2 <-
          c(v.extra.index[1] + 1, v.extra.index[2] + 1, v.extra.index[3] + 1)
        for (i in 4:(n - 1)) {
          v.extra.index[i] = v.extra.index[i - 1] + 2 * (n - i + 1)
          v.extra.index2[i] = v.extra.index[i] + 1
        }
      } else if (n == 2) {
        v.extra.index[1] <- 1
        v.extra.index2 <- v.extra.index[1] + 1
      } else if (n == 3) {
        v.extra.index[1] <- 1
        v.extra.index[2] <- v.extra.index[1] + 2 * n
        v.extra.index2 <- c(v.extra.index[1] + 1, v.extra.index[2] + 1)
      }
      v.extra.index
      v.extra.index2
      final.v <- sort(c(v.extra.index, v.extra.index2))
      new.df <- df[, setdiff(c(1:(n * (n + 1) - 2)), final.v)]
      
      #************************************************************
      #**************************Mutation**************************
      #************************************************************
      # mutation.prob <- 1/pop.size
      # mutation.prob <- 0.03
      rand.0.1 <- runif(1, 0, 1)
      if (rand.0.1 < mutation.prob)
      {
        chr.to.mutate <- sample(1:dim(new.df)[2], 1)
        gene.to.mutate <- sample(1:9, 1)
        new.gene.val <-
          mutation(df.mutate = new.df, index.gene = gene.to.mutate)
        new.df[gene.to.mutate, chr.to.mutate] <- new.gene.val
      }
      t.unique.new.df <- t(new.df)[!duplicated(t(new.df)), ]
      unique.new.df <- t(t.unique.new.df)
      
      #********************************
      if(is.null(dim(unique.new.df))){
        next.generation[,(direct.n + 1):2] <- unique.new.df
      } else {
        colNumbers <- (direct.n + 1):(dim(unique.new.df)[2] + direct.n )
        next.generation[,colNumbers] <- unique.new.df
        #   next.generation[,(direct.n + 1):(dim(unique.new.df)[2] + 1)] <- unique.new.df
      }
      dim(unique.new.df)[2]== n*(n-1) # number of crossover result of n  must be equal to length of
      rand.from.intial <- sample(setdiff((1:pop.size), n.min.int), N)
      N.crossover <- initial.generation[, rand.from.intial]
      new.cross.over <-
        apply(N.crossover, 2, crossover.chromosome, initial.generation[, n.min.int[length(n.min.int)]])
      N.df <- as.data.frame(matrix(data = unlist(new.cross.over),ncol =  N * 2 ,nrow = 9))
      #we do not want to trasfer all the crossover results
      odd.cols <- seq(1, ncol(N.df), by = 2)
      t.unique.N.df <- t(N.df[, odd.cols])[!duplicated(t(N.df[, odd.cols])), ]
      unique.N.df <- t(t.unique.N.df)
      next.generation[, (dim(unique.new.df)[2] + 2):(dim(unique.new.df)[2] + 1 + length(odd.cols))] <- N.df[, odd.cols]
      temp.gen <-t(t(next.generation)[!duplicated(t(next.generation)), ])
      temp.gen <-temp.gen[, colnames(temp.gen)[colSums(!is.na(temp.gen)) > 0]]
      next.generation[, 1:dim(temp.gen)[2]] <- temp.gen
      next.generation[, (dim(temp.gen)[2] + 1):pop.size] <- NA
      
      #*****************************************************************************************
      #*******Number of columns in next.generation that still are empty and we need to fill*****
      #*****************************************************************************************
      if (sum(is.na(next.generation) != 0)) {
        empty.col.size <- sum(is.na(next.generation)) / 9
        index.no <- (pop.size - empty.col.size) + 1
        next.generation[, index.no:pop.size] <-
          replicate(empty.col.size, create.random.genes())
      }
      GA.fitness <- sapply(as.list(next.generation),  run.unifiedmodel,filterTable=dfGL2)
      initial.generation <- next.generation
      out.result <- vector()
      out.result <- replicate(direct.n, "direct")
      cross.len <- length((direct.n + 1):(n * (n - 1) + 1))
      out.result[(direct.n + 1):(direct.n + cross.len)] <- replicate(cross.len, "crossover_1st")
      out.result[(direct.n + cross.len + 1):(direct.n + cross.len + N)] <- replicate(N, "crossover_2nd")
      out.result[(direct.n + cross.len + N):pop.size] <- replicate(length((direct.n + cross.len + N):pop.size), "Random_genes")
      generation.fit <- rbind(initial.generation, fitness = GA.fitness)
      gene.fit.df <- as.data.frame(rbind(generation.fit, result_of = out.result))
      t.gene.fit.df <- as.data.frame(t(gene.fit.df))
      
      write.table(t(t.gene.fit.df[order(t.gene.fit.df$fitness,decreasing = FALSE), ]),
                  f.path,sep = "\t",append = TRUE,col.names = FALSE,row.names = FALSE)
      #******************************************
      is.init.generation <- 0
    }
  } # ******** End of MAIN for loop ******************
  #write.table((ls.fitness),f.path, col.names = FALSE,row.names = FALSE)
} # ********** End of Function GeneticAlgo ***********

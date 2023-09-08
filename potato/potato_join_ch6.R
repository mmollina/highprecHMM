require(mappoly)
require(tidyverse)
require(ggsci)
require(highprecHMM)
ploidy <- 4

#### Reading multi-parental data ####
input.data <- read.csv("~/repos/highprecHMM/potato/TableS2_dose.csv")
input.data <- input.data[input.data$Chrom == 6, ]
pedigree1 <- read.csv("~/repos/highprecHMM/potato/TableS3_ped.csv")
pedigree <- pedigree1[,c(1,3,4)]
pedigree[1:3, 2:3][]<-NA
colnames(pedigree) <- c("Name", "Parent1", "Parent2")
colnames(input.data) <- stringr::str_replace_all(string = colnames(input.data), pattern = "\\.", replacement = "-")
colnames(input.data)[1:3] <- c("marker", "chromosome", "position")

input.data[1:3, 1:10]
pedigree[1:6,]

dat <- full_pop_to_full_sib(input.data, pedigree, ploidy)
names(dat)

#### Full-sib maps ####
fake.dat = dat[[1]]
system.time(map1 <- mapping_test(fake.dat, P1 = "W6511-1R", P2 = "VillettaRose"))
## Comparing with PolyOrigin paper: OK
#map1.err <- est_full_hmm_with_global_error(map1, error = 0.02, verbose = TRUE)
#plot(map1.err)
#g<-calc_genoprob_error(map1.err, error = 0.02)
#a<-calc_homoprob(g)
#plot(a, ind = "W15268.27R")

fake.dat = dat[[2]]
system.time(map2 <- mapping_test(fake.dat, P1 = "W6511-1R", P2 = "W9914-1R"))
fake.dat = dat[[3]]
system.time(map3 <- mapping_test(fake.dat, P1 = "VillettaRose", P2 = "W9914-1R"))

MAPs <- list("W6511-1RxVillettaRose" = map1, "W6511-1RxW9914-1R" = map2, "VillettaRosexW9914-1R" = map3)
plot_map_list(MAPs)
A <- cor_homolog(MAPs, pop1 = 1, pop2 = 2, Ppop1 = 1, Ppop2 = 1)
B <- cor_homolog(MAPs, pop1 = 1, pop2 = 3, Ppop1 = 2, Ppop2 = 1)
C <- cor_homolog(MAPs, pop1 = 2, pop2 = 3, Ppop1 = 2, Ppop2 = 2)

#### Building pedigree matrix ####
rownames(input.data) <- input.data[,"marker"]
input.data[1:3, 1:8]
rownames(pedigree) <- pedigree[,1]
ind.vec <- pedigree$Name[-c(1:3)]
z <- pedigree1$population
z <- as.numeric(stringr::str_remove_all(z, "pop"))
names(z) <- pedigree1$individual
pedigree$pop <- z[pedigree$Name]
head(pedigree)
#### States to visit ####
ngam <- choose(ploidy, ploidy/2)
A<-as.matrix(expand.grid(0:(ngam-1),
                         0:(ngam-1))[,2:1])
M <- vector("list", length(ind.vec))
names(M) <- ind.vec
for(j in ind.vec){
  cur.ind <- j
  cat(j, "\n")
  pop.id <- pedigree[cur.ind, "pop"]
  P1 <- MAPs[[pop.id]]$maps[[1]]$seq.ph$P
  P2 <- MAPs[[pop.id]]$maps[[1]]$seq.ph$Q
  P1M <- ph_list_to_matrix(P1, ploidy)
  P2M <- ph_list_to_matrix(P2, ploidy)
  rownames(P1M) <- rownames(P2M) <- MAPs[[pop.id]]$info$mrk.names
  d <- input.data[,cur.ind]
  names(d) <- rownames(input.data)
  d[setdiff(names(d), rownames(P1M))] <- NA
  I <- vector("list", length(d))
  names(I) <- names(d)
  for(i in rownames(P1M)){
    a <- kronecker(apply(combn(P1M[i,], ploidy/2), 2, sum),
                   apply(combn(P2M[i,], ploidy/2), 2, sum), "+")
    I[[i]] <- A[d[i] == a, , drop = FALSE]
  }
  for(i in which(sapply(I, is.null))){
    I[[i]] <- A
  }
  M[[j]] <- I
}
A <- vector("list", length(ind.vec))
names(A) <- ind.vec
H <- vector("list", nrow(input.data))
names(H) <- rownames(input.data)
for(i in names(H))
  H[[i]] <- A
for(i in names(H))
  for(j in names(A))
    H[[i]][[j]] <- M[[j]][[i]]

#### Joint mapping - no error modeling ####
map <- mappoly:::est_haplo_hmm(m = ploidy,
                               n.mrk = length(H),
                               n.ind = length(H[[1]]),
                               haplo = H,
                               rf_vec = rep(0.01, length(H)-1),
                               verbose = TRUE,
                               use_H0 = FALSE,
                               tol = 10e-4)
pos.map1 <- cumsum(c(0, imf_h(map[[2]])))
plot(pos.map1)
#### Joint mapping - error modeling ####
A<-as.matrix(expand.grid(0:(ngam-1),
                         0:(ngam-1))[,2:1])
B<-apply(A, 1, paste0, collapse = "-")
err <- 0.05
##  H[[mrk]]     H[[mrk]][[ind]]  H[[mrk]][[ind]][state,]
log_E <- numeric(length(H) * length(H[[1]]) * nrow(A))
count <- 1
for(j in 1:length(H[[1]])){
  for(i in 1:length(H)){
    id <- match(apply(H[[i]][[j]], 1, paste0, collapse = "-"), B)
    z <- rep(NA, ngam^2)
    z[id] <- (1 - err)/length(id)
    z[is.na(z)] <- err/(ngam^2 - length(id))
    log_E[count:(count+nrow(A)-1)] <- log(z)
    count <- count + nrow(A)
  }
}
n.mrk <- length(H)
n.ind <- length(H[[1]])
rf.vec <- rep(0.001, n.mrk-1)
system.time(map.err <- est_map_R(m = 4,
                                 n.mrk = n.mrk,
                                 n.ind = n.ind,
                                 emit = log_E,
                                 rf_vec = rf.vec,
                                 tol = 10e-4))
pos.map2 <- cumsum(c(0, imf_h(map.err[[2]])))
plot(pos.map2)
plot(pos.map1 ~ pos.map2)

E<-exp(log_E)
system.time(map.err2 <- est_map_R_no_log(m = 4,
                                 n.mrk = n.mrk,
                                 n.ind = n.ind,
                                 emit = E,
                                 rf_vec = rf.vec,
                                 tol = 10e-4))
pos.map3 <- cumsum(c(0, imf_h(map.err2[[2]])))
plot(pos.map3)
plot(pos.map2 ~ pos.map3)


#### Computation of genotype probability ####
system.time(res2 <- calc_genoprob_R(m = 4,
                                    mrknames = names(H),
                                    indnames = names(H[[1]]),
                                    n.mrk = n.mrk,
                                    n.ind = n.ind,
                                    emit = log_E,
                                    rf_vec = map.err[[2]]))
image(res2$probs[,,5])
res2$map
res3 <- calc_homoprob(res2)
head(res3$homoprob)
head(pedigree)
u <- pedigree$pop
names(u) <- pedigree$Name
v<-u[as.character(res3$homoprob$individual)]
x <- paste(as.character(res3$homoprob$homolog), v, sep="")
y <- sort(unique(x))
z <- matrix(c("a1", "a2", "a3", "a4", "b1", "b2", "b3", "b4",
               "a1", "a2", "a3", "a4", "c1", "c2", "c3", "c4",
               "b1", "b2", "b3", "b4", "c1", "c2", "c3", "c4"), 3, 8, byrow = T)
z <- as.vector(z)
names(z) <- y
res3$homoprob$homolog <- z[x]
head(res3$homoprob)
par(ask = T)
for(i in 1:420){
  plot(res3, ind = i, use.plotly = FALSE)
}

plot(res3, ind = 3, use.plotly = FALSE) ## axb
plot(res3, ind = 64, use.plotly = FALSE) ## bxc
plot(res3, ind = 132, use.plotly = FALSE) ## axc

save.image("~/repos/highprecHMM/potato/result_ch6.rda")


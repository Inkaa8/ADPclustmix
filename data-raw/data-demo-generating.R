# Demo dataset 1
set.seed(123)
new_dat1 <- genComplexMixDat(n = 100, K = 2, p.con = 2, p.cat = 2, MaxOmega = 0,
                             p.noise.con = 0, p.noise.cat = 0, n.out = 0,
                             sph=FALSE, hom=FALSE)
dat_demo1 <- new_dat1$SimDat

usethis::use_data(dat_demo1, compress = "xz", overwrite = TRUE)


# Demo dataset 2
set.seed(666)
new_dat2 <- genComplexMixDat_bc(n = 500, K = 4, p.con = 5, p.cat = 5, MaxOmega = 0,
                                p.noise.con = 1, p.noise.cat = 1, lambda = rep(2,12),
                                sph=FALSE, hom=FALSE)
dat_demo2 <- new_dat2$SimDat

usethis::use_data(dat_demo2, compress = "xz", overwrite = TRUE)
